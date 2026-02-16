### Spatial transcriptomic data analysis
## The follwoing code is written for data that was aquired using the Nanostring GeoMX panel "Immune pathways" and read out on the MAX/FLEX nCounter
# The analysis steps are based on the manual "MAN-10154-01, GeoMx DSP Data Analysis User Manual, section Data QC for nCounter Readout".

##Installing packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
library(pheatmap)

## In our case each .RCC files corresponds to 8 wells, resulting in a total of 24 .RCC files for 2x96-well plates
#Fields of view (FOV) quality control (QC) and binding density (BD) using QC .RCC files

files <- list.files("C:/ST", pattern = "^[^~].*\\.RCC$", full.names = TRUE)
ImagingQC <- function(file_path) {
    lines <- readLines(file_path)
    fov_line <- lines[grep("^FovCount,", lines)]
    fov_count <- as.numeric(sub("FovCount,", "", fov_line))
    fovc_line <- lines[grep("^FovCounted,", lines)]
    fov_counted <- as.numeric(sub("FovCounted,", "", fovc_line))
    
    if(length(fov_line) == 0 | length(fovc_line) == 0){
      imaging_qc <- "Missing FOV info"
    } else {
      fov_ratio <- fov_counted / fov_count
      imaging_qc <- ifelse(fov_ratio >= 0.75, "PASS", "FAIL")
    }
    bd_line <- lines[grep("^BindingDensity,", lines)]
    if(length(bd_line) == 0){
      binding_qc <- "Missing BindingDensity"
    } else {
      binding_density <- as.numeric(sub("BindingDensity,", "", bd_line))
      binding_qc <- ifelse(binding_density >= 0.1 & binding_density <= 2.25, "PASS", "FAIL")
    }
    
    return(data.frame(
      File = basename(file_path),
      ImagingQC = imaging_qc,
      BindingDensityQC = binding_qc,
      stringsAsFactors = FALSE
    ))
}

QC_imaging_binding <- do.call(rbind, lapply(files, ImagingQC))
print(QC_imaging_binding)
write.csv(QC_imaging_binding, "RCC_QC_report.csv", row.names = FALSE)

##Preparing raw data matrix .txt
# **this untangled data was provided by our collaborator, using DSPDA suite**
matrix = read.table("C:/ST/matrix.txt", sep = "\t", header = FALSE)
matrix_num = matrix[-c(1:7),-c(1:2)]
matrix_num = apply(matrix_num, 2, as.numeric)
matrix_num = as.matrix(matrix_num)
cn = matrix[c(1,2,5),-c(1:2)]
cn = t(cn)
cn = apply(cn,1,function(x) paste(x[1:3], collapse = ("_")))
rn = matrix[-c(1:7),c(2)]
rownames(matrix_num) = rn
colnames(matrix_num) = cn
neg_control = matrix_num[c("HYB-NEG"), ]
neg_control = as.numeric(neg_control)
neg_control[neg_control == 0] <- 0.1
pos_control = matrix_num[c("HYB-POS"), ]

## Positive control normalization QC
pos_ctrl_norm_factor <- function(pos_control) {
  ref_median <- median(pos_control[pos_control > 0])
  factor <- ifelse(pos_control == 0, NA, ref_median / pos_control)
  qc <- ifelse(!is.na(factor) & factor >= 0.3 & factor <= 3, "PASS", "FAIL")
  return(list(Factor = factor, QC = qc))
}
result <- pos_ctrl_norm_factor(pos_control)
df <- data.frame(Raw = pos_control, Factor = result$Factor, QC = result$QC)
print(df)
write.csv(df, "Positive_control_norm_factor.csv", row.names = FALSE)

## Positive control data normalization
matrix_norm <- matrix_num[-86, ]
factor_vec <- result$Factor[colnames(matrix_norm)]
matrix_norm <- sweep(matrix_norm, 2, factor_vec, `*`)
matrix_norm <- matrix_norm[, result$QC == "PASS"]
pheatmap(matrix_norm,
         scale = "none", color = colorRampPalette(c("blue","white","green"))(100),
         breaks = seq(min(matrix_sorted),
                      max(matrix_sorted), length.out = 101),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         width = 15,
         height = 10)

##Limit of detection (LOD) control
mean_neg = mean(neg_control)
sd_neg = sd(neg_control)
LOD_threshold = mean_neg + 2 * sd_neg
LOD_pass = pos_control > LOD_threshold
prop_pass = mean(LOD_pass)
geo_mean_neg = exp(mean(log(neg_control)))

LOD_list <- data.frame(
  ID = IDs,
  LOD_QC = ifelse(LOD_pass, "PASS", "FAIL"),
  stringsAsFactors = FALSE
)

print(LOD_list)

print(LOD_list)


#Assay efficiency
reference = matrix_num[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), ]
background = matrix_num[c("NegPrb1","NegPrb2","NegPrb3","NegPrb4","NegPrb5","NegPrb6"), ]


geo_mean_pos = exp(mean(log(pos_control)))

all_values = as.numeric(matrix_num)
all_values[all_values <= 0] <- 0.1
geo_mean_all = exp(mean(log(all_values)))
assay_efficiency = geo_mean_pos/geo_mean_all




#Visualising first part just for fun
sorted_cols <- rownames(annotation)[order(annotation$Segment.Tags)]
matrix_sorted <- matrix_num[, sorted_cols]

pheatmap(matrix_sorted,
         scale = "none", color = colorRampPalette(c("blue","white","green"))(100),
         breaks = seq(min(matrix_sorted),
                      max(matrix_sorted), length.out = 101),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         width = 15,
         height = 10)


#Background and Normalization


#Ratios and Differential Expression

# Annotation
annotation = read.table("C:/ST/annotation.txt", sep = "\t", header = TRUE, colClasses="character")
rn2 = annotation[,c(1:3)]
rn2 = apply(rn2,1,function(x) paste(x[1:3], collapse = ("_")))
annotation = annotation[,-c(1:3,5:8)]
rownames(annotation) = rn2