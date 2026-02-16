#Installing packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(pheatmap)

#Imaging quality control - RCC files
files <- list.files("C:/ST", pattern = "\\.RCC$", full.names = TRUE)
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
qc_results <- do.call(rbind, lapply(files, ImagingQC))
print(qc_results)
write.csv(qc_results, "RCC_QC_report.csv", row.names = FALSE)

#Preparing data
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
annotation = read.table("C:/ST/annotation.txt", sep = "\t", header = TRUE, colClasses="character")
rn2 = annotation[,c(1:3)]
rn2 = apply(rn2,1,function(x) paste(x[1:3], collapse = ("_")))
annotation = annotation[,-c(1:3,5:8)]
rownames(annotation) = rn2

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

#Limit of detection berechnen
neg_control = matrix_num[c("HYB-NEG"), ]
neg_control = as.numeric(neg_control)
neg_control[neg_control == 0] <- 0.1
pos_control = matrix_num[c("HYB-POS"), ]
pos_control = as.numeric(pos_control)
mean_neg = mean(neg_control)
sd_neg = sd(neg_control)
LOD_threshold = mean_neg + 2 * sd_neg
LOD_pass = pos_control > LOD_threshold
prop_pass = mean(LOD_pass)
geo_mean_neg = exp(mean(log(neg_control)))

#Assay efficiency berechnen
reference = matrix_num[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), ]
background = matrix_num[c("NegPrb1","NegPrb2","NegPrb3","NegPrb4","NegPrb5","NegPrb6"), ]


geo_mean_pos = exp(mean(log(pos_control)))

all_values = as.numeric(matrix_num)
all_values[all_values <= 0] <- 0.1
geo_mean_all = exp(mean(log(all_values)))
assay_efficiency = geo_mean_pos/geo_mean_all



#Background and Normalization


#Ratios and Differential Expression

