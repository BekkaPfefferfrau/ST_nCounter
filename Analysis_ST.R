### Spatial transcriptomic data analysis
## The follwoing code is written for data that was aquired using the Nanostring GeoMX panel "Immune pathways" and read out on the MAX/FLEX nCounter
# The analysis steps are based on the manual "MAN-10154-01, GeoMx DSP Data Analysis User Manual, section Data QC for nCounter Readout".


##Installing packages
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
write.csv(QC_imaging_binding, "RCC_QC_report.csv", row.names = FALSE)

##Preparing raw data matrix .txt
# **this untangled data was provided by our collaborator, using DSPDA suite**
matrix1 = read.table("C:/ST/matrix.txt", sep = "\t", header = FALSE)
matrix2 = matrix1[-c(1:7),-c(1:2)]
matrix2 = apply(matrix2, 2, as.numeric)
matrix2 = as.matrix(matrix2)
cn = matrix1[c(1,2,5),-c(1:2)]
cn = t(cn)
cn = apply(cn,1,function(x) paste(x[1:3], collapse = ("_")))
rn = matrix1[-c(1:7),c(2)]
rownames(matrix2) = rn
colnames(matrix2) = cn


## Positive control normalization QC
pos_control = matrix2[c("HYB-POS"), ]
pos_ctrl_norm_factor <- function(pos_control) {
  ref_median <- median(pos_control[pos_control > 0])
  factor <- ifelse(pos_control == 0, NA, ref_median / pos_control)
  qc <- ifelse(!is.na(factor) & factor >= 0.3 & factor <= 3, "PASS", "FAIL")
  return(list(Factor = factor, QC = qc))
}
result <- pos_ctrl_norm_factor(pos_control)
df <- data.frame(Raw = pos_control, Factor = result$Factor, QC = result$QC)
df$ROI <- names(pos_control)
write.csv(df, "HYB-POS_QC_report.csv", row.names = FALSE)


## Positive control data normalization
matrix3 <- matrix2[!rownames(matrix2) %in% "HYB-POS", ]
factor_vec <- result$Factor[colnames(matrix3)]
matrix3 <- sweep(matrix3, 2, factor_vec, `*`)
matrix3 <- matrix3[, result$QC == "PASS"]


#Housekeeping genes normalization
HK_control = matrix2[c("OAZ1", "POLR2A", "RAB7A", "SDHA", "UBB"), ]
HK_matrix <- matrix3[HK_control, , drop = FALSE]
HK_sum_per_roi <- colSums(HK_matrix)
HK_median <- median(HK_sum_per_roi)
HK_factor <- HK_median / HK_sum_per_roi
HK_QC <- ifelse(HK_factor >= 0.1 & HK_factor <= 10, "PASS", "FAIL")
hk_qc <- data.frame(Raw = HK_control, HK_QC)
write.csv(hk_qc, "HK_QC_report.csv", row.names = FALSE)
matrix_hk_norm <- sweep(matrix3, 2, hk_factor, `*`)


# Negative probes
neg_control = matrix2[c("NegPrb1", "NegPrb2", "NegPrb3", "NegPrb4", "NegPrb5"), ]

# Annotation
annotation1 = read.table("C:/ST/annotation.txt", sep = "\t", header = TRUE, colClasses="character")
rn2 = annotation1[,c(1:3)]
rn2 = apply(rn2,1,function(x) paste(x[1:3], collapse = ("_")))
annotation1 = annotation1[,-c(1:3,5:8)]
rownames(annotation1) = rn2
common_rois <- intersect(colnames(matrix3), rownames(annotation1))
annotation2 <- annotation1[common_rois, ]

#Visualizing first part just for fun
sorted_cols <- rownames(annotation2)[order(annotation2$distance)]
matrix3 <- matrix3[, sorted_cols]
pheatmap(matrix3,
         scale = "none", color = colorRampPalette(c("blue","white","green"))(100),
         breaks = seq(min(matrix3),
                      max(matrix3), length.out = 101),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         width = 15,
         height = 10)
