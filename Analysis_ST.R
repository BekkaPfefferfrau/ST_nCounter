### Spatial transcriptomic data analysis - GeoMX DSP nCounter readout
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
#creating object "matrix1" from the raw data file

matrix2 = as.matrix(apply(matrix1[-c(1:7),-c(1:2)], 2, as.numeric))
#creating a numeric table from matrix1 that contains only the numbers

rownames(matrix2) = matrix1[-c(1:7),c(2)]
#setting the rownames, using the geneID string from matrix1

colnames(matrix2) = apply(t(matrix1[c(1,2,5),-c(1:2)]),1,function(x) paste(x[1:3], collapse = ("_")))
#setting the column names, using a string combination of sampleID_ROI_morphologymarker from matrix1


## Positive control QC

write.csv(
  data.frame(
    Sample_ID = names(matrix2["HYB-POS", ]),
    Pos_ctrl = matrix2["HYB-POS", ],
    Factor = (function(x) { med <- median(x[x > 0]); ifelse(x == 0, NA, med / x) })(matrix2["HYB-POS", ]),
    QC = (function(x) { 
      med <- median(x[x > 0]); 
      fac <- ifelse(x == 0, NA, med / x)
      ifelse(!is.na(fac) & fac >= 0.3 & fac <= 3, "PASS", "FAIL")
    })(matrix2["HYB-POS", ])
  ),
  "HYB-POS_QC_report.csv",
  row.names = FALSE
)
#creating a .csv file containing the HYB_POS_QC_report
#calculating the median of all positive controls, using only numerical larger than 0
#calculating the pos_ctrl_factor by dividing the pos_ctrl_median trough each value, if larger than 0; in case of 0 return NA
#testing if the pos_ctrl_factor meets the requirement lager or equal 0.3 and smaller or equal 3
#creating a list with the column name "Factor" for the pos_ctrl_factor and QC for the pos_ctrl_qc


## Positive control data normalization

matrix3 <- (function(x) {
  pos_ctrl <- x["HYB-POS", ]
  pos_ctrl_median <- median(pos_ctrl[pos_ctrl > 0])
  pos_ctrl_factor <- ifelse(pos_ctrl == 0, NA, pos_ctrl_median / pos_ctrl)
  pos_ctrl_qc <- ifelse(!is.na(pos_ctrl_factor) & pos_ctrl_factor >= 0.3 & pos_ctrl_factor <= 3, "PASS", "FAIL")
  cols <- intersect(colnames(x), names(pos_ctrl_factor))
  sweep(
    x[rownames(x) != "HYB-POS", cols, drop = FALSE],
    2,
    pos_ctrl_factor[cols],
    `*`
  )[, pos_ctrl_qc[cols] == "PASS", drop = FALSE]
})(matrix2)

#creating matrix3 by removing "HYB-POS" row from matrix2, making sure that the sampleIDs of matrix2 and the pos_ctrl_function results are matching
#Multiplying the positive_ctrl factor from the function pos_ctrl_function with the respective values from each column with the respective pos_ctrl_factor with sampleID matching
#removing all samples that failed the pos_ctrl_qc with sampleID matching


##Housekeeping QC

write.csv(
  data.frame(
    Sample_ID = colnames(matrix3),
    t(matrix3[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), , drop = FALSE]),
    HK_sum = colSums(matrix3[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), , drop = FALSE]),
    HK_median = median(colSums(matrix3[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), , drop = FALSE])),
    Factor = (function(x) { med <- median(x); med / x })(
      colSums(matrix3[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), , drop = FALSE])
    ),
    QC = (function(x) { fac <- median(x) / x; ifelse(fac >= 0.1 & fac <= 10, "PASS", "FAIL") })(
      colSums(matrix3[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), , drop = FALSE])
    )
  ),
  "HK_QC_report.csv",
  row.names = FALSE
)
#creating a .csv file with the HK_report
#Rownames are taken from matrix3


## Housekeeping data normalization (not necessary for our samples as they all passed but I will still write the code for the sake of completeness)








## Negative control QC
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
