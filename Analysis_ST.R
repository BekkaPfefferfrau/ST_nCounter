### Spatial transcriptomic data analysis - GeoMX DSP nCounter readout
## The follwoing code is written for data that was aquired using the Nanostring GeoMX panel "Immune pathways" and read out on the MAX/FLEX nCounter
# The analysis steps are based on the manual "MAN-10154-01, GeoMx DSP Data Analysis User Manual, section Data QC for nCounter Readout".


##Installing packages
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
library(pheatmap)


## In our case each .RCC files corresponds to 8 wells, resulting in a total of 24 .RCC files for 2x96-well plates
#Fields of view (FOV) quality control (QC) and binding density (BD) using QC .RCC files

QC_imaging_binding <- (function() {

  files <- list.files(
    "C:/ST",
    pattern = "^[^~].*\\.RCC$",
    full.names = TRUE
  )

  ImagingQC <- function(file_path) {

    lines <- readLines(file_path)

    fov_line  <- lines[grep("^FovCount,", lines)]
    fovc_line <- lines[grep("^FovCounted,", lines)]

    if (length(fov_line) == 0 || length(fovc_line) == 0) {
      imaging_qc <- "Missing FOV info"
    } else {
      fov_count   <- as.numeric(sub("FovCount,", "", fov_line))
      fov_counted <- as.numeric(sub("FovCounted,", "", fovc_line))
      imaging_qc  <- ifelse(fov_counted / fov_count >= 0.75, "PASS", "FAIL")
    }

    bd_line <- lines[grep("^BindingDensity,", lines)]
    if (length(bd_line) == 0) {
      binding_qc <- "Missing BindingDensity"
    } else {
      bd <- as.numeric(sub("BindingDensity,", "", bd_line))
      binding_qc <- ifelse(bd >= 0.1 & bd <= 2.25, "PASS", "FAIL")
    }

    data.frame(
      File = basename(file_path),
      ImagingQC = imaging_qc,
      BindingDensityQC = binding_qc,
      stringsAsFactors = FALSE
    )
  }

  qc_df <- do.call(rbind, lapply(files, ImagingQC))

  write.csv(
    qc_df,
    "C:/ST/RCC_QC_report.csv",
    row.names = FALSE
  )

  qc_df
})()


##Preparing raw data matrix .txt
# **this untangled data was provided by our collaborator, using DSPDA suite**

matrix1 = read.table("C:/ST/matrix.txt", sep = "\t", header = FALSE)
#creating "matrix1" from the raw data file

matrix2 <- (function(x){
  m <- as.matrix(apply(x[-c(1:7), -c(1:2)], 2, as.numeric))
  rownames(m) <- x[-c(1:7), 2]
  colnames(m) <- apply(t(x[c(1, 2, 5), -c(1:2)]), 1, function(y) paste(y[1:3], collapse = "_"))
  return(m)  # <- entscheidend
})(matrix1)
#creating a numeric table from matrix1 that contains only the numbers
#setting the rownames, using the geneID string from matrix1
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

pos_ctrl_raw_boxplot <- (function(x){
  pos_data <- x["HYB-POS", , drop = FALSE]
  boxplot(
    pos_data,
    names = colnames(x),
    las = 2,
    main = "HYB-POS per sample",
    ylab = "Signal",
    outline = FALSE
  )
})(matrix2)

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

pos_ctrl_ratio_boxplot <- (function(x){
  pos_ctrl <- x["HYB-POS", , drop = FALSE]
  pos_ctrl_median <- median(pos_ctrl[pos_ctrl > 0])
  pos_ctrl_factor <- ifelse(pos_ctrl == 0, NA, pos_ctrl_median / pos_ctrl)
  boxplot(
    pos_ctrl_factor,
    names = colnames(x),
    las = 2,
    main = "HYB-POS factor per sample",
    ylab = "Signal",
    outline = FALSE
  )
  abline(h = 3, col = "red", lty = 1, lwd = 1)
  abline(h = 0.3, col = "red", lty = 1, lwd = 1)
})(matrix2)

##Scale to area
#here we should do the first normalization step according to the ROI are but I do not have this information
#--> need to ask collaborator to provide this


##Housekeeping QC

write.csv(
  data.frame(
    Sample_ID = colnames(matrix3),
    t(matrix3[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), , drop = FALSE]),
    Factor = (function(x) { 
      med <- median(x[x > 0])
      ifelse(x == 0, NA, med / x)
    })(colSums(matrix3[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), , drop = FALSE])),
    QC = (function(x) { 
      med <- median(x[x > 0])
      fac <- ifelse(x == 0, NA, med / x)
      ifelse(!is.na(fac) & fac >= 0.1 & fac <= 10, "PASS", "FAIL")
    })(colSums(matrix3[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), , drop = FALSE]))
  ),
  "HK_QC_report.csv",
  row.names = FALSE
)
#creating a .csv file with the HK_report
#Rownames are taken from matrix3


## Housekeeping data normalization

matrix4 <- (function(x) {
  HK_genes <- matrix3[c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), , drop = FALSE]
  HK_ctrl <- colSums(HK_genes)
  HK_ctrl_median <- median(HK_ctrl)
  HK_ctrl_factor <- HK_ctrl_median / HK_ctrl
  HK_ctrl_qc <- ifelse(HK_ctrl_factor >= 0.1 & HK_ctrl_factor <= 10, "PASS", "FAIL")
  cols <- intersect(colnames(x), names(HK_ctrl_factor))
  sweep(
    x[!rownames(x) %in% c("OAZ1","POLR2A","RAB7A","SDHA","UBB"), cols, drop = FALSE],
    2,
    HK_ctrl_factor[cols],
    `*`
  )[, HK_ctrl_qc[cols] == "PASS", drop = FALSE]
})(matrix3)


## Negative control data normalization

matrix5 <- (function(x) {
  neg_ctrl <- matrix4[c("NegPrb1", "NegPrb2", "NegPrb3", "NegPrb4", "NegPrb5"), ]
  neg_ctrl <- colSums(neg_ctrl)
  neg_ctrl_median <- median(neg_ctrl)
  neg_ctrl_factor <- neg_ctrl_median / neg_ctrl
  cols <- intersect(colnames(x), names(neg_ctrl_factor))
  sweep(
    x[!rownames(x) %in% c("NegPrb1", "NegPrb2", "NegPrb3", "NegPrb4", "NegPrb5"), cols, drop = FALSE],
    2,
    neg_ctrl_factor[cols],
    '*')
})(matrix4)
# calculating a negative control factor by using the 5 isotype controls
# then removing the isotype controls and multiplying each by its own factor

negative_ctrl_boxplot <- (function(x){
  neg_probes <- c("NegPrb1","NegPrb2","NegPrb3","NegPrb4","NegPrb5")
  neg_data <- matrix4[neg_probes, , drop = FALSE]
  boxplot(
    neg_data,
    las = 2,
    main = "Negative controls per sample",
    ylab = "Signal",
    outline = FALSE)
})(neg_probes)


## Negative hybridization control


# Annotation
annotation1 <- read.table("C:/ST/annotation.txt", sep = "\t", header = TRUE, colClasses="character")
annotation2 <- (function(x){
  rn <- apply(x[, 1:3], 1, paste, collapse = "_")
  rownames(x) <- rn
  x <- x[, -c(1:3, 5:8), drop = FALSE]
  cols <- intersect(rownames(x), colnames(matrix5))
  x <- x[cols, , drop = FALSE]
  return(x)
})(annotation1)


#Visualizing first part just for fun

matrix_sort_heatmap <- (function(x){
  sorted_cols <- rownames(annotation2)[order(annotation2$distance)]
  matrix <- matrix5[, sorted_cols]
  pheatmap(matrix,
         scale = "none",
         color = colorRampPalette(c("white","blue"))(100),
         breaks = seq(min(matrix5),
                      max(matrix5),
                      length.out = 100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         width = 15,
         height = 10)
  })(matrix5)
