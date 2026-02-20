## In our case each .RCC files corresponds to 8 wells, resulting in a total of 24 .RCC files for 2x96-well plates
#Fields of view (FOV) quality control (QC) and binding density (BD) using QC .RCC files

QC_imaging_binding <- (function() {
  
  files <- list.files(
    "C:/Users/Linfomi2/Documents/ST_nCounter/.RCC",
    pattern = "\\.RCC$",
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
      .RCC_file = basename(file_path),
      FoVRaw = fov_counted,
      FoVRegistrationQC = imaging_qc,
      BindingDensityRaw = bd,
      BindingDensityQC = binding_qc,
      stringsAsFactors = FALSE
    )
  }
  
  qc_df <- do.call(rbind, lapply(files, ImagingQC))
  
  
  write.csv(
    qc_df,
    "C:/Users/Linfomi2/Documents/ST_nCounter/Reports/1.RCC_QC_report.csv",
    row.names = FALSE
  )
  
  qc_df
})()


##Preparing raw data matrix .txt
# **this untangled data was provided by our collaborator, using DSPDA suite**

matrix1 = read.table("C:/Users/Linfomi2/Documents/ST_nCounter/Normalisation_R/matrix.txt", sep = "\t", header = FALSE)
#creating "matrix1" from the raw data file

matrix2 <- (function(x){
  m <- as.matrix(apply(x[-c(1:7), -c(1:2)], 2, as.numeric))
  rownames(m) <- x[-c(1:7), 2]
  colnames(m) <- apply(t(x[c(1, 2, 5), -c(1:2)]), 1, function(y) paste(y[1:3], collapse = "_"))
  return(m)
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
  "C:/Users/Linfomi2/Documents/ST_nCounter/Reports/2.HYB-POS_QC_report.csv",
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