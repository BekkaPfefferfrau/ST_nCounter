### Spatial transcriptomic data analysis - GeoMX DSP nCounter readout
## The follwoing code is written for data that was aquired using the Nanostring GeoMX panel "Immune pathways" and read out on the MAX/FLEX nCounter
# The analysis steps are based on the manual "MAN-10154-01, GeoMx DSP Data Analysis User Manual, section Data QC for nCounter Readout".


##Installing packages
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
library(pheatmap)


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
