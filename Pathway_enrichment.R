library(pheatmap)
library(readxl)

files <- list.files(
  "C:/ST/reanalisincounter",
  pattern = "^[^~].*\\.xlsx$",
  full.names = TRUE
)

for (file in files) {

  raw_data <- read_excel(file, col_names = FALSE)

  data_matrix <- as.matrix(apply(raw_data[-c(1:7), -c(1:3)], 2, as.numeric))
  

  colnames(data_matrix) <- as.character(unlist(raw_data[7, -c(1:3)]))

  rownames(data_matrix) <- as.character(unlist(raw_data[-c(1:7), 2]))

  data_matrix[!is.finite(data_matrix)] <- 0
  
  pheatmap(
    data_matrix,
    scale = "none",
    color = colorRampPalette(c("white", "blue"))(100),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    fontsize_row = 6,
    fontsize_col = 6,
    main = basename(file)
  )
}
