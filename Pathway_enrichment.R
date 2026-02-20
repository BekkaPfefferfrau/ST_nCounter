library(openxlsx)
library(pheatmap)

# Base directory and normalization methods
base_dir <- "C:/ST/reanalisincounter"
norm_methods <- c("backgorund corr", "HK", "SCALINGarea", "SCALINGnuclei")

# Annotations
annotations <- c(
  "distance CD3 near vs baseline far", 
  "distance CD30 near vs baseline far", 
  "distance CD68 near vs baseline far", 
  "histology CD3 unk vs baseline common", 
  "histology CD30 unk vs baseline common", 
  "histology CD68 unk vs baseline common", 
  "relapse CD3 relapse vs baseline CR", 
  "relapse CD30 relapse vs baseline CR", 
  "relapse CD68 relapse vs baseline CR"
)

heatmap_list <- list()

for(annot in annotations) {
  
  annot_data <- list()
  row_labels <- NULL
  
  for(norm in norm_methods) {
    # Match files starting with normalization method and containing annotation
    file_path <- list.files(
      base_dir,
      pattern = paste0("^", norm, ".*", annot, ".*\\.xlsx$"),
      full.names = TRUE
    )
    
    if(length(file_path) == 0) {
      warning(paste("No file found for", annot, "with normalization", norm))
      annot_data[[norm]] <- NA
      next
    }
    
    temp <- read.xlsx(file_path)
    
    if(ncol(temp) < 4 || nrow(temp) <= 7) {
      warning(paste("File too small or missing required columns in", file_path))
      annot_data[[norm]] <- NA
      next
    }
    
    # Column 4 = Log2 values, Column 3 = gene names
    annot_data[[norm]] <- as.numeric(temp[-c(1:7), 4])
    
    # Store gene names from first valid file
    if(is.null(row_labels)) {
      row_labels <- temp[-c(1:7), 3]
    }
  }
  
  # Pad all columns to same length
  max_len <- max(sapply(annot_data, length))
  for(i in seq_along(annot_data)) {
    length(annot_data[[i]]) <- max_len
  }
  
  # Pad row names to match max_len
  if(length(row_labels) < max_len){
    row_labels <- c(row_labels, rep("", max_len - length(row_labels)))
  }
  
  annot_data_df <- as.data.frame(annot_data)
  colnames(annot_data_df) <- paste0(norm_methods[1:ncol(annot_data_df)], "_Log2")
  rownames(annot_data_df) <- row_labels
  
  heatmap_list[[annot]] <- annot_data_df
}

# Save all heatmaps to PDF
pdf(file = "normalisation_comparison.pdf", width = 12, height = 10)

for(annot in names(heatmap_list)) {
  mat <- as.matrix(heatmap_list[[annot]])
  
  # Define symmetric color palette around 0
  max_val <- max(abs(mat), na.rm = TRUE)  # largest absolute value
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  print(
    pheatmap(
      mat,
      main = annot,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      color = color_palette,
      breaks = seq(-max_val, max_val, length.out = 101),  # center at 0
      na_col = "grey",
      show_rownames = TRUE,
      fontsize_row = 6
    )
  )
}

dev.off()