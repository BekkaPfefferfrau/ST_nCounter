library(openxlsx)
library(dplyr)
library(tidyr)
library(pheatmap)

base_dir <- "C:/ST/reanalisincounter"
norm_methods <- c("backgorund corr", "HK", "SCALINGarea", "SCALINGnuclei")

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

pdf(file = "normalisation_comparison_4heatmaps.pdf", width = 12, height = 10)

# Loop over normalization methods
for(norm in norm_methods){
  
  all_annotations_df <- list()
  
  for(annot in annotations){
    file_path <- list.files(
      base_dir,
      pattern = paste0("^", norm, ".*", annot, ".*\\.xlsx$"),
      full.names = TRUE
    )
    
    if(length(file_path) == 0){
      warning(paste("No file found for", annot, "with normalization", norm))
      next
    }
    
    temp <- read.xlsx(file_path)
    if(ncol(temp) < 4 || nrow(temp) <= 7){
      warning(paste("File too small or missing required columns in", file_path))
      next
    }
    
    log2_values <- as.numeric(temp[-c(1:7), 4])
    pathways <- temp[-c(1:7), 2]
    
    df <- data.frame(Pathways = pathways,
                     Value = log2_values,
                     Annotation = annot,
                     stringsAsFactors = FALSE) %>%
      mutate(Pathways = strsplit(as.character(Pathways), ",\\s*")) %>%
      tidyr::unnest(Pathways)
    
    all_annotations_df[[annot]] <- df
  }
  
  combined_df <- bind_rows(all_annotations_df)
  
  # Compute mean Log2 per pathway x annotation
  mat <- combined_df %>%
    group_by(Pathways, Annotation) %>%
    summarize(MeanLog2 = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Annotation, values_from = MeanLog2)
  
  rownames_mat <- mat$Pathways
  mat <- as.matrix(mat[,-1])
  rownames(mat) <- rownames_mat
  
  # Symmetric color scale
  max_val <- max(abs(mat), na.rm = TRUE)
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  pheatmap(mat,
           main = paste("Normalization:", norm),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           color = color_palette,
           breaks = seq(-max_val, max_val, length.out = 101),
           na_col = "grey",
           fontsize_row = 6,
           fontsize_col = 8)
}

dev.off()