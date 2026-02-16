#Installing packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(pheatmap)

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

