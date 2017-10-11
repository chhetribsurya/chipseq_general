
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)
library(gplots)
library(xlsx)


args <-  commandArgs(TRUE)
input_file <- args[1]
tf_name <- args[2]
out_dir_name <- args[3]

read_file <- fread("~/Dropbox/encode_3/jund_test/meth_analysis/final_wgbs_tf_ideas_intersect_heatmap_data.bed", sep="\t", header=TRUE)
read_file <- fread("~/Dropbox/encode_3/jund_test/meth_analysis/final_wgbs_tf_ideas_intersect_heatmap_data_for_HA_IDR_final.bed", sep="\t", header=TRUE)
read_file <- fread("~/Dropbox/encode_3/jund_test/meth_analysis/final_wgbs_tf_ideas_intersect_heatmap_data_for_ENC_IDR_final.bed", sep="\t", header=TRUE)

read_df <- as.data.frame(read_file)
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)

output_file_name <- paste0("~/Dropbox", "/", "JUND_methylation_contribution_from_different_regions.pdf")        
output_file_name <- paste0("~/Dropbox", "/", "JUND_methylation_contribution_from_different_regions_HA_IDR_final.pdf")        
output_file_name <- paste0("~/Dropbox", "/", "JUND_methylation_contribution_from_different_regions_ENC_IDR_final.pdf")        
pdf(output_file_name)

summary(read_df) #max_val = 325
par(cex.main=0.6)
heatmap.2(mat_data,
  main = "JUND methylation contribution from various states", 
  xlab = "IDEAS States",
  ylab = " ",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.45,
  cexCol = 0.65,
  cex.main=0.3,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    

dev.off()

