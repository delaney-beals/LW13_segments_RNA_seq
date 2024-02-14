library(tidyverse)
library(ggplot2)
library(gplots)
library(pheatmap)

# making a heatmap of log2 fold change of all genes in line vs above and line vs below comparisons
segments = read.csv("data/fold_change_heatmap.csv", head = T)
data <- data.frame(segments)
datas <-data %>% remove_rownames %>% column_to_rownames(var="gene.ID")

pale <- colorRampPalette(c("navy","white","firebrick3"))

pheatmap(datas,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         color = pale(50))

#end