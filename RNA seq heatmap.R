setwd("/Users/Delaney/OneDrive/unknown/Documents/University of Utah/Puri Lab/1. Gradient syringe validation/~RNA-seq/LW13 segment")
getwd()


library(tidyverse)
library(ggplot2)
library(gplots)
library(pheatmap)



segments = read.csv("fold_change_heatmap.csv", head = T)
data <- data.frame(segments)
datas <-data %>% remove_rownames %>% column_to_rownames(var="gene.ID")
dataz <- as.matrix(data)


t<- t(dataz)
heatmap(as.matrix(datas))

heatmap(t)

heatmap(t, scale="column")

#heatmap.2 package
heatmap.2(as.matrix(datas), 
          trace="none",
          main="log change",
          cluster_rows = FALSE,
          dendrogram = "none",
          margins = c(5,20))

heatmap.2(as.matrix(datas[, -1]),
          trace="none", 
          main="FPKM of segments",
          #ColSideColors=col.cell,
          density.info="none",
          Colv=FALSE,
          Rowv=FALSE,
          cexRow = TRUE,
          scale="row")

dev.off()

length(unique(data[,1]))==length(data[,1])
rownames(data)=paste(c(1:nrow(data)),data[,1])
data$Gene.Function <- NULL
pale <- colorRampPalette(c("navy","white","firebrick3"))

pheatmap(datas,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         color = pale(50))
pheatmap(datas)
# Create data matrix


# Plot a heatmap 
heatmap(t,Rowv=NA,Colv=NA, scale="column")

# Plot a corresponding legend
legend(x="right", legend=c("min", "med", "max"),fill=heat.colors(3))

ggplot(t, aes(X,Y, fill=Z) + geom_title())

ggplot(data, aes(x)) +geom_tile()


#PCA plot
library(stats)
library(ggplot2)
library(factoextra)
library(ggfortify)
library(devtools)
install.packages("ggfortify")


df <- data.frame(read.csv("segments_FPKM_expressionmatrix.csv", head = T))
prcomp(df[,1])
pcaResults <- prcomp(df)
autoplot(prcomp(pca))
print(pcaResults)
ggplot()

matrix = read.csv("segments_FPKM_expressionmatrix.csv", head = T)
df <- data.frame(matrix)
dfs <- df %>% remove_rownames %>% column_to_rownames(var="gene.ID")
pca <- as.matrix(dfs)
autoplot(prcomp(pca))
pcares <- prcomp(pca)


pca_plot <- pca$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table 
  full_join(sample_info, by = "sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(minute), shape = strain)) +
  geom_point()

pcaData <- as.data.frame(pcares$x[, 1:2]) # extract first two PCs
 # add species to df
colnames(pcaData) <- c("PC1", "PC2") # change column names
