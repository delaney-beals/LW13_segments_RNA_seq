##Attempt at volcano plot for Line vs. Below

# remove rows that contain NA values
lb <- diff.exp.matrix.line.trimmed.VS.below.trimmed[complete.cases(diff.exp.matrix.line.trimmed.VS.below.trimmed), ]


ggplot(data=lb, aes(x=log2_fold_change, y=-log10(p_value))) + 
  geom_point() + theme_minimal() + 
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.001), col="red")

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
lb$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.001, set as "UP" 
lb$diffexpressed[lb$log2_fold_change > 1 & lb$p_value < 0.001] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.001, set as "DOWN"
lb$diffexpressed[lb$log2_fold_change < -1 & lb$p_value < 0.001] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
ggplot(data=lb, aes(x=log2_fold_change, y=-log10(p_value), col=diffexpressed)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.001), col="red")


#add gene function labels instead of 292383843_gene names
install.packages("tidyverse")
library(tidyverse)
library(rvest)
library(dplyr)
library(ggrepel)

write.csv(lb,"C:/Users/Delaney/OneDrive/unknown/Documents/University of Utah/Puri Lab/1. Gradient syringe validation/~RNA-seq/diff-exp-matrix-line-trimmed-VS-below-trimmed.TSV/diff_exp_mat_line_below_names.csv", row.names = TRUE)

df3 = read.csv('diff_exp_mat_line_below_names.csv', head = T)
df4 = read.csv('line-below-features-all.csv', head = T, )

lb = df3 %>% left_join(df4, by=c("gene.ID"))

# option to save melted table as separate file
write.csv(lb,"C:/Users/Delaney/OneDrive/unknown/Documents/University of Utah/Puri Lab/1. Gradient syringe validation/~RNA-seq/LW13 segment/line-below-diff-exp-summary.csv", row.names = TRUE)


# add a column of NAs
lb$sigexpress <- "NO"
# if log2Foldchange > 1 and pvalue < 0.001, set as "UP" 
lb$sigexpress[lb$log2_fold_change > 1 & lb$p_value < 0.00001] <- "high"
# if log2Foldchange < -1 and pvalue < 0.001, set as "DOWN"
lb$sigexpress[lb$log2_fold_change < -1 & lb$p_value < 0.00001] <- "low"

lb$relabel <- NA
lb$relabel[lb$sigexpress != "NO"] <- lb$Function[lb$sigexpress != "NO"]

# final plot

ggplot(data=lb, aes(x=log2_fold_change, y=-log10(p_value), label=relabel)) + 
  geom_point(aes(colour = diffexpressed), size = 1.5) + 
  scale_x_continuous(limits = c(-5,5), breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5)) +  
  scale_colour_manual(values=c("#0000FF","#FF0000","#808080"), name="diffexpressed", breaks=c("DOWN","UP","NO")) + 
  geom_text_repel(size = 3) +  
  ggtitle("Differential Expression of Line vs. Below") +  
  geom_vline(xintercept=c(-1, 1), col="red") +  
  geom_hline(yintercept=-log10(0.001), col="red") + 
  theme_classic()


final.lb <- ggplot(data=lb, aes(x=log2_fold_change, y=-log10(p_value), label=relabel)) + 
  geom_point(aes(colour = diffexpressed), size = 1.5) + 
  scale_x_continuous(limits = c(-5,5), breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5)) +  
  scale_colour_manual(values=c("#0000FF","#FF0000","#808080"), name="diffexpressed", breaks=c("DOWN","UP","NO")) + 
  geom_text_repel(size = 3) +  
  ggtitle("Differential Expression of Line vs. Below") +  
  geom_vline(xintercept=c(-1, 1), col="red") +  
  geom_hline(yintercept=-log10(0.001), col="red") + 
  theme_classic()


print(final.lb)
