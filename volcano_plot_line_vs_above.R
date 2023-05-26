##Attempt at volcano plot for Line vs. Above


# remove rows that contain NA values
la <- diff.exp.matrix.line.trimmed.VS.above.trimmed[complete.cases(diff.exp.matrix.line.trimmed.VS.above.trimmed), ]

ggplot(data=la, aes(x=log2_fold_change, y=-log10(p_value))) + 
  geom_point() + 
  theme_minimal() + 
  geom_vline(xintercept=c(-1, 1), col="red") +  
  geom_hline(yintercept=-log10(0.001), col="red")


# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
la$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.001, set as "UP" 
la$diffexpressed[la$log2_fold_change > 1 & la$p_value < 0.001] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.001, set as "DOWN"
la$diffexpressed[la$log2_fold_change < -1 & la$p_value < 0.001] <- "DOWN"


# Re-plot but this time color the points with "diffexpressed"
ggplot(data=la, aes(x=log2_fold_change, y=-log10(p_value), col=diffexpressed)) + 
  geom_point() + theme_minimal() + 
  geom_vline(xintercept=c(-1, 1), col="red") +  
  geom_hline(yintercept=-log10(0.001), col="red")


# Now write down the name of genes beside the points...

#add gene function labels instead of 292383843_gene names
install.packages("tidyverse")
library(tidyverse)
library(rvest)
library(dplyr)
library(ggrepel)

# need to melt our two tables together so that the gene.ID is connected to the significance values and the functions assigned to each gene
write.csv(la,"C:/Users/Delaney/OneDrive/unknown/Documents/University of Utah/Puri Lab/1. Gradient syringe validation/~RNA-seq/diff-exp-matrix-line-trimmed-VS-above-trimmed.TSV/diff_exp_mat_line_above_names.csv", row.names = TRUE)

df1 = read.csv('diff_exp_mat_line_above_names.csv', head = T)
df2 = read.csv('line-above-features-all.csv', head = T, )

la = df1 %>% left_join(df2, by=c("gene.ID"))

write.csv(la,"C:/Users/Delaney/OneDrive/unknown/Documents/University of Utah/Puri Lab/1. Gradient syringe validation/~RNA-seq/LW13 segment/line-above-diff-exp-summary.csv", row.names = TRUE)

#start here for making plots again
la = read.csv("line-above-diff-exp-summary.csv", head = T)
# make a new column that will be based on an even smaller selection of the significant/colored dots; these ones will be labeled with their functions
# add a column of NAs
la$sigexpress <- "NO"
# if log2Foldchange > 1 and pvalue < 0.001, set as "UP" 
la$sigexpress[la$log2_fold_change > 1 & la$p_value < 0.000000000001] <- "high"
# if log2Foldchange < -1 and pvalue < 0.001, set as "DOWN"
la$sigexpress[la$log2_fold_change < -1 & la$p_value < 0.000000000001] <- "low"

la$zelabel <- NA
la$zelabel[la$sigexpress != "NO"] <- la$Function[la$sigexpress != "NO"]



# final plot

ggplot(data=la, aes(x=log2_fold_change, y=-log10(p_value), label=zelabel)) + 
  geom_point(aes(colour = diffexpressed), size = 1.5) + 
  scale_x_continuous(limits = c(-5,5), breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5)) +  
  scale_colour_manual(values=c("#0000FF","#FF0000","#808080"), name="diffexpressed", breaks=c("DOWN","UP","NO")) + 
  geom_text_repel(size = 3) +  
  ggtitle("Differential Expression of Line vs. Above") +  
  geom_vline(xintercept=c(-1, 1), col="red") +  
  geom_hline(yintercept=-log10(0.001), col="red") + 
  theme_classic()


final.la <- ggplot(data=la, aes(x=log2_fold_change, y=-log10(p_value), label=zelabel)) + 
  geom_point(aes(colour = diffexpressed), size = 1.5) + 
  scale_x_continuous(limits = c(-5,5), breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5)) +  
  scale_colour_manual(values=c("#0000FF","#FF0000","#808080"), name="diffexpressed", breaks=c("DOWN","UP","NO")) + 
  geom_text_repel(size = 3) +  
  ggtitle("Differential Expression of Line vs. Above") +  
  geom_vline(xintercept=c(-1, 1), col="red") +  
  geom_hline(yintercept=-log10(0.001), col="red") + 
  theme_classic()

print(final.la)

## trying to get diff fonts; none of this works
install.packages("extrafont")
library(extrafont)
font_import("Avenir Next LT Pro Regular")


library(extrafont)
library(remotes)
remotes::install_version("Rttf2pt1", version = "1.3.8")
extrafont::font_import()
font_import(paths = c("C:/Users/Delaney/OneDrive/unknown/Documents/University of Utah/Puri Lab/1. Gradient syringe validation/~RNA-seq/LW13 segment", prompt = F))

library(ggplot2)
qplot(1:10)+theme(text=element_text(family="Trebuchet MS"))


ggplot(data=la, aes(x=log2_fold_change, y=-log10(p_value), label=zelabel)) + 
  geom_point(aes(colour = diffexpressed), size = 1.5) + 
  scale_x_continuous(limits = c(-5,5), breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5)) +  
  scale_colour_manual(values=c("#0000FF","#FF0000","#808080"), name="diffexpressed", breaks=c("DOWN","UP","NO")) + 
  geom_text_repel(size = 3) +  
  ggtitle("Differential Expression of Line vs. Above") +  
  geom_vline(xintercept=c(-1, 1), col="red") +  
  geom_hline(yintercept=-log10(0.001), col="red") + 
  theme_classic()+
  theme(text = element_text(family = "Times New Roman"))


ggplot(data=la, aes(x=log2_fold_change, y=-log10(p_value))) + 
  geom_point(aes(colour = diffexpressed), size = 1.5) + 
  scale_x_continuous(limits = c(-5,5), breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5)) +  
  scale_colour_manual(values=c("#0000FF","#FF0000","#808080"), name="diffexpressed", breaks=c("DOWN","UP","NO")) + 
  ggtitle("Differential Expression of Line vs. Above") +  
  theme_classic()


#trying to color mutated genes by a diff color (https://www.biostars.org/p/437487/)

ggplot(data=la, aes(x=log2_fold_change, y=-log10(p_value), label=belabel)) + 
  geom_point(aes(colour = colorz), size = 1.5) + 
  scale_x_continuous(limits = c(-5,5), breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5)) +   
  ggtitle("Differential Expression of Line vs. Above") +  
  geom_vline(xintercept=c(-1, 1), col="red") +  
  geom_hline(yintercept=-log10(0.001), col="red") + 
  geom_text_repel(size = 3) +
  theme_classic()

la$interesting <- "NO"
la$interesting[la$colorz > 0] <- "yes"

la$belabel <- NA
la$belabel[la$interesting != "NO"] <- la$Function[la$interesting != "NO"]
