library(ggpubr)
library(ggrepel)

# Load in table with differential expression information
diffexp <- read.csv("data/LW13_gradient_syringe_RNA_seq.csv")

# Function to create volcano plots, including optional labels
createVolcanoPlot <- function(data, x, y, colorLabel, ylims, xlims, yBreaks, xBreaks, colorTitle, labelColumn = NULL) {
  plot <- ggplot(data, aes(x = .data[[x]], y = -log10(.data[[y]]), color = .data[[colorLabel]])) +
    geom_point(size = 2) +
    scale_color_manual(values = c("#0571b0", "grey", "#ca0020"), 
                       labels = c("downregulated", "not significant", "upregulated")) +
    coord_cartesian(ylim = ylims, xlim = xlims) +
    labs(color = colorTitle, 
         x = expression(log[2] * " fold change"), 
         y = expression("-log"[10] * " p value")) +
    scale_x_continuous(breaks = seq(xlims[1], xlims[2], xBreaks)) +
    scale_y_continuous(breaks = seq(ylims[1], ylims[2], yBreaks)) +
    theme_classic(base_size = 15) +
    theme(
      axis.title.y = element_text(face = "bold", margin = margin(0, 2, 0, 0), size = rel(1.1), color = 'black'),
      axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(2, 0, 0, 0), size = rel(1.1), color = 'black'),
      plot.title = element_text(hjust = 0.5, family = "Arial"),
      axis.text = element_text(color = "black", size = rel(1.1))
    ) +
    geom_vline(xintercept = c(-1.5, 1.5), col = "black", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.00001), col = "black", linetype = 'dashed')
  
  if (!is.null(labelColumn)) {
    plot <- plot + geom_text_repel(aes(label = .data[[labelColumn]]), max.overlaps = Inf, size = 3)
  }
  
  return(plot)
}  

# Create plots without labels
line_below <- createVolcanoPlot(diffexp, "log2FC_band_vs_below", "pvalue_band_vs_below", "diffexp_band_vs_below", c(0, 15), c(-4, 4), 5, 1, 'Line vs. below')
line_above <- createVolcanoPlot(diffexp, "log2FC_band_vs_above", "pvalue_band_vs_above", "diffexp_band_vs_above", c(0, 20), c(-4, 4), 5, 1, 'Line vs. above')
below_above <- createVolcanoPlot(diffexp, "log2FC_below_vs_above", "pvalue_below_vs_above", "diffexp_below_vs_above", c(0, 60), c(-6, 6), 10, 3, 'Below vs. above')

# Create plots with labels
line_below_lab <- createVolcanoPlot(diffexp, "log2FC_band_vs_below", "pvalue_band_vs_below", "diffexp_band_vs_below", c(0, 15), c(-4, 4), 5, 1, 'Line vs. below', "labs_line_below")
line_above_lab <- createVolcanoPlot(diffexp, "log2FC_band_vs_above", "pvalue_band_vs_above", "diffexp_band_vs_above", c(0, 20), c(-4, 4), 5, 1, 'Line vs. above', "labs_line_above")
below_above_lab <- createVolcanoPlot(diffexp, "log2FC_below_vs_above", "pvalue_below_vs_above", "diffexp_below_vs_above", c(0, 60), c(-6, 6), 10, 3, 'Below vs. above', "labs_below_above")

#arrange plots that you just made into a grid on one page
ggarrange(line_below, below_above, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

#re-arranging plots
ggarrange(line_above,                                                 # First row with scatter plot
          ggarrange(line_below, below_above, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A")   