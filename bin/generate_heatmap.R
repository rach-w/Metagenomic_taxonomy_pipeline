# R script to generate heatmap from read information and gene 

library(tidyverse)
library(viridis)
library(tidyheatmaps)
library(optparse)
library(ggplot2)
library(stringr)


option_list = list(
    make_option(c("-i", "--input"), type="character", help="input virus taxa matrix to generate heatmap of"),
    make_option(c("-c", "--color"), type="character", default= "viridis", help="change color of heatmap" ),
    make_option(c("-f", "--family"), type="logical", default=FALSE, help= "plot family level" )
)
# Parse the arguments
opt = parse_args(OptionParser(option_list=option_list))

taxa_matrix = opt$input
color = opt$color
plot_family = opt$family
#read tally matrix file

tally_vals <- read.table(taxa_matrix, header = TRUE, stringsAsFactors = FALSE, sep = "\t")



# Assuming you have a column in your data called 'Tool' that indicates whether it's DIAMOND or BLASTn
# If not, you can create one based on the 'Blastn_color' column

# Generate heatmap
# Generate heatmap
if (plot_family) {
  p <- ggplot(tally_vals, aes(x = Family, y = Barcode)) +
    geom_tile(aes(fill = Normalized_tally)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 5)) +
    scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 90, color = tally_vals$Blastn_color, size = 3)) +
    theme(axis.text.y = element_text(angle = 70, size = 4)) +
    theme(legend.position = "bottom", legend.key.width = unit(2, 'cm')) +
    labs(caption = "Red: DIAMOND | Blue: BLASTn") +  # Add caption below the plot
    theme(plot.caption = element_text(hjust = 0.5, size = 10, color = "black"))  # Customize caption
} else {
  p <- ggplot(tally_vals, aes(x = Scientific_Name, y = Barcode)) +
    geom_tile(aes(fill = Normalized_tally)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 5)) +
    scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 90, color = tally_vals$Blastn_color, size = 3)) +
    theme(axis.text.y = element_text(angle = 0, size = 4)) +
    theme(legend.position = "bottom", legend.key.width = unit(2, 'cm')) +
    labs(caption = "Red: DIAMOND | Blue: BLASTn") +  # Add caption below the plot
    theme(plot.caption = element_text(hjust = 0.5, size = 10, color = "black"))  # Customize caption
}



# Print the plot
print(p)




