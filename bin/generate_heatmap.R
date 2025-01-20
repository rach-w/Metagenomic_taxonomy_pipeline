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
    make_option(c("-f", "--family"), type="logical", default=FALSE, help"plot family level")
)
# Parse the arguments
opt = parse_args(OptionParser(option_list=option_list))

taxa_matrix = opt$input
color = opt$color
plot_family = opt$family
#read tally matrix file

tally_vals <- read.table(taxa_matrix, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

#generate heatmao
#tidy_heatmap(tally_vals, rows = TAXID, columns = Barcode, values = Normalized_tally , filename = "heatmap.pdf")  + scale_fill_viridis()                    
if(plot_family){
  ggplot(tally_vals, aes(x = Family, Barcode)) +
    geom_tile(aes(fill = Normalized_tally)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 5)) +
    scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 90, color=tally_vals$Blastn_color,size=3)) +
    theme(axis.text.y = element_text(angle = 70, size=4))+
    theme(legend.position = "bottom")
}else{
  ggplot(tally_vals, aes(x = Scientific_Name, Barcode)) +
    geom_tile(aes(fill = Normalized_tally)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 5)) +
    scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 90, color=tally_vals$Blastn_color,size=3)) +
    theme(axis.text.y = element_text(angle = 70, size=4))+
    theme(legend.position = "bottom")
}


