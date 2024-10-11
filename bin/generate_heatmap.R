# R script to generate heatmap from read information and gene 

library(tidyverse)
library(viridis)
library(tidyheatmaps)
library(optparse)

option_list = list(
    make_option(c("-i", "--input"), type="character", help="input virus taxa matrix to generate heatmap of"),
    make_option(c("-c", "--color"), type="character", default= "viridis", help="change color of heatmap" )
)
# Parse the arguments
opt = parse_args(OptionParser(option_list=option_list))

taxa_matrix = opt$input
color = opt$color

process_tally_file <- function(file){
    tally_vals <- read.table(taxa_matrix) %>%
        pivot_longer(cols= -sample_id, names_to="taxid", values_to="tally_count")
    tidy_heatmap(tally_vals, rows = scientific_name, columns = sample_id, values = tally_count , colors =viridis(), filename = "heatmap.pdf")                      
}

