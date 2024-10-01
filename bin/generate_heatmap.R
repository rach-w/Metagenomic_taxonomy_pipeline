# R script to generate heatmap from read information and gene 

library(tidyverse)
library(viridis)
library(tidyheatmaps)
library(optparse)

option_list = list(
    make_option(c("-i", "--input"), type="character", help="input tally table to generate heatmap of"),
    mkae_option(c("-c", "--color"), type="character", default= "viridis", help="change color of heatmap" )
)
# Parse the arguments
opt = parse_args(OptionParser(option_list=option_list))

tally_file = opt$input
color = opt$color


