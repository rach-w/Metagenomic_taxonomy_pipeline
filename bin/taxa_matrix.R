#!/usr/bin/env Rscript

# Load necessary libraries
library(argparse)
library(taxonomizr)
# Create an ArgumentParser object
parser <- ArgumentParser(description = "Process input tally files")

# Set up command line options
parser$add_argument("-r", "--rank", action = "store_true", default = FALSE,
                    help = "Input tally files contain an extra rank column")
parser$add_argument("-t", "--taxids_file", type = "character", default = NULL,
                    help = "Limit output to the taxids listed in this file")
parser$add_argument("-c", "--tally_cutoff", type = "numeric", default = 0,
                    help = "Only output a read count value if value is > this cutoff")
parser$add_argument("-v", "--virus_only", action = "store_true", default = FALSE,
                    help = "Only include virus taxids in matrix")
parser$add_argument("-p", "--no_phage", action = "store_true", default = FALSE,
                    help = "Exclude phage taxids")
parser$add_argument("-e", "--exclude_family", type = "character", default = NULL,
                    help = "Exclude a family of taxids from output")
parser$add_argument("input_tally", nargs = "+", type = "character", 
                    help = "Space-separated list of input tally matrices")
parser$add_argument("-o", "--output_file", type = "character", default="taxa_matrix.tsv", 
                    help="Name of output file to output to")
parser$add_argument("-n", "--ncbi_database", type="character",
                    help = "taxonomy database for family level lookup")


# Parse the arguments 
args <- parser$parse_args()

taxonomyDatabase <- args$ncbi_database

# Read taxids file if provided
taxids <- character()
taxid_name_map <- list()

if (!is.null(args$taxids_file)) {
  taxid_data <- read.table(args$taxids_file, header = FALSE, stringsAsFactors = FALSE)
  taxids <- taxid_data$V1
  taxid_name_map <- setNames(taxid_data$V2, taxid_data$V1)
}

# Initialize data structures
matrix_data <- list()
barcodes <- character()
output <- data.frame(
                        Barcode = character(),
                        TAXID = character(),
                        Scientific_Name = character(),
                        Common_Name = character(),
                        Family = character(),
                        Kingdom = character(),
                        Tally = numeric(),
                        Blastn_color = character()
                          )

for (tally_file in args$input_tally) {

  if(grepl("bx_nr", tally_file, fixed = TRUE) ){
    blastn <- "blastx"
  } else{
    blastn <- "blastn"
  }
  # Extract barcode from filename
  barcode <- sub("[_-].*$", "", tally_file)
  if (barcode == tally_file) {
    barcode <- tally_file
  }
  
  if (!(barcode %in% barcodes)) {
    barcodes <- c(barcodes, barcode)
  }

  # Read tally file
  print(barcode)
  tally_data <- read.table(tally_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t", fill=TRUE)


  for (i in 1:nrow(tally_data)) {
    if (args$rank) {
      rank <- tally_data[i, 1]
      taxid <- tally_data[i, 2]
      scientific_name <- tally_data[i, 3]
      common_name <- tally_data[i, 4]
      kingdom <- tally_data[i, 5]
      tally <- tally_data[i, 6]
      normalized_tally <- tally_data[i, 7]
      median_evalue<- tally_data[i,8]
      min_evalue <- tally_data[i,9]
      max_evalue <- tally_data[i,10]
      median_pct_id <- tally_data[i,11]
    } else {
      taxid <- tally_data[i, 1]
      scientific_name <- tally_data[i, 2]
      common_name <- tally_data[i, 3]
      family <- tally_data[i, 4]
      kingdom <- tally_data[i, 5]
      tally <- tally_data[i, 6]
      normalized_tally <- tally_data[i,7]
      median_evalue<- tally_data[i,8]
      min_evalue <- tally_data[i,9]
      max_evalue <- tally_data[i,10]
      median_pct_id <- tally_data[i,11]
    }

    if (nchar(taxid) == 0) {
      taxid <- "X"
      scientific_name <- "Unknown Taxid"
      kingdom <- "X"
    }

    # Filter based on options
    if( !is.na(kingdom) && !is.null(kingdom)){
      if (args$virus_only && kingdom != "Viruses") {
      next
      }
    }
    
    if (args$no_phage && grepl("phage", scientific_name, ignore.case = TRUE)) {
      next
    }
    if (!is.null(args$exclude_family) && grepl(args$exclude_family, scientific_name, ignore.case = TRUE)){
        next
    }

    if (!(taxid %in% names(taxid_name_map))) {
      taxid_name_map[[taxid]] <- scientific_name
      if (length(taxids) == 0) {
        taxids <- c(taxids, taxid)
      }
    }
    # Populate data frame
    output <- rbind(output, data.frame(
            Barcode = barcode,
            TAXID = taxid,
            Scientific_Name = scientific_name,
            Common_Name = common_name,
            Family = family,
            Kingdom = kingdom,
            Normalized_tally = normalized_tally,
            Blastn_Blastx = blastn
            ))
    
  }
}

output_file <- write.table(output, file = args$output_file, sep= "\t", row.names=FALSE, quote = FALSE)






                        

