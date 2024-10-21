#!/usr/bin/env Rscript

# Load necessary libraries
library(optparse)

# Set up command line options
option_list <- list(
  make_option(c("-r", "--rank"),type = "logical", default = FALSE, help = "Input tally files contain an extra rank column"),
  make_option(c("-t", "--taxids_file"), type = "character", default = NULL, help = "Limit output to the taxids listed in this file"),
  make_option(c("-c", "--tally_cutoff"), type = "numeric", default = 0, help = "Only output a read count value if value is > this cutoff"),
  make_option(c("-v", "--virus_only"), type= "logical", default = FALSE, help = "Only include virus taxids in matrix"),
  make_option(c("-p", "--no_phage"), type = "logical", default = FALSE, help = "Exclude phage taxids"),
  make_option(c("-i", "--input_tally"), type = "character", help = "input tally matrices")
) 

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)



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
                        Kingdom = character(),
                        Tally = numeric()
                          )

# Process each tally file
tally_files <- unlist(strsplit(args$input_tally, ","))
for (tally_file in tally_files) {
  # Extract barcode from filename
  barcode <- sub("\\..*$", "", tally_file)
  if (barcode == tally_file) {
    barcode <- tally_file
  }
  
  if (!(barcode %in% barcodes)) {
    barcodes <- c(barcodes, barcode)
  }

  # Read tally file
  tally_data <- read.table(tally_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

  for (i in 1:nrow(tally_data)) {
    if (args$rank) {
      rank <- tally_data[i, 1]
      taxid <- tally_data[i, 2]
      scientific_name <- tally_data[i, 4]
      common_name <- tally_data[i, 5]
      kingdom <- tally_data[i, 6]
      tally <- tally_data[i, 7]
      normalized_tally <- tally_data[i, 8]
      median_evalue<- tally_data[i,9]
      min_evalue <- tally_data[i,10]
      max_evalue <- tally_data[i,11]
      median_pct_id <- tally_data[i,12]
    } else {
      taxid <- tally_data[i, 1]
      scientific_name <- tally_data[i, 2]
      common_name <- tally_data[i, 3]
      kingdom <- tally_data[i, 4]
      tally <- tally_data[i, 5]
      normalized_tally <- tally_data[i,6]
      median_evalue<- tally_data[i,7]
      min_evalue <- tally_data[i,8]
      max_evalue <- tally_data[i,9]
      median_pct_id <- tally_data[i,10]
    }

    if (nchar(taxid) == 0) {
      taxid <- "X"
      scientific_name <- "Unknown Taxid"
      kingdom <- "X"
    }

    # Filter based on options
    if (args$virus_only && kingdom != "Viruses") {
      next
    }
    if (args$no_phage && grepl("phage", scientific_name, ignore.case = TRUE)) {
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
            Kingdom = kingdom,
            Normalized_tally = normalized_tally
            ))
    
  }
}

output_file <- write.table(output, file = "taxa_matrix.tsv", sep= "\t", row.names=FALSE, quote = FALSE)






                        

