#!/usr/bin/env Rscript

# A program to split up reads based on taxid of best blast hit(s)
#
# Idea is based on Sharon Chao Wooton's distributeReads.pl script
#
# Takes 2 files as input:
# (1) a fasta file used as a blast input
# (2) the output of the blast (-m8 format)
#
# Program flow:
# (1) parse blast output, creating a map of blast query id -> taxid of best hit
#     and a map of taxid->description
#
# (2) parse fasta file, for each record with a hit, output it to a new file for the
#     taxid of the best hit
#
# Original perl script created by Mark Stenglein August 13, 2011
# Updated: 4/17/2014
# Converted to R by Rachel Wu 7/20/2024

library(optparse)
library(RSQLite)
library(Biostrings)
library(taxonomizr)
library(parallel)


# Define command-line options
option_list = list(
  make_option(c("-k", "--kingdom"), type="character", help="Filter for a kingdom"),
  make_option(c("-x", "--taxid_file"), type="character", help="File containing taxids to filter"),
  make_option(c("-c", "--min_tally"), type="numeric", default=0, help="Tally cutoff for output (default 0)"),
  make_option(c("-n", "--ncbi_tax_db"), type="character", help="Path to NCBI Taxonomy sqlite database"),
  make_option(c("-f", "--fasta_file"), type="character", help="File containing fasta sequences"),
  make_option(c("-b", "--blast_file"), type="character", help="File containing blast output"),
  make_option(c("-t", "--taxids"), type="character", help="Specific taxids to output"),
  make_option(c("-d", "--diamond"), type="logical", default=FALSE, help="Blastx vs Blastn")
)

opt = parse_args(OptionParser(option_list=option_list))
num_cores = detectCores() -1
ncbi_tax_db = opt$ncbi_tax_db
filter = opt$kingdom
# Handle input files
fasta_file = opt$fasta_file
blast_file = opt$blast_file

# Handle taxids to filter
taxids_to_filter <- c()
if (!is.null(opt$taxid_file)) {
  taxids_to_filter <- readLines(opt$taxid_file)
}
if (!is.null(opt$taxids)) {
  taxids_to_filter <- c(taxids_to_filter, strsplit(opt$taxids, ",")[[1]])
}
taxids_to_filter <- as.numeric(taxids_to_filter)
output_taxid_subset <- length(taxids_to_filter) > 0

# Function to get descendants of a taxid
get_descendants <- function(taxid) {
  taxid_info::getDescendants(taxid, ncbi_tax_db)
}

higher_level_taxids <- sapply(taxids_to_filter, function(t) {
  descendants <- get_descendants(t)
  length(descendants) > 0
})

# Parse BLAST output and create a map of queries to best hits
blast_results <- read.table(blast_file, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
if(opt$diamond != FALSE){
colnames(blast_results) <- c("query_id", "subject_id", "percent_identity", "alignment_length", "mismatches", 
                             "gap_openings", "query_start", "query_end", "subject_start", "subject_end", 
                             "evalue", "bit_score", "taxid", "scientific_name", "super_kingdom")
}else{
colnames(blast_results) <- c("query_id", "subject_id", "percent_identity", "alignment_length", "mismatches", 
                             "gap_openings", "query_start", "query_end", "subject_start", "subject_end", 
                             "evalue", "bit_score", "taxid", "scientific_name", "common_name", "blast_name", 
                             "super_kingdom")
}

queries <- list()

for (i in 1:nrow(blast_results)) {
  row <- blast_results[i, ]
  query <- row$query_id
  acc <- row$subject_id
  bit_score <- row$bit_score
  scientific_name <- row$bit_score
  scientific_name <- gsub(" ", "_", scientific_name, fixed=TRUE)
  scientific_name <- gsub("/", "_", scientific_name, fixed = TRUE)
  scientific_name <- gsub("[", "_", scientific_name, fixed = TRUE)
  scientific_name <- gsub("]", "_", scientific_name, fixed = TRUE)
  scientific_name <- gsub("(", "_", scientific_name, fixed = TRUE)
  scientific_name <- gsub(")", "_", scientific_name, fixed = TRUE)
  scientific_name <- gsub("#", "_", scientific_name, fixed = TRUE)
  super_kingdom <- row$super_kingdom
  if(!is.na(as.numeric(row$taxid))){
      taxid <- row$taxid
  }else{
      next  
  }
  
  if (!is.null(queries[[query]])) {
    if (bit_score > queries[[query]]$best_bitscore) {
      queries[[query]]$best_accs <- list(acc)
      queries[[query]]$best_taxids <- list(taxid)
      queries[[query]]$best_bitscore <- bit_score
      queries[[query]]$scientific_name <- scientific_name
      queries[[query]]$super_kingdom <- super_kingdom
    } else if (bit_score == queries[[query]]$best_bitscore) {
      queries[[query]]$best_accs <- c(queries[[query]]$best_accs, acc)
      queries[[query]]$best_taxids <- c(queries[[query]]$best_taxids, taxid)
    }
  } else {
    queries[[query]] <- list(best_accs = list(acc), best_taxids = list(taxid), best_bitscore = bit_score)
  }
}

# Tally scores for each taxid
taxid_tally <- table(unlist(mclapply(queries, function(q) q$best_taxids, mc.cores=num_cores)))

# Initialize an empty cache
taxonomy_cache <- readRDS("tax_cache.Rds")

# Cached version of getTaxonomy
get_cached_taxonomy <- function(taxid, ncbi_tax_db) {
  # Check if the result is already in the cache
  if (taxonomy_cache$has(as.character(taxid))) {
    return(taxonomy_cache$get(as.character(taxid)))
  }
  
  # If not in the cache, query the database
  taxonomy_info <- getTaxonomy(taxid, ncbi_tax_db)
  
  # Store the result in the cache
  taxonomy_cache$set(as.character(taxid), taxonomy_info)
  
  # Return the result
  return(taxonomy_info)
}
# Read FASTA file and filter based on taxid
fasta_sequences <- readDNAStringSet(fasta_file, format="fasta")
for (i in 1:length(fasta_sequences)) {
  query <- names(fasta_sequences)[i]
  if (!is.null(queries[[query]])) {
    for (taxid in queries[[query]]$best_taxids) {
      if ((output_taxid_subset && taxid %in% taxids_to_filter) || !output_taxid_subset) {
        if (opt$min_tally != 0 && taxid_tally[taxid] < opt$min_tally) next
        #get names from taxonomizr
        #faster way to do this?
        taxonomy_info <- get_cached_taxonomy(taxid, ncbi_tax_db)
        scientific_name <- ifelse(is.null(taxonomy_info[1,7]), NA, taxonomy_info[1,7])
        kingdom <- ifelse(is.null(taxonomy_info[1,1]), NA, taxonomy_info[1,1])
        if (is.na(kingdom)|| is.na(scientific_name)) next
        
        print(paste("Kingdom:", kingdom))
        if (!is.null(filter) && kingdom != filter) next
        
        #remove potential characters that would break file name
        #TODO: get this in a nicer written way
        scientific_name <- gsub(" ", "_", scientific_name, fixed=TRUE)
        scientific_name <- gsub("/", "_", scientific_name, fixed = TRUE)
        scientific_name <- gsub("[", "_", scientific_name, fixed = TRUE)
        scientific_name <- gsub("]", "_", scientific_name, fixed = TRUE)
        scientific_name <- gsub("(", "_", scientific_name, fixed = TRUE)
        scientific_name <- gsub(")", "_", scientific_name, fixed = TRUE)
        scientific_name <- gsub("#", "_", scientific_name, fixed = TRUE)
        filename <- paste0(fasta_file,"_",taxid,"_",  scientific_name, ".fa" )
        filename <- gsub("'", "", filename, fixed = TRUE)
        writeXStringSet(fasta_sequences[i], filepath = filename, append = TRUE)
      }
    }
  }
}


