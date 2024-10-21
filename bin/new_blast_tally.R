#!/usr/bin/env R
# This script parses a blast (m8) or sam output file and 
# tallies the hits by taxonomy.  #
# The best hit (highest bit score) for each query is 
# attributed to the NCBI TAXID corresponding
# to the hit's GI. 
# This R script is based off of the perl script written by Mark Stenglein for the same purpose
# 
# In case of a tie (equal bitscores), the attribution is given to the LCA node 
#
# The output is tab-delimited with these fields 
# 
# TAXID scientific_name common_name kingdom tally median_evalue min_evalue max_evalue
# 
# The output will be written to stdout and will be 
# sorted by tally
#
# This script utilitzes various R libraries
# (taxonomizr and taxize) to assign common, scientific names and kingdom
#   Rachel Wu, August 8, 2024
#


# Load necessary libraries
library(optparse)
library(taxonomizr)
library(taxize)
# Command line options
option_list = list(
  make_option(c("-e", "--max_eval"), type="numeric", default=Inf, help="Maximum e-value cutoff"),
  make_option(c("-i", "--input"), type="character", help="Annotated BLAST output file"),
  make_option(c("-p", "--output_pct_id"), type="logical", default=TRUE, help="Output % Identity statistics (default TRUE)"),
  make_option(c("-w", "--query_weight_file"), type="character", help="File with query weights"),
  make_option(c("-c", "--out_tally_cutoff"), type="numeric", default=0, help="Tally cutoff for output (default 0)"),
  make_option(c("-n", "--ncbi_tax_db"), type="character", help="Path to sqlite database file with info from the NCBI Taxonomy database"),
  make_option(c("-f", "--filter"), type="character", help="Kingdom to filter for"),
  make_option(c("-u", "--unique_reads"), type = "numeric", help="Unique reads from sample")
)

# Parse the arguments
opt = parse_args(OptionParser(option_list=option_list))

# Variables initialized based on options
max_evalue = opt$max_eval
output_suffix = "tally"
annotated_blast_output = opt$input
query_weight_file = opt$query_weight_file
ncbi_tax_db = opt$ncbi_tax_db
tally_cutoff = opt$out_tally_cutoff
filter = opt$filter
unique_reads = opt$unique_reads


# Function to read query weight file (if provided)
read_query_weights <- function(file) {
  if (is.null(file)) return(NULL)
  weights <- read.table(opt$query_weights, header=FALSE, sep="\t", col.names=c("query", "weight"))
  return(weights)
}

# Read BLAST output and process
# Function to read and process each file
process_blast_file <- function(input_file) {
  # Read lines from the file
  blast_lines <- readLines(input_file)

  #Skip header line
  blast_lines <- blast_lines[-1]
  
  # Initialize necessary variables
  queries <- list()
  gi_taxid_map <- list()
  taxid_tally <- list()
  query_weights <- list()
  
  # Process each line
  for (line in blast_lines) {
    fields <- strsplit(line, "\t")[[1]]
    
    if (length(fields) < 13) {
      warning(sprintf("Ignoring line with unexpected number of fields: %s", line))
      next
    }
    
    query <- fields[1]
    gi <- fields[2]
    evalue <- sprintf("%0.1e", as.numeric(fields[11]))
    pct_id <- as.numeric(fields[3])
    taxid <- fields[13]
    
    if (is.null(taxid)) {
      warning(sprintf("Warning: TAXID undefined for accession: %s\n", gi))
      next
    }
    
    # Store GI to TaxID mapping
    gi_taxid_map[[gi]] <- taxid
    
    # Process best hits based on criteria
    if (is.null(queries[[query]])) {
      queries[[query]] <- list(best_hits = list(), best_bitscore = NULL)
    }
    
    
      best_bitscore <- queries[[query]]$best_bitscore
      bit_score <- as.numeric(fields[11])
      
      if (is.null(queries[[query]]$best_bitscore) || bit_score >= queries[[query]]$best_bitscore) {
        queries[[query]]$best_bitscore <- bit_score
        if (!is.null(queries[[query]]$best_hits)) {
            queries[[query]]$best_hits <- list()
        }
    queries[[query]]$best_hits <- append(queries[[query]]$best_hits, list(line))
    }

  }
  
  # Iterate through queries and calculate tallies
  for (query in names(queries)) {
    hits <- queries[[query]]$best_hits
    number_hits <- length(hits)
    evalue_sum <- 0
    pct_id_sum <- 0
    hit_taxids <- list()
    
    for (hit in hits) {
      fields <- strsplit(hit, "\t")[[1]]
      gi <- fields[2]
      taxid <- gi_taxid_map[[gi]]
      
      if (is.na(taxid)) {
        warning(sprintf("Warning: TAXID undefined for accession: %s. Setting this taxid to root", gi))
        taxid <- 1
      }
      
      if (is.null(hit_taxids[[taxid]])) {
        hit_taxids[[taxid]] <- 1
      } else {
        hit_taxids[[taxid]] <- hit_taxids[[taxid]] + 1
      }
      
      evalue_sum <- evalue_sum + as.numeric(fields[11])
      pct_id_sum <- pct_id_sum + as.numeric(fields[3])
    }
    
    mean_evalue <- if (number_hits > 0) sprintf("%.1e", evalue_sum / number_hits) else NA
    mean_pct_id <- if (number_hits > 0) pct_id_sum / number_hits else NA
    
    lca_taxid <- identify_lca(names(hit_taxids))
    taxid_tally[[lca_taxid]]$lca <-lca_taxid
    # Non-norm tally
    if (is.null(taxid_tally[[lca_taxid]])) {
      taxid_tally[[lca_taxid]] <- list(tally = 0, queries = list(), evalues = list(), pct_ids = list())
    }
    taxid_tally[[lca_taxid]]$tally <- taxid_tally[[lca_taxid]]$tally + query_weights[[query]]
    taxid_tally[[lca_taxid]]$queries <- c(taxid_tally[[lca_taxid]]$queries, query)
    taxid_tally[[lca_taxid]]$evalues <- c(taxid_tally[[lca_taxid]]$evalues, mean_evalue)
    taxid_tally[[lca_taxid]]$pct_ids <- c(taxid_tally[[lca_taxid]]$pct_ids, mean_pct_id)
  }
  #calculate statistics for taxid_tally
  #does this work with lca??
  for (taxid in names(taxid_tally)) {
    evalues <- as.numeric(taxid_tally[[taxid]]$evalues)
    pct_ids <- as.numeric(taxid_tally[[taxid]]$pct_ids)
    
    if (length(evalues) > 0) {
      taxid_tally[[taxid]]$median_evalue <- median(evalues, na.rm = TRUE)
      taxid_tally[[taxid]]$min_evalue <- min(evalues, na.rm = TRUE)
      taxid_tally[[taxid]]$max_evalue <- max(evalues, na.rm = TRUE)
    } else {
      taxid_tally[[taxid]]$median_evalue <- NA
      taxid_tally[[taxid]]$min_evalue <- NA
      taxid_tally[[taxid]]$max_evalue <- NA
    }
    
    if (length(pct_ids) > 0) {
      taxid_tally[[taxid]]$median_pct_id <- median(pct_ids, na.rm = TRUE)
      taxid_tally[[taxid]]$min_pct_id <- min(pct_ids, na.rm = TRUE)
      taxid_tally[[taxid]]$max_pct_id <- max(pct_ids, na.rm = TRUE)
    } else {
      taxid_tally[[taxid]]$median_pct_id <- NA
      taxid_tally[[taxid]]$min_pct_id <- NA
      taxid_tally[[taxid]]$max_pct_id <- NA
    }
  }
  
  return(taxid_tally)
}

# Function to get lineage of a taxon
get_lineage <- function(taxon_id) {
  lineage <- classification(taxon_id, db = ncbi_tax_db)
  lineage <- lineage[[1]]$id
  return(lineage)
}

# Function to calculate LCA of two nodes
calculate_lca_of_2 <- function(node1, node2) {
  
  # Get lineages for both nodes
  lineage1 <- get_lineage(node1)
  lineage2 <- get_lineage(node2)
  
  # Find common ancestors
  common_ancestors <- intersect(lineage1, lineage2)
  
  # The LCA is the most recent common ancestor (last common in lineage)
  if (length(common_ancestors) > 0) {
    lca <- tail(common_ancestors, 1)
    return(lca)
  } else {
    warning("No common ancestor found.")
    return(NULL)
  }
}

# Recursive function to calculate LCA of a list of nodes
identify_lca <- function(taxids) {
  if (length(taxids) == 0) {
    return(NULL)
  } else if (length(taxids) == 1) {
    return(taxids[1])
  }
  
  # Start with the first two taxids
  lca <- calculate_lca_of_2(taxids[1], taxids[2])
  
  # Recursively calculate LCA with the remaining nodes
  for (i in 3:length(taxids)) {
    lca <- calculate_lca_of_2(lca, taxids[i])
  }
  
  return(lca)
}

# Process the file
taxid_tally <- process_blast_file(annotated_blast_output)



# Output results to file
output_results <- function(taxid_tally, output_suffix) {
  output_file <- paste0(annotated_blast_output, ".", output_suffix)
  output_df <- data.frame(TAXID = character(),
                          Scientific_Name = character(),
                          Common_Name = character(),
                          Kingdom = character(),
                          Tally = numeric(),
                          Normalized_tally = numeric(),
                          Median_evalue = numeric(),
                          Min_evalue = numeric(),
                          Max_evalue = numeric(),
                          Median_pct_id = numeric()
                          #Min_pct_id = numeric(),
                          #ax_pct_id = numeric()
                          )
  str(taxid_tally)
  for (taxid in names(taxid_tally)) {
    taxonomy_info <- getTaxonomy(taxid, ncbi_tax_db )
    scientific_name <- ifelse(is.null(taxonomy_info[1,7]), NA, taxonomy_info[1,7])
    scientific_name <- gsub("'", "", scientific_name, fixed = TRUE)
    scientific_name <- gsub(" ", "_", scientific_name, fixed = TRUE)
    if(is.na(scientific_name)) next
    common_info <- getCommon(taxid, ncbi_tax_db, c("genbank common name", "common name"))
    
    # Handle cases where common_info might be empty
    if (length(common_info) > 0 && !is.null(common_info[[1]]$name)) {
      common_name <- common_info[[1]]$name[1]
      common_name <- gsub("'", "", common_name, fixed = TRUE)
    } else {
      common_name <- NA
    }
    if (is.na(common_name)) common_name <- "Unknown"
    
    kingdom <- ifelse(is.null(taxonomy_info[1,1]), NA, taxonomy_info[1,1])
    if (!is.null(filter)){
      if (kingdom != filter){
        next
      }
    }
    
    # Check if any of these are NA, and if so, handle accordingly
    tally <- length(taxid_tally[[taxid]]$queries)
    normalized_tally <- round(tally/unique_reads * 1000000)
    median_evalue <- taxid_tally[[taxid]]$median_evalue
    min_evalue <- taxid_tally[[taxid]]$min_evalue
    max_evalue <- taxid_tally[[taxid]]$max_evalue
    median_pct_id <- ifelse(!is.null(taxid_tally[[taxid]]$median_pct_id), taxid_tally[[taxid]]$median_pct_id, NA)
    #min_pct_id <- ifelse(!is.null(taxid_tally[[taxid]]$min_pct_id), taxid_tally[[taxid]]$min_pct_id, NA)
    #max_pct_id <- ifelse(!is.null(taxid_tally[[taxid]]$max_pct_id), taxid_tally[[taxid]]$max_pct_id, NA)
    print(paste("Scientific Name:", scientific_name))
    print(paste("Common Name:", common_name))
    print(paste("Kingdom:", kingdom))
    #print(paste("LCA:", lca))
    
    output_df <- rbind(output_df, data.frame(
      TAXID = taxid,
      Scientific_Name = scientific_name,
      Common_Name = common_name,
      Kingdom = kingdom,
      Tally = tally,
      Normalized_tally = normalized_tally,
      Median_evalue = median_evalue,
      Min_evalue = min_evalue,
      Max_evalue = max_evalue,
      Median_pct_id = median_pct_id
      #Min_pct_id = min_pct_id,
      #Max_pct_id = max_pct_id
    ))
  }
  # Sort by Tally (number of hits) in descending order
  output_df <- output_df[order(-output_df$Tally), ]

  output_file <- write.table(output_df, file=output_file, sep="\t", row.names=FALSE, quote=FALSE)
}


# Call output_results function
output_results(taxid_tally, output_suffix)

