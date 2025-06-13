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
# # header line with informative metadata
# TAXID scientific_name common_name kingdom tally median_evalue min_evalue max_evalue
# 
# The output will be written to stdout and will be 
# sorted by tally
#
# This script utilitzes various R libraries
# (taxonomizr and taxize) to assign common, scientific names and kingdom
#   Rachel Wu, August 8, 2024
#


#!/usr/bin/env R
# This script parses a blast (m8) or sam output file and 
# tallies the hits by taxonomy using fastmap for efficient lookups
# The best hit (highest bit score) for each query is 
# attributed to the NCBI TAXID corresponding
# to the hit's GI. 
# This R script is based off of the perl script written by Mark Stenglein for the same purpose
# 
# In case of a tie (equal bitscores), the attribution is given to the LCA node 
#
# The output is tab-delimited with these fields 
# # header line with informative metadata
# TAXID scientific_name common_name kingdom tally median_evalue min_evalue max_evalue
# 
# The output will be written to stdout and will be 
# sorted by tally
#
# This script utilitzes various R libraries
# (taxonomizr and taxize) to assign common, scientific names and kingdom
#   Rachel Wu, August 8, 2024


#!/usr/bin/env R
# Load necessary libraries
library(optparse)
library(taxonomizr)
library(taxize)
library(parallel)
library(fastmap)

# Command line options
option_list = list(
  make_option(c("-e", "--max_eval"), type="numeric", default=Inf, help="Maximum e-value cutoff"),
  make_option(c("-i", "--input"), type="character", help="Annotated BLAST output file"),
  make_option(c("-p", "--output_pct_id"), type="logical", default=TRUE, help="Output % Identity statistics (default TRUE)"),
  make_option(c("-w", "--query_weight_file"), type="character", help="File with query weights"),
  make_option(c("-c", "--out_tally_cutoff"), type="numeric", default=0, help="Tally cutoff for output (default 0)"),
  make_option(c("-n", "--ncbi_tax_db"), type="character", help="Path to sqlite database file with info from the NCBI Taxonomy database"),
  make_option(c("-f", "--filter"), type="character", help="Kingdom to filter for"),
  make_option(c("-u", "--unique_reads"), type = "numeric", help="Unique reads from sample"),
  make_option(c("-o", "--num_cores"), type = "numeric", help="Number of cores to use"),
  make_option(c("-t", "--taxonomic_level"), type= "character", help="tally at different taxonomic level such as family"),
  make_option(c("-s", "--species_tally"), type="logical", default=FALSE, help="tally at species level vs lca")
)

# Parse the arguments
opt = parse_args(OptionParser(option_list=option_list))

# Initialize variables
max_evalue <- opt$max_eval
output_suffix <- "tally"
annotated_blast_output <- opt$input
query_weight_file <- opt$query_weight_file
ncbi_tax_db <- opt$ncbi_tax_db
tally_cutoff <- opt$out_tally_cutoff
filter <- opt$filter
unique_reads <- opt$unique_reads
num_cores <- ifelse(is.null(opt$num_cores), detectCores() - 1, opt$num_cores)


# Initialize fastmap objects
queries <- fastmap()
gi_taxid_map <- fastmap()
taxid_tally <- fastmap()
query_weights <- fastmap()
taxonomy_cache <- fastmap()

# Function to read query weight file
read_query_weights <- function(file) {
  if (is.null(file)) return(NULL)
  weights <- read.table(file, header=FALSE, sep="\t", col.names=c("query", "weight"))
  for (i in 1:nrow(weights)) {
    query_weights$set(weights$query[i], weights$weight[i])
  }
  return(TRUE)
}

# Process BLAST file with fastmap
process_blast_file <- function(input_file) {
  blast_lines <- readLines(input_file)
  if (length(blast_lines) < 2) stop("Input file is empty or has no data rows")
  blast_lines <- blast_lines[-1]  # Skip header
  
  # Process each line
  for (line in blast_lines) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 13) {
      warning(sprintf("Ignoring line with %d fields (expected 13+): %s", length(fields), substr(line, 1, 50)))
      next
    }
    
    query <- fields[1]
    gi <- fields[2]
    taxid <- sub(";.*", "", fields[13])
    
    # Store GI to TaxID mapping
    gi_taxid_map$set(gi, taxid)
    
    # Initialize query if not exists
    if (!queries$has(query)) {
      queries$set(query, fastmap())
      queries$get(query)$set("best_bitscore", -Inf)
      queries$get(query)$set("best_hits", list())
    }
    
    current_query <- queries$get(query)
    bit_score <- as.numeric(fields[12])
    
    # Update best hits
    if (bit_score > current_query$get("best_bitscore")) {
      current_query$set("best_bitscore", bit_score)
      current_query$set("best_hits", list(line))
    } else if (bit_score == current_query$get("best_bitscore")) {
      current_hits <- current_query$get("best_hits")
      current_query$set("best_hits", c(current_hits, list(line)))
    }
  }
  message(paste("Processed blast file"))
  return(TRUE)
}

# Improved LCA functions with fastmap compatibility
get_lineage <- function(taxon_id, ncbi_tax_db) {
  if (is.na(taxon_id)) return(NA)
  
  cache_key <- paste0("lineage_", taxon_id)
  if (taxonomy_cache$has(cache_key)) {
    return(taxonomy_cache$get(cache_key))
  }
  
  lineage <- tryCatch({
    classification(taxon_id, db = ncbi_tax_db)[[1]]$id
  }, error = function(e) {
    warning(paste("Error getting lineage for taxon_id:", taxon_id, "Error:", e$message))
    return(NA)
  })
  
  taxonomy_cache$set(cache_key, lineage)
  return(lineage)
}

calculate_lca_of_2 <- function(node1, node2, ncbi_tax_db) {
  if (is.na(node1) || is.na(node2)) return(NA)
  if (node1 == node2) return(node1)
  
  lineage1 <- get_lineage(node1, ncbi_tax_db)
  lineage2 <- get_lineage(node2, ncbi_tax_db)
  
  if (identical(lineage1, NA) || identical(lineage2, NA)) return(NA)
  
  common_ancestors <- intersect(lineage1, lineage2)
  if (length(common_ancestors) == 0) return(NA)
  
  return(tail(common_ancestors, 1))
}

identify_lca <- function(taxids, ncbi_tax_db) {
  taxids <- unique(taxids)
  taxids <- taxids[!is.na(taxids)]
  
  if (length(taxids) == 0) return(NA)
  if (length(taxids) == 1) return(taxids[1])
  
  lca <- calculate_lca_of_2(taxids[1], taxids[2], ncbi_tax_db)
  if (is.na(lca)) return(NA)
  
  for (i in 3:length(taxids)) {
    lca <- calculate_lca_of_2(lca, taxids[i], ncbi_tax_db)
    if (is.na(lca)) return(NA)
  }
  
  return(lca)
}

# Process queries in parallel with fastmap
process_queries <- function(ncbi_tax_db) {

  query_list <- queries$keys()
  message(paste("Processing", length(query_list), "queries"))
  
  process_single_query <- function(query) {
      query_data <- queries$get(query)
      hits <- query_data$get("best_hits")
      if (length(hits) == 0) return(NULL)

      if (opt$species_tally){
        best_hit <- hits[[1]]  # Takes first best hit if multiple
        fields <- strsplit(best_hit, "\t")[[1]]
        gi <- fields[2]
        taxid <- gi_taxid_map$get(gi)
        message(paste("Processing", taxid, "taxid"))
        if (is.null(taxid) || taxid == "") return(NULL)
        
        # NEW: Get species-level taxid from lineage
        #lineage <- get_lineage(taxid, ncbi_tax_db)
        #if (is.na(lineage)) return(NULL)
        
        return(list(
          taxid = taxid,
          query = query,
          tally = if (query_weights$has(query)) query_weights$get(query) else 1,
          evalues = as.numeric(fields[11]),
          pct_ids = as.numeric(fields[3])
        ))
    }else{
      hit_taxids <- fastmap()
      evalues <- numeric()
      pct_ids <- numeric()
      
      for (hit in hits) {
        fields <- strsplit(hit, "\t")[[1]]
        gi <- fields[2]
        taxid <- gi_taxid_map$get(gi)
        if (is.null(taxid) || taxid == "") taxid <- NA
        
        current_count <- if (hit_taxids$has(taxid)) hit_taxids$get(taxid) else 0
        hit_taxids$set(taxid, current_count + 1)
        
        evalues <- c(evalues, as.numeric(fields[11]))
        pct_ids <- c(pct_ids, as.numeric(fields[3]))
      }
      
      lca_taxid <- identify_lca(hit_taxids$keys(), ncbi_tax_db)
      if (is.na(lca_taxid)) return(NULL)
      
      return(list(
        taxid = lca_taxid,
        query = query,
        tally = if (query_weights$has(query)) query_weights$get(query) else 1,
        evalues = evalues,
        pct_ids = pct_ids
      ))
    }
  }
  
  # Parallel processing
  query_results <- mclapply(query_list, function(query) {
    tryCatch({
      process_single_query(query)
    }, error = function(e) {
      warning(sprintf("Error processing query %s: %s", query, e$message))
      return(NULL)
    })
  }, mc.cores = num_cores)
  
  # Filter NULL results and aggregate
  query_results <- Filter(Negate(is.null), query_results)
  message(paste("Processing", length(query_results), "queries"))

  # Aggregate results into taxid_tally
  for (result in query_results) {
    taxid <- as.character(result$taxid)
    if (!taxid_tally$has(taxid)) {
      taxid_tally$set(taxid, fastmap())
      taxid_tally$get(taxid)$set("tally", 0)
      taxid_tally$get(taxid)$set("queries", list())
      taxid_tally$get(taxid)$set("evalues", numeric())
      taxid_tally$get(taxid)$set("pct_ids", numeric())
    }
    
    current <- taxid_tally$get(taxid)
    current$set("tally", current$get("tally") + result$tally)
    current$set("queries", c(current$get("queries"), result$query))
    current$set("evalues", c(current$get("evalues"), result$evalues))
    current$set("pct_ids", c(current$get("pct_ids"), result$pct_ids))
  }
  
  return(TRUE)
}

# Get taxonomy info with caching
get_cached_taxonomy <- function(taxid, ncbi_tax_db) {
  if (is.na(taxid)) return(NULL)
  cache_key <- paste0("taxinfo_", taxid)
  
  if (taxonomy_cache$has(cache_key)) {
    return(taxonomy_cache$get(cache_key))
  }
  
  tax_info <- tryCatch({
    getTaxonomy(taxid, ncbi_tax_db)
  }, error = function(e) {
    warning(paste("Error getting taxonomy for taxid:", taxid, "Error:", e$message))
    return(NULL)
  })
  
  if (!is.null(tax_info)) {
    taxonomy_cache$set(cache_key, tax_info)
  }
  
  return(tax_info)
}

# Output results with fastmap
output_results <- function(output_suffix) {
  output_file <- paste0(annotated_blast_output, ".", output_suffix)
  
  # Prepare output data
  output_data <- list()
  taxids <- taxid_tally$keys()
  
  for (taxid in taxids) {
    tax_info <- get_cached_taxonomy(taxid, ncbi_tax_db)
    if (is.null(tax_info)) next
    
    current <- taxid_tally$get(taxid)
    tally_data <- list(
      TAXID = taxid,
      Scientific_Name = if (!is.na(tax_info[1,7])) tax_info[1,7] else "Unknown",
      Common_Name = "Unknown",  # Will be updated below
      Family = if (!is.na(tax_info[1,5])) tax_info[1,5] else "Unknown",
      Kingdom = if (!is.na(tax_info[1,1])) tax_info[1,1] else "Unknown",
      Tally = current$get("tally"),
      Normalized_tally = round(current$get("tally")/unique_reads * 1000000),
      Median_evalue = median(current$get("evalues"), na.rm = TRUE),
      Min_evalue = min(current$get("evalues"), na.rm = TRUE),
      Max_evalue = max(current$get("evalues"), na.rm = TRUE),
      Median_pct_id = median(current$get("pct_ids"), na.rm = TRUE)
    )
    
    # Get common name if needed
    if (taxid != "1") {  # Skip for root
      common_name <- tryCatch({
        com <- getCommon(taxid, ncbi_tax_db, c("genbank common name", "common name"))
        if (length(com) > 0 && !is.null(com[[1]]$name)) com[[1]]$name[1] else "Unknown"
      }, error = function(e) "Unknown")
      tally_data$Common_Name <- common_name
    }
    
    # Apply filter if specified
    if (!is.null(filter) && tally_data$Kingdom != filter) next
    
    output_data[[length(output_data) + 1]] <- tally_data
  }
  
  # Convert to data frame and sort
  if (length(output_data) == 0) {
    warning("No results to output after filtering")
    return(FALSE)
  }
  
  output_df <- do.call(rbind, lapply(output_data, as.data.frame))
  output_df <- output_df[order(-output_df$Tally), ]
  
  # Write output
  write.table(output_df, file=output_file, sep="\t", row.names=FALSE, quote=FALSE)
  return(TRUE)
}

# Main execution flow
tryCatch({
  # Read query weights if provided
  if (!is.null(query_weight_file)) {
    read_query_weights(query_weight_file)
  }
  
  # Process BLAST file
  process_blast_file(annotated_blast_output)
  
  # Process queries and calculate LCAs
  process_queries(ncbi_tax_db)

  #save taxonomy cache as new rds file to use in distribution step
  saveRDS(taxonomy_cache, file="tax_cache.Rds")

  # Output results
  output_results(output_suffix)
})