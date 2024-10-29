library(BiocManager)
library(BiocGenerics)
library(parallel)
library(Biostrings)
library(tidyr)

setwd("")

# Function to calculate pairwise identity
calculate_identity <- function(query_seq, target_seq) {
  alignment <- pairwiseAlignment(query_seq, target_seq)
  identity <- pid(alignment, type = "PID1") # Percent identity
  return(identity)
}

# Function to perform pairwise identity comparison
pairwise_identity <- function(fasta_file, query_seq_id) {
  # Read the sequences from the FASTA file
  sequences <- readAAStringSet(fasta_file)
  
  # Extract the query sequence
  query_seq <- sequences[[query_seq_id]]
  
  # Initialize a vector to store results
  identity_results <- vector("list", length(sequences))
  names(identity_results) <- names(sequences)
  
  # Perform pairwise comparison
  for (i in seq_along(sequences)) {
    target_seq <- sequences[[i]]
    identity <- calculate_identity(query_seq, target_seq)
    identity_results[[i]] <- identity
  }
  
  return(identity_results)
}

# Usage
fasta_file <- ".fasta"
query_seq_id <- ""  # Replace with the actual ID of your query sequence

# Run the pairwise identity comparison
identity_results <- pairwise_identity(fasta_file, query_seq_id)
identity_results <- data.frame(identity_results) %>% t()

#save as CSV
write.csv(identity_results, ".csv", row.names = TRUE)

