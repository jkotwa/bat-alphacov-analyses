library(BiocManager)
library(BiocGenerics)
library(parallel)
library(Biostrings)

# Define the path to the folder and the list of genes
folder_path <- ""
genes <- c("wg", "orf1ab", "s", "orf3", "e", "m", "n", "orf7") # Ensure genes here match genes of interest

# Function to remove non-ATCG characters
remove_non_ATCG <- function(dna_string) {
  valid_chars <- c("A", "T", "C", "G")
  dna_chars <- as.character(dna_string)
  cleaned_chars <- sapply(strsplit(dna_chars, ""), function(x) {
    paste(x[x %in% valid_chars], collapse = "")
  })
  DNAStringSet(cleaned_chars)
}

# Initialize a data frame to store results
results <- data.frame(
  identity = character(),
  gene = character(),
  percent_completeness = numeric(),
  percent_identity = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each gene
for (gene in genes) {
  # Define file paths for reference and query
  ref_file <- file.path(folder_path, paste0("ref_", gene, ".fasta")) # Reference fasta is noted as ref_gene.fasta, where gene is gene you want to compare
  query_files <- list.files(folder_path, pattern = paste0("eb_.*_", gene, ".fasta"), full.names = TRUE) # Query fasta is noted as eb_gene.fasta
  
  # Read reference sequence
  s1 <- readDNAStringSet(ref_file, format="fasta", use.names = FALSE)
  cs1 <- remove_non_ATCG(s1)
  
  # Loop through each query file
  for (query_file in query_files) {
    # Extract identity from file name
    identity <- sub(paste0(".*eb_(.*)_", gene, ".fasta"), "\\1", query_file)
    
    # Read query sequence
    s2 <- readDNAStringSet(query_file, format="fasta", use.names = FALSE)
    cs2 <- remove_non_ATCG(s2)
    
    # Perform pairwise alignment
    palign1 <- pairwiseAlignment(cs1, cs2)
    
    # Calculate % completeness
    cs1len <- as.numeric(letterFrequency(cs1, letters="ATGC"))
    cs2_completeness <- as.numeric(letterFrequency(cs2, letters="ACGT"))
    percent_completeness <- cs2_completeness/cs1len*100
    
    # Calculate percent identity
    percent_identity <- pid(palign1, type="PID2")
    
    # Append results to data frame
    results <- rbind(results, data.frame(identity = identity, gene = gene, percent_completeness = percent_completeness, percent_identity = percent_identity, stringsAsFactors = FALSE))
  }
}


# Write results to a CSV file
write.csv(results, file = file.path(folder_path, "percent_identity_results.csv"), row.names = FALSE)

# Print results
print(results)
