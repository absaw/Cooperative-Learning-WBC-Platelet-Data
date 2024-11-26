# Required Packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")

library(biomaRt)
library(dplyr)
library(readr)
library(data.table)

# Function to Replace ENSG IDs with Gene Symbols
replace_ensg_with_symbols <- function(input_csv, output_csv) {
  # Read the CSV file
  data <- read_csv(input_csv)
  
  # Clean ENSG IDs (Remove version numbers)
  ensg_ids <- colnames(data)
  cleaned_ensg_ids <- sub("[.][0-9]*", "", ensg_ids)
  
  # Connect to Ensembl database using biomaRt
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  # Query biomaRt to get gene symbols
  gene_IDs <- getBM(
    filters = "ensembl_gene_id", 
    attributes = c("ensembl_gene_id", "hgnc_symbol"), 
    values = cleaned_ensg_ids, 
    mart = mart
  )
  
  # Create a mapping data frame
  mapping <- data.frame(
    ensg_cleaned = cleaned_ensg_ids,
    original = ensg_ids
  ) %>%
    left_join(gene_IDs, by = c("ensg_cleaned" = "ensembl_gene_id")) %>%
    mutate(final_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "", ensg_cleaned, hgnc_symbol))
  
  # Replace column names in the original data
  colnames(data) <- mapping$final_symbol
  
  # Save the updated data to a new CSV file
  write_csv(data, output_csv)
}

# Usage
input_file <- "InputNormCountFiles/NormCountsPlatelet.csv"  # Input CSV file
output_file <- "InputNormCountFiles/NormCountsPlateletGeneSymbol.csv"  # Output CSV file
data <- read_csv(input_file)

replace_ensg_with_symbols(input_file, output_file)

# Function to Check for Null Values in Columns
check_null_columns <- function(file_path) {
  # Read the CSV file
  data <- read_csv(file_path)
  
  # Identify columns with null (NA) values
  null_columns <- colnames(data)[colSums(is.na(data)) > 0]
  
  # Output results
  if (length(null_columns) > 0) {
    message("Columns with null values:")
    print(null_columns)
  } else {
    message("No columns with null values found.")
  }
}
check_null_columns(output_file)
