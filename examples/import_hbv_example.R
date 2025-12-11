# =============================================================================
# Example: Import HBV Test Data to HAV Database
# =============================================================================
# This script demonstrates how to import sequences with metadata

# Set working directory to hav_db folder
setwd("c:/Users/jonr/OneDrive - Folkehelseinstituttet/Prosjekter/ngs_scripts/hav_db")

# Load the helper functions
source("R/db_functions.R")
source("R/sequence_functions.R")

# =============================================================================
# Step 1: Validate the metadata file
# =============================================================================

# First, let's validate our metadata file to check for issues
validation <- validate_metadata(
  metadata_file = "examples/hbv_test_metadata.csv",
  fasta_file = "examples/hbv_test_sequences.fasta"
)

cat("=== Metadata Validation ===\n")
cat("Valid:", validation$valid, "\n")
cat("Number of samples:", validation$n_samples, "\n")
cat("Columns:", paste(validation$columns, collapse = ", "), "\n")

if (length(validation$errors) > 0) {
  cat("\nErrors:\n")
  for (e in validation$errors) cat("  - ", e, "\n")
}

if (length(validation$warnings) > 0) {
  cat("\nWarnings:\n")
  for (w in validation$warnings) cat("  - ", w, "\n")
}

# =============================================================================
# Step 2: Connect to the database
# =============================================================================

# Get connection to the database
con <- get_db_connection("data/hav.db")

# =============================================================================
# Step 3: Import sequences with metadata
# =============================================================================

# Import using the import_with_metadata function
# This will ONLY import sequences that have a matching entry in the metadata file

cat("\n=== Starting Import ===\n")

result <- import_with_metadata(
  con = con,
  fasta_file = "examples/hbv_test_sequences.fasta",
  metadata_file = "examples/hbv_test_metadata.csv",
  skip_existing = TRUE,  # Skip if sample_id already exists
  verbose = TRUE         # Print progress
)

# =============================================================================
# Step 4: Verify the import
# =============================================================================

cat("\n=== Verifying Import ===\n")

# Count sequences
n_seq <- dbGetQuery(con, "SELECT COUNT(*) as n FROM sequences")$n
cat("Total sequences in database:", n_seq, "\n")

# Count by virus
by_virus <- dbGetQuery(con, "
  SELECT virus_name, COUNT(*) as n 
  FROM metadata 
  WHERE virus_name IS NOT NULL 
  GROUP BY virus_name
")
cat("\nSequences by virus:\n")
print(by_virus)

# Count by county
by_county <- dbGetQuery(con, "
  SELECT geo_county, COUNT(*) as n 
  FROM metadata 
  GROUP BY geo_county
")
cat("\nSequences by county:\n")
print(by_county)

# Show all sequences with metadata
cat("\n=== All Sequences ===\n")
all_data <- get_sequences_with_metadata(con)
# Don't print full sequences
all_data$sequence <- substr(all_data$sequence, 1, 50)
all_data$sequence <- paste0(all_data$sequence, "...")
print(all_data[, c("sample_id", "sequence_length", "virus_name", "geo_location", "sampling_date")])

# =============================================================================
# Step 5: Clean up
# =============================================================================

close_db_connection(con)

cat("\n=== Import Complete ===\n")
cat("Successfully imported", result$imported, "sequences with metadata.\n")
