# =============================================================================
# Example: Importing Sequences with Metadata
# =============================================================================
# This script demonstrates how to format metadata and import sequences
# into the HAV database.

library(DBI)
library(RSQLite)

# Set working directory to hav_db folder
# setwd("path/to/hav_db")

# Source the helper functions
source("R/db_functions.R")
source("R/sequence_functions.R")

# Connect to database
con <- get_db_connection("data/hav.db")

# =============================================================================
# STEP 1: Prepare your FASTA file
# =============================================================================
# Your multi-FASTA file should look like this:
#
# >SAMPLE_001
# ATGCGATCGATCGATCGATCGATCGATCGATCG...
# >SAMPLE_002
# ATGCGATCGATCGATCGATCGATCGATCGATCG...
# >SAMPLE_003
# ATGCGATCGATCGATCGATCGATCGATCGATCG...
# >SAMPLE_004
# ATGCGATCGATCGATCGATCGATCGATCGATCG...
#
# The header (after >) is the sample_id that must match your metadata.

# =============================================================================
# STEP 2: Prepare your metadata file (CSV)
# =============================================================================
# The metadata file is a CSV with sample_id as the only REQUIRED column.
# All other columns are optional.
#
# Column names:
#   sample_id       - REQUIRED - must match FASTA headers exactly
#   sampling_date   - Optional - format: YYYY-MM-DD (e.g., 2025-01-15)
#   sequence_length - Optional - informational only, actual length calculated from sequence
#   geo_location    - Optional - city/place name
#   geo_county      - Optional - county/region
#   geo_country     - Optional - country (defaults to "Norway" if not provided)
#   organism        - Optional - species/strain name
#   gene_target     - Optional - gene name (e.g., HA, NA, NS1)
#   submitter       - Optional - who submitted the sample
#   sequencing_run  - Optional - sequencing run identifier
#   notes           - Optional - any additional notes

# Example: Create a sample metadata data frame
metadata <- data.frame(
  sample_id = c("SAMPLE_001", "SAMPLE_002", "SAMPLE_003"),
  sampling_date = c("2025-01-15", "2025-01-20", "2025-02-01"),
  sequence_length = c(1500, 1480, 1520),  # informational
  geo_location = c("Oslo", "Bergen", "Trondheim"),
  geo_county = c("Oslo", "Vestland", "Trøndelag"),
  geo_country = c("Norway", "Norway", "Norway"),
  organism = c("Influenza A", "Influenza A", "Influenza B"),
  gene_target = c("HA", "HA", "HA"),
  submitter = c("Lab1", "Lab1", "Lab2"),
  sequencing_run = c("RUN_2025_001", "RUN_2025_001", "RUN_2025_002"),
  notes = c("", "Good quality", "Low coverage region at 3' end"),
  stringsAsFactors = FALSE
)

# View the structure
print(metadata)

# Save to CSV (this is what your actual metadata file should look like)
write.csv(metadata, "example_metadata.csv", row.names = FALSE)
cat("Example metadata saved to: example_metadata.csv\n")

# =============================================================================
# STEP 3: Create example FASTA file (for demonstration)
# =============================================================================
# In practice, you would already have your FASTA file from sequencing

example_fasta <- c(
  ">SAMPLE_001",
  "ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAATGCAGACACATTATGTATAGGTTATCATGCG",
  ">SAMPLE_002", 
  "ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTACAACCGCAAATGCAGACACATTATGTATAGGTTATCATGCG",
  ">SAMPLE_003",
  "ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAATGCAGACACATTATGTATAGGTTATCATGCG",
  ">SAMPLE_004",
  "ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAATGCAGACACATTATGTATAGGTTATCATGCG",
  ">SAMPLE_005",
  "ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAATGCAGACACATTATGTATAGGTTATCATGCG"
)

writeLines(example_fasta, "example_sequences.fasta")
cat("Example FASTA saved to: example_sequences.fasta\n")
cat("Note: FASTA has 5 sequences, but metadata only has 3 - only those 3 will be imported\n\n")

# =============================================================================
# STEP 4: Validate metadata before import
# =============================================================================
cat("Validating metadata...\n")
validation <- validate_metadata("example_metadata.csv", "example_sequences.fasta")

if (validation$valid) {
  cat("✓ Metadata is valid!\n")
} else {
  cat("✗ Metadata has errors:\n")
  for (e in validation$errors) cat("  ERROR:", e, "\n")
}

if (length(validation$warnings) > 0) {
  cat("Warnings:\n")
  for (w in validation$warnings) cat("  WARNING:", w, "\n")
}

cat("\nSamples in metadata:", validation$n_samples, "\n")
cat("Columns:", paste(validation$columns, collapse = ", "), "\n\n")

# =============================================================================
# STEP 5: Import sequences
# =============================================================================
cat("Importing sequences...\n\n")

result <- import_with_metadata(
  con = con,
  fasta_file = "example_sequences.fasta",
  metadata_file = "example_metadata.csv",
  skip_existing = TRUE,
  verbose = TRUE
)

# =============================================================================
# STEP 6: Verify import
# =============================================================================
cat("\n\nVerifying import...\n")

# Check what's in the database
all_seqs <- get_all_sequences(con)
cat("Total sequences in database:", nrow(all_seqs), "\n")

# View sequences with metadata
data <- get_sequences_with_metadata(con)
print(data[, c("sample_id", "sequence_length", "organism", "gene_target", "sampling_date", "geo_location")])

# =============================================================================
# STEP 7: Clean up (optional - for this demo only)
# =============================================================================
# Uncomment to delete the example data:
# delete_sequence(con, "SAMPLE_001")
# delete_sequence(con, "SAMPLE_002")
# delete_sequence(con, "SAMPLE_003")

# Close connection
close_db_connection(con)

cat("\nDone! You can now open the Shiny app to view the imported sequences.\n")
cat("Run: shiny::runApp('path/to/hav_db')\n")

# =============================================================================
# ALTERNATIVE: Minimal metadata example
# =============================================================================
# If you only have sample IDs and one piece of metadata, that's fine!
# Here's a minimal example:

minimal_metadata <- data.frame(
  sample_id = c("SEQ_A", "SEQ_B", "SEQ_C"),
  virus_name = c("HAV", "HAV", "HAV")
)

# This is perfectly valid - only sample_id is required!
# write.csv(minimal_metadata, "minimal_metadata.csv", row.names = FALSE)

# =============================================================================
# Create a blank template
# =============================================================================
# You can also export a blank template to fill in:

template <- create_metadata_template()
write.csv(template, "metadata_template.csv", row.names = FALSE)
cat("\nBlank template saved to: metadata_template.csv\n")
