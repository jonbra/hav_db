# =============================================================================
# Import Clean Data into HAV Database
# =============================================================================
# This script imports the cleaned and parsed metadata into the database.
# Run parse_metadata.R first to generate the cleaned CSV files.
#
# Input:  data/cleaned_sequences.csv, data/cleaned_metadata.csv
# Output: Populated database at data/hav.db
# =============================================================================

library(DBI)
library(RSQLite)

# Source the parsing script to get the function
source("R/parse_metadata.R")

# =============================================================================
# Main Import Function
# =============================================================================

import_to_database <- function(db_path = "data/hav.db",
                                fasta_file = "export/2PA.fa",
                                csv_file = "export/export.csv") {
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("HAV DATABASE IMPORT\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  # Step 1: Parse and clean the data
  cat("Step 1: Parsing and cleaning data...\n")
  data <- prepare_metadata(fasta_file, csv_file)
  sequences_df <- data$sequences
  metadata_df <- data$metadata
  
  # Step 2: Connect to database
  cat("\nStep 2: Connecting to database...\n")
  if (!file.exists(db_path)) {
    stop("Database not found at: ", db_path, "\nRun setup_database.R first.")
  }
  
  con <- dbConnect(RSQLite::SQLite(), db_path)
  on.exit(dbDisconnect(con))
  
  # Step 3: Clear existing data
  cat("Step 3: Clearing existing data...\n")
  dbExecute(con, "DELETE FROM metadata")
  dbExecute(con, "DELETE FROM sequences")
  dbExecute(con, "DELETE FROM analysis_results")
  cat("  Database cleared.\n")
  
  # Step 4: Import sequences
  cat("\nStep 4: Importing sequences...\n")
  
  success_seq <- 0
  errors_seq <- 0
  
  for (i in 1:nrow(sequences_df)) {
    tryCatch({
      row <- sequences_df[i, ]
      dbExecute(con, 
        "INSERT INTO sequences (sample_id, sequence, sequence_length) VALUES (?, ?, ?)",
        params = list(row$sample_id, row$sequence, row$sequence_length)
      )
      success_seq <- success_seq + 1
    }, error = function(e) {
      errors_seq <<- errors_seq + 1
      if (errors_seq <= 5) {
        cat("  Error importing sequence", sequences_df$sample_id[i], ":", e$message, "\n")
      }
    })
    
    if (i %% 100 == 0) cat("  Imported", i, "sequences...\n")
  }
  
  cat("  Sequences imported:", success_seq, "| Errors:", errors_seq, "\n")
  
  # Step 5: Import metadata
  cat("\nStep 5: Importing metadata...\n")
  
  success_meta <- 0
  errors_meta <- 0
  
  for (i in 1:nrow(metadata_df)) {
    tryCatch({
      row <- metadata_df[i, ]
      
      # Convert NA to NULL-friendly format
      sampling_date <- if (is.na(row$sampling_date) || row$sampling_date == "") NA else row$sampling_date
      sample_year <- if (is.na(row$sample_year)) NA else as.integer(row$sample_year)
      genotype <- if (is.na(row$genotype) || row$genotype == "") NA else row$genotype
      variant <- if (is.na(row$variant) || row$variant == "") NA else row$variant
      patient_id <- if (is.na(row$patient_id) || row$patient_id == "") NA else row$patient_id
      geo_location <- if (is.na(row$geo_location) || row$geo_location == "") NA else row$geo_location
      geo_country <- if (is.na(row$geo_country) || row$geo_country == "") NA else row$geo_country
      source <- if (is.na(row$source) || row$source == "") NA else row$source
      transmission_route <- if (is.na(row$transmission_route) || row$transmission_route == "") NA else row$transmission_route
      comment <- if (is.na(row$comment) || row$comment == "") NA else row$comment
      
      dbExecute(con, "
        INSERT INTO metadata 
        (sample_id, sampling_date, sample_year, genotype, variant, patient_id, 
         geo_location, geo_country, source, transmission_route, comment)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        params = list(
          row$sample_id, sampling_date, sample_year, genotype, variant, patient_id,
          geo_location, geo_country, source, transmission_route, comment
        )
      )
      success_meta <- success_meta + 1
    }, error = function(e) {
      errors_meta <<- errors_meta + 1
      if (errors_meta <= 5) {
        cat("  Error importing metadata", metadata_df$sample_id[i], ":", e$message, "\n")
      }
    })
    
    if (i %% 100 == 0) cat("  Imported", i, "metadata records...\n")
  }
  
  cat("  Metadata imported:", success_meta, "| Errors:", errors_meta, "\n")
  
  # Step 6: Verify and show summary
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("IMPORT COMPLETE\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  # Count totals
  total_seq <- dbGetQuery(con, "SELECT COUNT(*) as n FROM sequences")$n
  total_meta <- dbGetQuery(con, "SELECT COUNT(*) as n FROM metadata")$n
  
  cat("Database statistics:\n")
  cat("  Total sequences:", total_seq, "\n")
  cat("  Total metadata records:", total_meta, "\n")
  
  # Show sample data
  cat("\nSample data (first 10 rows):\n")
  sample_data <- dbGetQuery(con, "
    SELECT s.sample_id, s.sequence_length, 
           m.sampling_date, m.sample_year, m.genotype, m.variant,
           m.patient_id, m.geo_country, m.geo_location, m.comment
    FROM sequences s
    LEFT JOIN metadata m ON s.sample_id = m.sample_id
    LIMIT 10
  ")
  print(sample_data)
  
  # Show genotype distribution
  cat("\nGenotype distribution:\n")
  genotypes <- dbGetQuery(con, "SELECT genotype, COUNT(*) as count FROM metadata GROUP BY genotype ORDER BY count DESC")
  print(genotypes)
  
  # Show country distribution
  cat("\nCountry distribution (top 15):\n")
  countries <- dbGetQuery(con, "SELECT geo_country, COUNT(*) as count FROM metadata GROUP BY geo_country ORDER BY count DESC LIMIT 15")
  print(countries)
  
  # Show year distribution
  cat("\nYear distribution:\n")
  years <- dbGetQuery(con, "SELECT sample_year, COUNT(*) as count FROM metadata WHERE sample_year IS NOT NULL GROUP BY sample_year ORDER BY sample_year")
  print(years)
  
  # Show records with comments
  cat("\nRecords with comments:\n")
  comments <- dbGetQuery(con, "SELECT COUNT(*) as count FROM metadata WHERE comment IS NOT NULL")
  cat("  Total:", comments$count, "\n")
  
  cat("\nDone!\n")
}

# =============================================================================
# Run if called directly
# =============================================================================
if (!interactive()) {
  import_to_database()
}
