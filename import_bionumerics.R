# =============================================================================
# Import BioNumerics data into HAV database
# =============================================================================

# Load database functions
source("R/db_functions.R")
source("R/sequence_functions.R")

# Input files
fasta_file <- "data/bionumerics_sequences.fasta"
metadata_file <- "data/bionumerics_metadata.csv"

# Connect to database
cat("Connecting to database...\n")
con <- get_db_connection("data/hav.db")

# Read metadata
cat("Reading metadata...\n")
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
cat("  Found", nrow(metadata), "metadata records\n")

# Read FASTA sequences
cat("Reading FASTA sequences...\n")
fasta_lines <- readLines(fasta_file)
sequences <- list()
current_id <- NULL
current_seq <- ""

for (line in fasta_lines) {
  if (startsWith(line, ">")) {
    if (!is.null(current_id)) {
      sequences[[current_id]] <- current_seq
    }
    current_id <- sub("^>", "", line)
    current_seq <- ""
  } else {
    current_seq <- paste0(current_seq, line)
  }
}
# Don't forget the last sequence
if (!is.null(current_id)) {
  sequences[[current_id]] <- current_seq
}
cat("  Found", length(sequences), "sequences\n")

# Check existing samples to avoid duplicates
cat("\nChecking for existing samples...\n")
existing <- dbGetQuery(con, "SELECT sample_id FROM sequences")$sample_id
new_samples <- setdiff(names(sequences), existing)
cat("  Existing samples in database:", length(existing), "\n")
cat("  New samples to import:", length(new_samples), "\n")

if (length(new_samples) == 0) {
  cat("\nNo new samples to import. All samples already exist in database.\n")
  close_db_connection(con)
  quit(save = "no")
}

# Import new samples
cat("\nImporting", length(new_samples), "new samples...\n")

success_count <- 0
error_count <- 0
errors <- list()

for (i in seq_along(new_samples)) {
  sample_id <- new_samples[i]
  
  tryCatch({
    # Get sequence
    sequence <- sequences[[sample_id]]
    
    # Get metadata for this sample
    meta_row <- metadata[metadata$sample_id == sample_id, ]
    
    # Add sequence
    add_sequence(con, sample_id, sequence)
    
    # Build notes field with extra info from BioNumerics
    notes_parts <- c()
    if (nrow(meta_row) > 0) {
      if (!is.na(meta_row$genotype) && meta_row$genotype != "") {
        notes_parts <- c(notes_parts, paste0("Genotype: ", meta_row$genotype))
      }
      if (!is.na(meta_row$outbreak_variant) && meta_row$outbreak_variant != "") {
        notes_parts <- c(notes_parts, paste0("Variant: ", meta_row$outbreak_variant))
      }
      if (!is.na(meta_row$patient_initials) && meta_row$patient_initials != "") {
        notes_parts <- c(notes_parts, paste0("Patient: ", meta_row$patient_initials))
      }
    }
    notes <- if (length(notes_parts) > 0) paste(notes_parts, collapse = "; ") else NA
    
    # Add metadata
    if (nrow(meta_row) > 0) {
      add_metadata(con, 
        sample_id = sample_id,
        sampling_date = if (!is.na(meta_row$sampling_date) && meta_row$sampling_date != "") meta_row$sampling_date else NA,
        geo_location = if (!is.na(meta_row$origin) && meta_row$origin != "") meta_row$origin else NA,
        geo_county = NA,
        geo_country = "Norway",
        virus_name = "Hepatitis A Virus",
        submitter = "BioNumerics Import",
        sequencing_run = if (!is.na(meta_row$outbreak_variant) && meta_row$outbreak_variant != "") meta_row$outbreak_variant else NA,
        notes = notes
      )
    } else {
      add_metadata(con,
        sample_id = sample_id,
        virus_name = "Hepatitis A Virus",
        submitter = "BioNumerics Import"
      )
    }
    
    success_count <- success_count + 1
    
    if (i %% 50 == 0) {
      cat("  Processed", i, "of", length(new_samples), "\n")
    }
    
  }, error = function(e) {
    error_count <<- error_count + 1
    errors[[sample_id]] <<- e$message
  })
}

# Final summary
cat("\n" , rep("=", 50), "\n", sep = "")
cat("IMPORT COMPLETE\n")
cat(rep("=", 50), "\n", sep = "")
cat("Successfully imported:", success_count, "samples\n")
cat("Errors:", error_count, "\n")

if (error_count > 0) {
  cat("\nErrors encountered:\n")
  for (sid in names(errors)) {
    cat("  ", sid, ":", errors[[sid]], "\n")
  }
}

# Show database stats
total_seqs <- dbGetQuery(con, "SELECT COUNT(*) as n FROM sequences")$n
cat("\nTotal sequences in database:", total_seqs, "\n")

# Show recent imports
cat("\nRecently imported samples (first 10):\n")
recent <- dbGetQuery(con, "
  SELECT s.sample_id, s.sequence_length, m.sampling_date, m.geo_location, m.notes
  FROM sequences s
  LEFT JOIN metadata m ON s.sample_id = m.sample_id
  ORDER BY s.created_at DESC
  LIMIT 10
")
print(recent)

# Close connection
close_db_connection(con)
cat("\nDone!\n")
