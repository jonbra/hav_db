# =============================================================================
# Sequence Functions for HAV Database
# =============================================================================
# Functions for importing/exporting sequences in FASTA format and handling
# DNA sequences with Biostrings.

library(seqinr)
library(Biostrings)

# =============================================================================
# FASTA Import Functions
# =============================================================================

#' Read sequences from a FASTA file
#' @param fasta_file Path to FASTA file
#' @return Data frame with sample_id and sequence columns
read_fasta_to_df <- function(fasta_file) {
  if (!file.exists(fasta_file)) {
    stop("FASTA file not found: ", fasta_file)
  }
  
  # Read using seqinr
  seqs <- seqinr::read.fasta(fasta_file, as.string = TRUE, forceDNAtolower = FALSE)
  
  data.frame(
    sample_id = names(seqs),
    sequence = sapply(seqs, function(x) toupper(as.character(x))),
    stringsAsFactors = FALSE
  )
}

#' Import FASTA file into database
#' @param con Database connection
#' @param fasta_file Path to FASTA file
#' @param metadata_df Optional data frame with metadata (must have sample_id column)
#' @param skip_existing Skip samples that already exist in database
#' @return List with counts of imported/skipped sequences
import_fasta_to_db <- function(con, fasta_file, metadata_df = NULL, skip_existing = TRUE) {
  # Read sequences
  seqs_df <- read_fasta_to_df(fasta_file)
  
  imported <- 0
  skipped <- 0
  errors <- character()
  
  # Check for existing sequences
  existing <- dbGetQuery(con, "SELECT sample_id FROM sequences")$sample_id
  
  for (i in seq_len(nrow(seqs_df))) {
    sample_id <- seqs_df$sample_id[i]
    sequence <- seqs_df$sequence[i]
    
    # Skip if exists
    if (skip_existing && sample_id %in% existing) {
      skipped <- skipped + 1
      next
    }
    
    tryCatch({
      # Add sequence
      add_sequence(con, sample_id, sequence)
      
      # Add metadata if provided
      if (!is.null(metadata_df) && sample_id %in% metadata_df$sample_id) {
        meta_row <- metadata_df[metadata_df$sample_id == sample_id, , drop = FALSE]
        
        add_metadata(con, sample_id,
          sampling_date = meta_row$sampling_date[1] %||% NULL,
          geo_location = meta_row$geo_location[1] %||% NULL,
          geo_county = meta_row$geo_county[1] %||% NULL,
          geo_country = meta_row$geo_country[1] %||% "Norway",
          virus_name = meta_row$virus_name[1] %||% NULL,
          submitter = meta_row$submitter[1] %||% NULL,
          sequencing_run = meta_row$sequencing_run[1] %||% NULL,
          notes = meta_row$notes[1] %||% NULL
        )
      }
      
      imported <- imported + 1
    }, error = function(e) {
      errors <<- c(errors, paste0(sample_id, ": ", e$message))
    })
  }
  
  list(
    imported = imported,
    skipped = skipped,
    errors = errors,
    total_in_file = nrow(seqs_df)
  )
}

#' Import sequences with metadata - ONLY imports sequences listed in metadata
#' 
#' This is the recommended function for importing sequences. It takes a multi-FASTA
#' file and a metadata file (CSV or data frame). Only sequences that have a matching
#' entry in the metadata file will be imported. The FASTA headers must match the
#' sample_id column in the metadata.
#' 
#' @param con Database connection
#' @param fasta_file Path to multi-FASTA file
#' @param metadata_file Path to metadata CSV file, or a data frame
#' @param skip_existing Skip samples that already exist in database (default TRUE)
#' @param verbose Print progress messages (default TRUE)
#' @return List with import summary
#' 
#' @details
#' Metadata file format (CSV):
#' - Required column: sample_id (must match FASTA headers)
#' - Optional columns: sampling_date, sequence_length, geo_location, geo_county,
#'   geo_country, organism, gene_target, submitter, sequencing_run, notes
#' 
#' The sequence_length column in metadata is informational only - actual length
#' is calculated from the sequence.
#' 
#' @examples
#' \\dontrun{
#' con <- get_db_connection("data/hav.db")
#' result <- import_with_metadata(con, "sequences.fasta", "metadata.csv")
#' print(result)
#' }
import_with_metadata <- function(con, fasta_file, metadata_file, 
                                  skip_existing = TRUE, verbose = TRUE) {
  
  # Load metadata
  if (is.character(metadata_file)) {
    if (!file.exists(metadata_file)) {
      stop("Metadata file not found: ", metadata_file)
    }
    metadata_df <- read.csv(metadata_file, stringsAsFactors = FALSE)
  } else if (is.data.frame(metadata_file)) {
    metadata_df <- metadata_file
  } else {
    stop("metadata_file must be a file path or data frame")
  }
  
  # Validate metadata has sample_id column
  if (!"sample_id" %in% names(metadata_df)) {
    stop("Metadata must contain 'sample_id' column")
  }
  
  # Ensure sample_id is character
  metadata_df$sample_id <- as.character(metadata_df$sample_id)
  
  # Check for duplicate sample_ids in metadata
  dup_ids <- metadata_df$sample_id[duplicated(metadata_df$sample_id)]
  if (length(dup_ids) > 0) {
    warning("Duplicate sample_ids in metadata (using first occurrence): ", 
            paste(head(dup_ids, 5), collapse = ", "),
            if (length(dup_ids) > 5) paste0(" ... and ", length(dup_ids) - 5, " more"))
    metadata_df <- metadata_df[!duplicated(metadata_df$sample_id), ]
  }
  
  # Read FASTA file
  if (verbose) cat("Reading FASTA file...\n")
  seqs_df <- read_fasta_to_df(fasta_file)
  
  if (verbose) {
    cat("  Found", nrow(seqs_df), "sequences in FASTA file\n")
    cat("  Found", nrow(metadata_df), "entries in metadata file\n")
  }
  
  # Find sequences that are in metadata
  samples_to_import <- intersect(metadata_df$sample_id, seqs_df$sample_id)
  samples_in_meta_not_fasta <- setdiff(metadata_df$sample_id, seqs_df$sample_id)
  samples_in_fasta_not_meta <- setdiff(seqs_df$sample_id, metadata_df$sample_id)
  
  if (verbose) {
    cat("  Sequences to import:", length(samples_to_import), "\n")
    if (length(samples_in_meta_not_fasta) > 0) {
      cat("  WARNING: Samples in metadata but not in FASTA:", 
          length(samples_in_meta_not_fasta), "\n")
    }
    cat("  Sequences in FASTA but not in metadata (will be skipped):", 
        length(samples_in_fasta_not_meta), "\n")
  }
  
  if (length(samples_to_import) == 0) {
    warning("No matching samples found between FASTA and metadata!")
    return(list(
      imported = 0,
      skipped_existing = 0,
      skipped_no_metadata = length(seqs_df$sample_id),
      not_in_fasta = samples_in_meta_not_fasta,
      errors = character(),
      total_in_fasta = nrow(seqs_df),
      total_in_metadata = nrow(metadata_df)
    ))
  }
  
  # Check for existing sequences
  existing <- dbGetQuery(con, "SELECT sample_id FROM sequences")$sample_id
  
  imported <- 0
  skipped_existing <- 0
  errors <- character()
  
  if (verbose) cat("\nImporting sequences...\n")
  
  for (sample_id in samples_to_import) {
    # Skip if exists
    if (skip_existing && sample_id %in% existing) {
      skipped_existing <- skipped_existing + 1
      next
    }
    
    # Get sequence
    sequence <- seqs_df$sequence[seqs_df$sample_id == sample_id]
    
    # Get metadata row
    meta_row <- metadata_df[metadata_df$sample_id == sample_id, , drop = FALSE]
    
    tryCatch({
      # Start transaction
      dbBegin(con)
      
      # Add sequence
      add_sequence(con, sample_id, sequence)
      
      # Add metadata (with safe column access)
      add_metadata(con, sample_id,
        sampling_date = safe_get(meta_row, "sampling_date"),
        geo_location = safe_get(meta_row, "geo_location"),
        geo_county = safe_get(meta_row, "geo_county"),
        geo_country = safe_get(meta_row, "geo_country", "Norway"),
        virus_name = safe_get(meta_row, "virus_name"),
        submitter = safe_get(meta_row, "submitter"),
        sequencing_run = safe_get(meta_row, "sequencing_run"),
        notes = safe_get(meta_row, "notes")
      )
      
      dbCommit(con)
      imported <- imported + 1
      
      if (verbose && imported %% 10 == 0) {
        cat("  Imported", imported, "sequences...\n")
      }
      
    }, error = function(e) {
      dbRollback(con)
      errors <<- c(errors, paste0(sample_id, ": ", e$message))
    })
  }
  
  if (verbose) {
    cat("\n=== Import Summary ===\n")
    cat("Imported:", imported, "\n")
    cat("Skipped (already exist):", skipped_existing, "\n")
    cat("Skipped (no metadata):", length(samples_in_fasta_not_meta), "\n")
    if (length(errors) > 0) {
      cat("Errors:", length(errors), "\n")
      for (e in head(errors, 5)) cat("  -", e, "\n")
    }
  }
  
  list(
    imported = imported,
    skipped_existing = skipped_existing,
    skipped_no_metadata = length(samples_in_fasta_not_meta),
    not_in_fasta = samples_in_meta_not_fasta,
    errors = errors,
    total_in_fasta = nrow(seqs_df),
    total_in_metadata = nrow(metadata_df)
  )
}

#' Safely get a column value from a data frame row
#' @param df Data frame (single row)
#' @param col Column name
#' @param default Default value if column doesn't exist or value is NA/empty
#' @return Column value or default
safe_get <- function(df, col, default = NULL) {
  if (!col %in% names(df)) return(default)
  val <- df[[col]][1]
  if (is.null(val) || is.na(val) || val == "") return(default)
  return(val)
}

#' Create a template metadata data frame
#' 
#' Helper function to create an empty metadata template with all columns
#' @param n Number of rows to create (default 0)
#' @return Data frame with correct column structure
#' 
#' @examples
#' template <- create_metadata_template()
#' write.csv(template, "metadata_template.csv", row.names = FALSE)
create_metadata_template <- function(n = 0) {
  data.frame(
    sample_id = character(n),
    sampling_date = as.Date(character(n)),
    sequence_length = integer(n),
    geo_location = character(n),
    geo_county = character(n),
    geo_country = character(n),
    virus_name = character(n),
    submitter = character(n),
    sequencing_run = character(n),
    notes = character(n),
    stringsAsFactors = FALSE
  )
}

#' Validate metadata file before import
#' 
#' Check metadata file for common issues before importing
#' @param metadata_file Path to CSV file or data frame
#' @param fasta_file Optional path to FASTA file to check for matching IDs
#' @return List with validation results
#' 
#' @examples
#' validation <- validate_metadata("metadata.csv", "sequences.fasta")
#' if (length(validation$errors) == 0) {
#'   cat("Metadata is valid!\n")
#' }
validate_metadata <- function(metadata_file, fasta_file = NULL) {
  # Load metadata
  if (is.character(metadata_file)) {
    if (!file.exists(metadata_file)) {
      return(list(valid = FALSE, errors = "Metadata file not found"))
    }
    metadata_df <- read.csv(metadata_file, stringsAsFactors = FALSE)
  } else {
    metadata_df <- metadata_file
  }
  
  errors <- character()
  warnings <- character()
  
  # Check for sample_id column
  if (!"sample_id" %in% names(metadata_df)) {
    errors <- c(errors, "Missing required column: sample_id")
    return(list(valid = FALSE, errors = errors, warnings = warnings))
  }
  
  # Check for empty sample_ids
  empty_ids <- sum(is.na(metadata_df$sample_id) | metadata_df$sample_id == "")
  if (empty_ids > 0) {
    errors <- c(errors, paste("Found", empty_ids, "rows with empty sample_id"))
  }
  
  # Check for duplicates
  dup_ids <- metadata_df$sample_id[duplicated(metadata_df$sample_id)]
  if (length(dup_ids) > 0) {
    warnings <- c(warnings, paste("Found", length(dup_ids), "duplicate sample_ids"))
  }
  
  # Check date format if column exists - STRICT YYYY-MM-DD format required
  if ("sampling_date" %in% names(metadata_df)) {
    dates <- metadata_df$sampling_date[!is.na(metadata_df$sampling_date) & metadata_df$sampling_date != ""]
    if (length(dates) > 0) {
      # Check for strict YYYY-MM-DD format using regex
      date_pattern <- "^\\d{4}-\\d{2}-\\d{2}$"
      bad_format <- !grepl(date_pattern, dates)
      if (any(bad_format)) {
        errors <- c(errors, paste(sum(bad_format), "dates do not match required YYYY-MM-DD format:",
                                  paste(head(dates[bad_format], 3), collapse = ", "),
                                  if (sum(bad_format) > 3) "..."))
      } else {
        # Format is correct, check if dates are valid
        parsed <- suppressWarnings(as.Date(dates, format = "%Y-%m-%d"))
        bad_dates <- sum(is.na(parsed))
        if (bad_dates > 0) {
          errors <- c(errors, paste(bad_dates, "dates are invalid (must be valid dates in YYYY-MM-DD format)"))
        }
      }
    }
  }
  
  # Check against FASTA if provided
  if (!is.null(fasta_file)) {
    if (file.exists(fasta_file)) {
      seqs_df <- read_fasta_to_df(fasta_file)
      matching <- sum(metadata_df$sample_id %in% seqs_df$sample_id)
      not_in_fasta <- sum(!metadata_df$sample_id %in% seqs_df$sample_id)
      
      if (matching == 0) {
        errors <- c(errors, "No sample_ids match between metadata and FASTA")
      } else if (not_in_fasta > 0) {
        warnings <- c(warnings, paste(not_in_fasta, "sample_ids in metadata not found in FASTA"))
      }
    }
  }
  
  # List recognized columns
  valid_cols <- c("sample_id", "sampling_date", "sequence_length", "geo_location", 
                  "geo_county", "geo_country", "virus_name", 
                  "submitter", "sequencing_run", "notes")
  extra_cols <- setdiff(names(metadata_df), valid_cols)
  if (length(extra_cols) > 0) {
    warnings <- c(warnings, paste("Extra columns will be ignored:", paste(extra_cols, collapse = ", ")))
  }
  
  list(
    valid = length(errors) == 0,
    errors = errors,
    warnings = warnings,
    n_samples = nrow(metadata_df),
    columns = names(metadata_df)
  )
}

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x) || is.na(x) || x == "") y else x

#' Import metadata from CSV file
#' @param con Database connection
#' @param csv_file Path to CSV file (must have sample_id column)
#' @param update_existing Update metadata for existing samples
#' @return Number of records updated
import_metadata_csv <- function(con, csv_file, update_existing = TRUE) {
  if (!file.exists(csv_file)) {
    stop("CSV file not found: ", csv_file)
  }
  
  metadata_df <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  if (!"sample_id" %in% names(metadata_df)) {
    stop("CSV must contain 'sample_id' column")
  }
  
  updated <- 0
  
  for (i in seq_len(nrow(metadata_df))) {
    row <- metadata_df[i, , drop = FALSE]
    sample_id <- row$sample_id
    
    # Check if sample exists
    exists <- dbGetQuery(con, "SELECT COUNT(*) as n FROM sequences WHERE sample_id = ?",
                         params = list(sample_id))$n > 0
    
    if (!exists) {
      warning("Sample not found in database: ", sample_id)
      next
    }
    
    tryCatch({
      add_metadata(con, sample_id,
        sampling_date = row$sampling_date %||% NULL,
        geo_location = row$geo_location %||% NULL,
        geo_county = row$geo_county %||% NULL,
        geo_country = row$geo_country %||% "Norway",
        virus_name = row$virus_name %||% NULL,
        submitter = row$submitter %||% NULL,
        sequencing_run = row$sequencing_run %||% NULL,
        notes = row$notes %||% NULL
      )
      updated <- updated + 1
    }, error = function(e) {
      warning("Failed to update ", sample_id, ": ", e$message)
    })
  }
  
  updated
}

# =============================================================================
# FASTA Export Functions
# =============================================================================

#' Export sequences to FASTA file
#' @param con Database connection
#' @param output_file Path for output FASTA file
#' @param sample_ids Optional vector of sample IDs to export (NULL = all)
#' @param include_metadata Include metadata in FASTA headers
#' @return Number of sequences exported
export_to_fasta <- function(con, output_file, sample_ids = NULL, include_metadata = FALSE) {
  # Get sequences
  if (!is.null(sample_ids) && length(sample_ids) > 0) {
    data <- get_sequences_with_metadata(con, sample_ids)
  } else {
    data <- get_sequences_with_metadata(con)
  }
  
  if (nrow(data) == 0) {
    warning("No sequences to export")
    return(0)
  }
  
  # Build FASTA content
  fasta_lines <- character()
  
  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    
    # Build header
    if (include_metadata) {
      header_parts <- row$sample_id
      if (!is.na(row$organism) && row$organism != "") {
        header_parts <- paste0(header_parts, "|", row$organism)
      }
      if (!is.na(row$gene_target) && row$gene_target != "") {
        header_parts <- paste0(header_parts, "|", row$gene_target)
      }
      if (!is.na(row$sampling_date) && row$sampling_date != "") {
        header_parts <- paste0(header_parts, "|", row$sampling_date)
      }
      if (!is.na(row$geo_location) && row$geo_location != "") {
        header_parts <- paste0(header_parts, "|", row$geo_location)
      }
      header <- header_parts
    } else {
      header <- row$sample_id
    }
    
    fasta_lines <- c(fasta_lines, paste0(">", header))
    
    # Add sequence (wrapped at 80 characters)
    seq_wrapped <- gsub("(.{80})", "\\1\n", row$sequence)
    seq_wrapped <- sub("\n$", "", seq_wrapped)  # Remove trailing newline if any
    fasta_lines <- c(fasta_lines, seq_wrapped)
  }
  
  # Write to file
  writeLines(fasta_lines, output_file)
  
  nrow(data)
}

#' Export metadata to CSV
#' @param con Database connection
#' @param output_file Path for output CSV file
#' @param sample_ids Optional vector of sample IDs to export
#' @return Number of records exported
export_metadata_csv <- function(con, output_file, sample_ids = NULL) {
  if (!is.null(sample_ids) && length(sample_ids) > 0) {
    data <- get_sequences_with_metadata(con, sample_ids)
  } else {
    data <- get_sequences_with_metadata(con)
  }
  
  # Select metadata columns
  export_cols <- c("sample_id", "sequence_length", "sampling_date", "geo_location",
                   "geo_county", "geo_country", "organism", "gene_target", 
                   "submitter", "sequencing_run", "notes", "created_at")
  
  export_data <- data[, intersect(export_cols, names(data)), drop = FALSE]
  
  write.csv(export_data, output_file, row.names = FALSE)
  
  nrow(export_data)
}

# =============================================================================
# Biostrings Conversion Functions
# =============================================================================

#' Convert database sequences to DNAStringSet
#' @param con Database connection
#' @param sample_ids Vector of sample IDs
#' @return DNAStringSet object with named sequences
db_to_dna_stringset <- function(con, sample_ids) {
  if (length(sample_ids) == 0) {
    return(Biostrings::DNAStringSet())
  }
  
  data <- get_sequences_by_ids(con, sample_ids)
  
  if (nrow(data) == 0) {
    return(Biostrings::DNAStringSet())
  }
  
  seqs <- Biostrings::DNAStringSet(data$sequence)
  names(seqs) <- data$sample_id
  
  seqs
}

#' Convert DNAStringSet to data frame for database import
#' @param dna_set DNAStringSet object
#' @return Data frame with sample_id and sequence columns
dna_stringset_to_df <- function(dna_set) {
  data.frame(
    sample_id = names(dna_set),
    sequence = as.character(dna_set),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Sequence Validation Functions
# =============================================================================

#' Validate DNA sequence
#' @param sequence Character string
#' @return List with is_valid, length, and issues
validate_sequence <- function(sequence) {
  # Clean sequence
  seq_clean <- toupper(gsub("[^A-Za-z]", "", sequence))
  
  # Check for valid characters (IUPAC nucleotide codes)
  valid_chars <- c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", 
                   "B", "D", "H", "V", "N", "-")
  seq_chars <- strsplit(seq_clean, "")[[1]]
  invalid_chars <- setdiff(unique(seq_chars), valid_chars)
  
  issues <- character()
  
  if (length(invalid_chars) > 0) {
    issues <- c(issues, paste("Invalid characters:", paste(invalid_chars, collapse = ", ")))
  }
  
  if (nchar(seq_clean) == 0) {
    issues <- c(issues, "Empty sequence")
  }
  
  if (nchar(seq_clean) < 50) {
    issues <- c(issues, "Sequence very short (< 50 bp)")
  }
  
  # Check for high N content
  n_count <- sum(seq_chars == "N")
  n_fraction <- n_count / length(seq_chars)
  if (n_fraction > 0.1) {
    issues <- c(issues, sprintf("High N content: %.1f%%", n_fraction * 100))
  }
  
  list(
    is_valid = length(invalid_chars) == 0 && nchar(seq_clean) > 0,
    length = nchar(seq_clean),
    gc_content = sum(seq_chars %in% c("G", "C")) / length(seq_chars),
    n_content = n_fraction,
    issues = issues
  )
}

#' Validate multiple sequences
#' @param sequences Named vector or list of sequences
#' @return Data frame with validation results
validate_sequences <- function(sequences) {
  results <- lapply(names(sequences), function(name) {
    val <- validate_sequence(sequences[[name]])
    data.frame(
      sample_id = name,
      is_valid = val$is_valid,
      length = val$length,
      gc_content = round(val$gc_content, 3),
      n_content = round(val$n_content, 3),
      issues = paste(val$issues, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, results)
}
