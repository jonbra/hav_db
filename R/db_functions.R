# =============================================================================
# Database Functions for HAV Database
# =============================================================================
# CRUD operations for sequences, metadata, and analysis results.

library(DBI)
library(RSQLite)
library(dplyr)
library(dbplyr)

# =============================================================================
# Connection Management
# =============================================================================

#' Get database connection
#' @param db_path Path to SQLite database file
#' @return Database connection object
get_db_connection <- function(db_path = NULL) {
  if (is.null(db_path)) {
    # Default path relative to this file's location
    db_path <- file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "data", "hav.db")
  }
  
  if (!file.exists(db_path)) {
    stop("Database not found at: ", db_path, "\nRun setup_database.R first.")
  }
  
  con <- dbConnect(RSQLite::SQLite(), db_path)
  dbExecute(con, "PRAGMA foreign_keys = ON;")
  return(con)
}

#' Close database connection safely
#' @param con Database connection object
close_db_connection <- function(con) {
  if (!is.null(con) && dbIsValid(con)) {
    dbDisconnect(con)
  }
}

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
# Sequence Operations
# =============================================================================

#' Add a new sequence to the database
#' @param con Database connection
#' @param sample_id Unique sample identifier
#' @param sequence DNA sequence string
#' @return TRUE if successful
add_sequence <- function(con, sample_id, sequence) {
  sequence <- toupper(gsub("[^A-Za-z]", "", sequence))  # Clean sequence
  seq_length <- nchar(sequence)
  
  dbExecute(con, 
    "INSERT INTO sequences (sample_id, sequence, sequence_length) VALUES (?, ?, ?)",
    params = list(sample_id, sequence, seq_length)
  )
  
  return(TRUE)
}

#' Get sequence by sample ID
#' @param con Database connection
#' @param sample_id Sample identifier
#' @return Data frame with sequence info
get_sequence <- function(con, sample_id) {
  dbGetQuery(con,
    "SELECT * FROM sequences WHERE sample_id = ?",
    params = list(sample_id)
  )
}

#' Get all sequences
#' @param con Database connection
#' @return Data frame with all sequences
get_all_sequences <- function(con) {
  dbGetQuery(con, "SELECT * FROM sequences ORDER BY created_at DESC")
}

#' Get sequences by multiple sample IDs
#' @param con Database connection
#' @param sample_ids Vector of sample identifiers
#' @return Data frame with sequences
get_sequences_by_ids <- function(con, sample_ids) {
  if (length(sample_ids) == 0) return(data.frame())
  
  placeholders <- paste(rep("?", length(sample_ids)), collapse = ", ")
  query <- sprintf("SELECT * FROM sequences WHERE sample_id IN (%s)", placeholders)
  
  dbGetQuery(con, query, params = as.list(sample_ids))
}

#' Delete a sequence (cascades to metadata)
#' @param con Database connection
#' @param sample_id Sample identifier
#' @return Number of rows deleted
delete_sequence <- function(con, sample_id) {
  dbExecute(con, 
    "DELETE FROM sequences WHERE sample_id = ?",
    params = list(sample_id)
  )
}

#' Update a sequence
#' @param con Database connection
#' @param sample_id Sample identifier
#' @param new_sequence New sequence string
#' @return Number of rows updated
update_sequence <- function(con, sample_id, new_sequence) {
  new_sequence <- toupper(gsub("[^A-Za-z]", "", new_sequence))
  seq_length <- nchar(new_sequence)
  
  dbExecute(con,
    "UPDATE sequences SET sequence = ?, sequence_length = ? WHERE sample_id = ?",
    params = list(new_sequence, seq_length, sample_id)
  )
}

# =============================================================================
# Metadata Operations
# =============================================================================

#' Convert empty strings to NA for database storage
#' RSQLite requires NA (not NULL) to store NULL values in the database
#' @param x Value to check
#' @return NA_character_ if empty/NULL, otherwise the value as character
to_na_if_empty <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  if (is.na(x)) return(NA_character_)
  if (is.character(x) && x == "") return(NA_character_)
  return(as.character(x))
}

#' Add or update metadata for a sequence
#' @param con Database connection
#' @param sample_id Sample identifier (must exist in sequences table)
#' @param ... Metadata fields
#' @return TRUE if successful
add_metadata <- function(con, sample_id, sampling_date = NULL, sample_year = NULL,
                         genotype = NULL, variant = NULL, patient_id = NULL,
                         geo_location = NULL, geo_country = NULL,
                         source = NULL, transmission_route = NULL, comment = NULL) {
  
  # Convert empty strings/NULL to NA for RSQLite
  sampling_date <- to_na_if_empty(sampling_date)
  sample_year <- if (is.null(sample_year) || is.na(sample_year)) NA else as.integer(sample_year)
  genotype <- to_na_if_empty(genotype)
  variant <- to_na_if_empty(variant)
  patient_id <- to_na_if_empty(patient_id)
  geo_location <- to_na_if_empty(geo_location)
  geo_country <- to_na_if_empty(geo_country)
  source <- to_na_if_empty(source)
  transmission_route <- to_na_if_empty(transmission_route)
  comment <- to_na_if_empty(comment)
  
  dbExecute(con, "
    INSERT OR REPLACE INTO metadata 
    (sample_id, sampling_date, sample_year, genotype, variant, patient_id,
     geo_location, geo_country, source, transmission_route, comment)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
    params = list(sample_id, sampling_date, sample_year, genotype, variant, patient_id,
                  geo_location, geo_country, source, transmission_route, comment)
  )
  
  return(TRUE)
}

#' Get metadata for a sample
#' @param con Database connection
#' @param sample_id Sample identifier
#' @return Data frame with metadata
get_metadata <- function(con, sample_id) {
  dbGetQuery(con,
    "SELECT * FROM metadata WHERE sample_id = ?",
    params = list(sample_id)
  )
}

#' Get all metadata
#' @param con Database connection
#' @return Data frame with all metadata
get_all_metadata <- function(con) {
  dbGetQuery(con, "SELECT * FROM metadata")
}

#' Update specific metadata fields
#' @param con Database connection
#' @param sample_id Sample identifier
#' @param ... Named fields to update
#' @return Number of rows updated
update_metadata <- function(con, sample_id, ...) {
  updates <- list(...)
  if (length(updates) == 0) return(0)
  
  valid_fields <- c("sampling_date", "sample_year", "genotype", "variant", "patient_id",
                    "geo_location", "geo_country", "source", "transmission_route", "comment")
  updates <- updates[names(updates) %in% valid_fields]
  
  if (length(updates) == 0) return(0)
  
  set_clause <- paste(names(updates), "= ?", sep = " ", collapse = ", ")
  query <- sprintf("UPDATE metadata SET %s WHERE sample_id = ?", set_clause)
  
  dbExecute(con, query, params = c(unname(updates), list(sample_id)))
}

# =============================================================================
# Combined Operations
# =============================================================================

#' Add sequence with metadata in one operation
#' @param con Database connection
#' @param sample_id Unique sample identifier
#' @param sequence DNA sequence string
#' @param ... Metadata fields
#' @return TRUE if successful
add_sequence_with_metadata <- function(con, sample_id, sequence, ...) {
  # Start transaction

  dbBegin(con)
  
  tryCatch({
    add_sequence(con, sample_id, sequence)
    add_metadata(con, sample_id, ...)
    dbCommit(con)
    return(TRUE)
  }, error = function(e) {
    dbRollback(con)
    stop("Failed to add sequence: ", e$message)
  })
}

#' Get sequences joined with metadata
#' @param con Database connection
#' @param sample_ids Optional vector of sample IDs to filter
#' @return Data frame with sequences and metadata
get_sequences_with_metadata <- function(con, sample_ids = NULL) {
  query <- "
    SELECT s.sample_id, s.sequence, s.sequence_length, s.created_at,
           m.sampling_date, m.sample_year, m.genotype, m.variant, m.patient_id,
           m.geo_location, m.geo_country, m.source, m.transmission_route, m.comment
    FROM sequences s
    LEFT JOIN metadata m ON s.sample_id = m.sample_id
  "
  

  if (!is.null(sample_ids) && length(sample_ids) > 0) {
    placeholders <- paste(rep("?", length(sample_ids)), collapse = ", ")
    query <- paste0(query, sprintf(" WHERE s.sample_id IN (%s)", placeholders))
    dbGetQuery(con, query, params = as.list(sample_ids))
  } else {
    dbGetQuery(con, query)
  }
}

# =============================================================================
# Search and Filter Operations
# =============================================================================

#' Search sequences by various criteria
#' @param con Database connection
#' @param genotype Filter by genotype
#' @param variant Filter by variant
#' @param patient_id Filter by patient ID
#' @param geo_location Filter by location
#' @param geo_country Filter by country
#' @param year_from Start year (inclusive)
#' @param year_to End year (inclusive)
#' @param date_from Start date (inclusive)
#' @param date_to End date (inclusive)
#' @param min_length Minimum sequence length
#' @param max_length Maximum sequence length
#' @return Data frame with matching sequences and metadata
search_sequences <- function(con, genotype = NULL, variant = NULL,
                             patient_id = NULL, geo_location = NULL, 
                             geo_country = NULL, year_from = NULL, year_to = NULL,
                             date_from = NULL, date_to = NULL,
                             min_length = NULL, max_length = NULL) {
  
  query <- "
    SELECT s.sample_id, s.sequence, s.sequence_length, s.created_at,
           m.sampling_date, m.sample_year, m.genotype, m.variant, m.patient_id,
           m.geo_location, m.geo_country, m.source, m.transmission_route, m.comment
    FROM sequences s
    LEFT JOIN metadata m ON s.sample_id = m.sample_id
    WHERE 1=1
  "
  
 params <- list()
  
  # Exact match for genotype and variant (selected from dropdown)
  if (!is.null(genotype) && genotype != "") {
    query <- paste0(query, " AND m.genotype = ?")
    params <- c(params, list(genotype))
  }
  
  if (!is.null(variant) && variant != "") {
    query <- paste0(query, " AND m.variant = ?")
    params <- c(params, list(variant))
  }
  
  # Partial match for patient_id (user types in text box)
  if (!is.null(patient_id) && patient_id != "") {
    query <- paste0(query, " AND m.patient_id LIKE ?")
    params <- c(params, list(paste0("%", patient_id, "%")))
  }
  
  # Exact match for location/country (selected from dropdown)
  if (!is.null(geo_location) && geo_location != "") {
    query <- paste0(query, " AND m.geo_location = ?")
    params <- c(params, list(geo_location))
  }
  
  if (!is.null(geo_country) && geo_country != "") {
    query <- paste0(query, " AND m.geo_country = ?")
    params <- c(params, list(geo_country))
  }
  
  if (!is.null(year_from) && !is.na(year_from)) {
    query <- paste0(query, " AND m.sample_year >= ?")
    params <- c(params, list(as.integer(year_from)))
  }
  
  if (!is.null(year_to) && !is.na(year_to)) {
    query <- paste0(query, " AND m.sample_year <= ?")
    params <- c(params, list(as.integer(year_to)))
  }
  
  if (!is.null(date_from) && date_from != "") {
    query <- paste0(query, " AND m.sampling_date >= ?")
    params <- c(params, list(as.character(date_from)))
  }
  
  if (!is.null(date_to) && date_to != "") {
    query <- paste0(query, " AND m.sampling_date <= ?")
    params <- c(params, list(as.character(date_to)))
  }
  
  if (!is.null(min_length)) {
    query <- paste0(query, " AND s.sequence_length >= ?")
    params <- c(params, list(min_length))
  }
  
  if (!is.null(max_length)) {
    query <- paste0(query, " AND s.sequence_length <= ?")
    params <- c(params, list(max_length))
  }
  
  query <- paste0(query, " ORDER BY s.created_at DESC")
  
  if (length(params) > 0) {
    dbGetQuery(con, query, params = params)
  } else {
    dbGetQuery(con, query)
  }
}

#' Get unique values for a metadata field (for dropdowns)
#' @param con Database connection
#' @param field Field name
#' @return Vector of unique values
get_unique_values <- function(con, field) {
  valid_fields <- c("genotype", "variant", "patient_id", "geo_location", 
                    "geo_country", "source", "transmission_route", "sample_year")
  
  if (!field %in% valid_fields) {
    stop("Invalid field: ", field)
  }
  
  query <- sprintf("SELECT DISTINCT %s FROM metadata WHERE %s IS NOT NULL ORDER BY %s",
                   field, field, field)
  result <- dbGetQuery(con, query)
  result[[1]]
}

# =============================================================================
# Analysis Results Operations
# =============================================================================
#' Save analysis results
#' @param con Database connection
#' @param analysis_name Name for this analysis
#' @param analysis_type Type (alignment, phylogeny, clustering)
#' @param sample_ids Vector of sample IDs included
#' @param result_data Serialized result data (JSON or RDS base64)
#' @param parameters Optional parameters used
#' @return ID of inserted record
save_analysis_result <- function(con, analysis_name, analysis_type, sample_ids,
                                 result_data, parameters = NULL) {
  sample_ids_str <- paste(sample_ids, collapse = ",")
  
  dbExecute(con, "
    INSERT INTO analysis_results (analysis_name, analysis_type, sample_ids, result_data, parameters)
    VALUES (?, ?, ?, ?, ?)",
    params = list(analysis_name, analysis_type, sample_ids_str, result_data, parameters)
  )
  
  # Return the ID of the inserted record
  dbGetQuery(con, "SELECT last_insert_rowid() as id")$id
}

#' Get analysis results
#' @param con Database connection
#' @param analysis_type Optional filter by type
#' @return Data frame with analysis results
get_analysis_results <- function(con, analysis_type = NULL) {
  if (!is.null(analysis_type)) {
    dbGetQuery(con, 
      "SELECT * FROM analysis_results WHERE analysis_type = ? ORDER BY created_at DESC",
      params = list(analysis_type)
    )
  } else {
    dbGetQuery(con, "SELECT * FROM analysis_results ORDER BY created_at DESC")
  }
}

#' Delete analysis result
#' @param con Database connection
#' @param id Analysis result ID
#' @return Number of rows deleted
delete_analysis_result <- function(con, id) {
  dbExecute(con, "DELETE FROM analysis_results WHERE id = ?", params = list(id))
}

# =============================================================================
# BLAST Results Operations
# =============================================================================

#' Save BLAST results to database
#' @param con Database connection
#' @param hits Data frame with BLAST hits (from run_blastn or run_blast_search)
#' @return Number of rows inserted
save_blast_results <- function(con, hits) {
  if (is.null(hits) || nrow(hits) == 0) return(0)
  
  # Delete existing results for this query
  for (qid in unique(hits$query_id)) {
    dbExecute(con, "DELETE FROM blast_results WHERE query_sample_id = ?", params = list(qid))
  }
  
  # Insert new results
  n_inserted <- 0
  for (i in seq_len(nrow(hits))) {
    row <- hits[i, ]
    dbExecute(con, "
      INSERT INTO blast_results 
      (query_sample_id, hit_sample_id, identity_pct, alignment_length, mismatches,
       gap_opens, query_start, query_end, subject_start, subject_end, evalue, bit_score, snp_count)
      VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
      params = list(
        row$query_id,
        row$subject_id,
        row$identity_pct,
        row$alignment_length,
        row$mismatches %||% NA,
        row$gap_opens %||% NA,
        row$query_start %||% NA,
        row$query_end %||% NA,
        row$subject_start %||% NA,
        row$subject_end %||% NA,
        row$evalue %||% NA,
        row$bit_score %||% NA,
        row$snp_count %||% NA
      )
    )
    n_inserted <- n_inserted + 1
  }
  
  n_inserted
}

#' Get BLAST results for a query sample
#' @param con Database connection
#' @param query_sample_id Query sample ID
#' @param max_hits Maximum hits to return (NULL for all)
#' @return Data frame with BLAST hits
get_blast_results <- function(con, query_sample_id, max_hits = NULL) {
  query <- "SELECT * FROM blast_results WHERE query_sample_id = ? ORDER BY identity_pct DESC, snp_count ASC"
  if (!is.null(max_hits)) {
    query <- paste(query, "LIMIT", as.integer(max_hits))
  }
  dbGetQuery(con, query, params = list(query_sample_id))
}

#' Get all BLAST results with metadata
#' @param con Database connection
#' @param query_sample_id Query sample ID
#' @return Data frame with BLAST hits joined with hit sequence metadata
get_blast_results_with_metadata <- function(con, query_sample_id) {
  dbGetQuery(con, "
    SELECT b.*, m.genotype, m.geo_country, m.sampling_date, m.sample_year
    FROM blast_results b
    LEFT JOIN metadata m ON b.hit_sample_id = m.sample_id
    WHERE b.query_sample_id = ?
    ORDER BY b.identity_pct DESC, b.snp_count ASC",
    params = list(query_sample_id)
  )
}

#' Delete BLAST results for a query sample
#' @param con Database connection
#' @param query_sample_id Query sample ID
#' @return Number of rows deleted
delete_blast_results <- function(con, query_sample_id) {
  dbExecute(con, "DELETE FROM blast_results WHERE query_sample_id = ?", 
            params = list(query_sample_id))
}

#' Get samples with stored BLAST results
#' @param con Database connection
#' @return Data frame with query sample IDs and result counts
get_samples_with_blast_results <- function(con) {
  dbGetQuery(con, "
    SELECT query_sample_id, COUNT(*) as hit_count, MAX(created_at) as last_run
    FROM blast_results 
    GROUP BY query_sample_id 
    ORDER BY last_run DESC")
}

# =============================================================================
# Database Statistics
# =============================================================================

#' Get database statistics
#' @param con Database connection
#' @return List with counts and summaries
get_db_stats <- function(con) {
  list(
    total_sequences = dbGetQuery(con, "SELECT COUNT(*) as n FROM sequences")$n,
    total_with_metadata = dbGetQuery(con, "SELECT COUNT(*) as n FROM metadata")$n,
    organisms = dbGetQuery(con, "SELECT organism, COUNT(*) as n FROM metadata GROUP BY organism"),
    gene_targets = dbGetQuery(con, "SELECT gene_target, COUNT(*) as n FROM metadata GROUP BY gene_target"),
    date_range = dbGetQuery(con, "SELECT MIN(sampling_date) as min_date, MAX(sampling_date) as max_date FROM metadata"),
    length_stats = dbGetQuery(con, "SELECT MIN(sequence_length) as min_len, MAX(sequence_length) as max_len, AVG(sequence_length) as avg_len FROM sequences")
  )
}
