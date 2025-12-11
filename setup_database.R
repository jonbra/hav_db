# =============================================================================
# Database Initialization Script
# =============================================================================
# Run this script once to create the SQLite database with the required schema.
# Usage: source("setup_database.R") or Rscript setup_database.R

library(DBI)
library(RSQLite)

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Get the directory where this script is located
get_script_dir <- function() {
  # Try multiple methods to find script location
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  # Fallback to working directory
  return(getwd())
}

script_dir <- get_script_dir()

# Database path
db_path <- file.path(script_dir, "data", "hav.db")

cat("Creating HAV database...\n")
cat("Database path:", db_path, "\n\n")

# Create data directory if it doesn't exist
dir.create(dirname(db_path), showWarnings = FALSE, recursive = TRUE)

# Connect to database (creates file if it doesn't exist)
con <- dbConnect(RSQLite::SQLite(), db_path)

# Enable foreign keys
dbExecute(con, "PRAGMA foreign_keys = OFF;")

# Drop existing tables to recreate with new schema
cat("Dropping existing tables...\n")
tryCatch(dbExecute(con, "DROP TABLE IF EXISTS metadata"), error = function(e) NULL)
tryCatch(dbExecute(con, "DROP TABLE IF EXISTS sequences"), error = function(e) NULL)
tryCatch(dbExecute(con, "DROP TABLE IF EXISTS analysis_results"), error = function(e) NULL)

dbExecute(con, "PRAGMA foreign_keys = ON;")

# =============================================================================
# Create Tables
# =============================================================================

cat("Creating sequences table...\n")
dbExecute(con, "
CREATE TABLE IF NOT EXISTS sequences (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT UNIQUE NOT NULL,
    sequence TEXT NOT NULL,
    sequence_length INTEGER,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP
);
")

cat("Creating metadata table...\n")
dbExecute(con, "
CREATE TABLE IF NOT EXISTS metadata (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT UNIQUE NOT NULL,
    sampling_date DATE,
    sample_year INTEGER,
    genotype TEXT,
    variant TEXT,
    patient_id TEXT,
    geo_location TEXT,
    geo_country TEXT DEFAULT 'Norway',
    source TEXT,
    transmission_route TEXT,
    comment TEXT,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (sample_id) REFERENCES sequences(sample_id) ON DELETE CASCADE
);
")

cat("Creating analysis_results table...\n")
dbExecute(con, "
CREATE TABLE IF NOT EXISTS analysis_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    analysis_name TEXT NOT NULL,
    analysis_type TEXT NOT NULL,
    sample_ids TEXT NOT NULL,
    result_data TEXT,
    parameters TEXT,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP
);
")

# =============================================================================
# Create Indexes for searchable fields
# =============================================================================

cat("Creating indexes...\n")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_sequences_sample_id ON sequences(sample_id);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_metadata_sample_id ON metadata(sample_id);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_metadata_sampling_date ON metadata(sampling_date);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_metadata_sample_year ON metadata(sample_year);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_metadata_genotype ON metadata(genotype);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_metadata_variant ON metadata(variant);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_metadata_patient_id ON metadata(patient_id);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_metadata_geo_country ON metadata(geo_country);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_metadata_geo_location ON metadata(geo_location);")

# =============================================================================
# Verify Setup
# =============================================================================

cat("\nVerifying database setup...\n")
tables <- dbListTables(con)
cat("Tables created:", paste(tables, collapse = ", "), "\n")

for (table in tables) {
  fields <- dbListFields(con, table)
  cat(sprintf("  %s: %s\n", table, paste(fields, collapse = ", ")))
}

# Close connection
dbDisconnect(con)

cat("\n✓ Database setup complete!\n")
cat("Database location:", db_path, "\n")
