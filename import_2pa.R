# =============================================================================
# Import HAV data from 2PA.fa and export.csv
# =============================================================================

library(DBI)
library(RSQLite)

# Paths
fasta_file <- "export/2PA.fa"
csv_file <- "export/export.csv"
db_path <- "data/hav.db"

# -----------------------------------------------------------------------------
# Step 1: Clear the database
# -----------------------------------------------------------------------------
cat("Clearing database...\n")
con <- dbConnect(RSQLite::SQLite(), db_path)
dbExecute(con, "DELETE FROM metadata")
dbExecute(con, "DELETE FROM sequences")
dbExecute(con, "DELETE FROM analysis_results")
cat("  Database cleared.\n")

# -----------------------------------------------------------------------------
# Step 2: Read and parse FASTA file
# -----------------------------------------------------------------------------
cat("\nReading FASTA file...\n")
fasta_lines <- readLines(fasta_file, encoding = "latin1", warn = FALSE)

sequences <- list()
current_id <- NULL
current_genotype <- NULL
current_variant <- NULL
current_seq <- ""

for (line in fasta_lines) {
  if (startsWith(line, ">")) {
    # Save previous sequence
    if (!is.null(current_id) && nchar(current_seq) > 0) {
      sequences[[current_id]] <- list(
        sequence = toupper(current_seq),
        genotype = current_genotype,
        variant = current_variant
      )
    }
    
    # Parse header: >2PA|sample_id|genotype|variant
    header <- sub("^>", "", line)
    parts <- strsplit(header, "\\|")[[1]]
    
    # parts[1] = "2PA" (not needed)
    # parts[2] = sample_id
    # parts[3] = genotype
    # parts[4] = variant (may be missing)
    current_id <- if (length(parts) >= 2) parts[2] else NA
    current_genotype <- if (length(parts) >= 3) parts[3] else NA
    current_variant <- if (length(parts) >= 4) parts[4] else NA
    current_seq <- ""
    
  } else {
    # Append sequence line
    current_seq <- paste0(current_seq, gsub("\\s", "", line))
  }
}

# Don't forget the last sequence
if (!is.null(current_id) && nchar(current_seq) > 0) {
  sequences[[current_id]] <- list(
    sequence = toupper(current_seq),
    genotype = current_genotype,
    variant = current_variant
  )
}

cat("  Found", length(sequences), "sequences in FASTA\n")

# -----------------------------------------------------------------------------
# Step 3: Read CSV metadata
# -----------------------------------------------------------------------------
cat("\nReading CSV metadata...\n")
metadata <- read.csv(csv_file, sep = ";", stringsAsFactors = FALSE, 
                     fileEncoding = "latin1", na.strings = c("", "NA"))

# Clean column names
names(metadata) <- trimws(names(metadata))
cat("  Found", nrow(metadata), "rows in CSV\n")
cat("  Columns:", paste(names(metadata), collapse = ", "), "\n")

# Create lookup by Key (sample_id)
metadata$Key <- as.character(metadata$Key)

# -----------------------------------------------------------------------------
# Step 4: Helper function to convert date format
# -----------------------------------------------------------------------------
convert_date <- function(date_str) {
  if (is.na(date_str) || date_str == "") return(NA)
  # Try DD.MM.YYYY format
  if (grepl("^\\d{2}\\.\\d{2}\\.\\d{4}$", date_str)) {
    parts <- strsplit(date_str, "\\.")[[1]]
    return(paste(parts[3], parts[2], parts[1], sep = "-"))
  }
  return(date_str)
}

# -----------------------------------------------------------------------------
# Step 5: Import data into database
# -----------------------------------------------------------------------------
cat("\nImporting into database...\n")

dbExecute(con, "PRAGMA foreign_keys = ON;")

success <- 0
errors <- 0

for (sample_id in names(sequences)) {
  tryCatch({
    seq_data <- sequences[[sample_id]]
    
    # Insert sequence
    dbExecute(con, 
      "INSERT INTO sequences (sample_id, sequence, sequence_length) VALUES (?, ?, ?)",
      params = list(sample_id, seq_data$sequence, nchar(seq_data$sequence))
    )
    
    # Find metadata row
    meta_row <- metadata[metadata$Key == sample_id, ]
    
    # Prepare metadata values
    sampling_date <- NA
    origin <- NA
    source <- NA
    year <- NA
    genotype <- seq_data$genotype
    variant <- seq_data$variant
    patient <- NA
    
    if (nrow(meta_row) > 0) {
      meta_row <- meta_row[1, ]  # Take first match
      sampling_date <- convert_date(meta_row$Sample.date)
      origin <- if (!is.na(meta_row$Origin)) meta_row$Origin else NA
      source <- if (!is.na(meta_row$Source)) meta_row$Source else NA
      year <- if (!is.na(meta_row$Year)) as.character(meta_row$Year) else NA
      patient <- if (!is.na(meta_row$Patient.initials)) meta_row$Patient.initials else NA
      # Use CSV genotype/variant if available (may be more accurate)
      if (!is.na(meta_row$Genotype)) genotype <- meta_row$Genotype
      if (!is.na(meta_row$OUTBREAK_VARIANT)) variant <- meta_row$OUTBREAK_VARIANT
    }
    
    # Build notes
    notes_parts <- c()
    if (!is.na(genotype)) notes_parts <- c(notes_parts, paste0("Genotype: ", genotype))
    if (!is.na(variant)) notes_parts <- c(notes_parts, paste0("Variant: ", variant))
    if (!is.na(patient)) notes_parts <- c(notes_parts, paste0("Patient: ", patient))
    if (!is.na(year)) notes_parts <- c(notes_parts, paste0("Year: ", year))
    notes <- if (length(notes_parts) > 0) paste(notes_parts, collapse = "; ") else NA
    
    # Insert metadata
    dbExecute(con, "
      INSERT INTO metadata 
      (sample_id, sampling_date, geo_location, geo_county, geo_country, 
       virus_name, submitter, sequencing_run, notes)
      VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
      params = list(
        sample_id,
        sampling_date,
        origin,
        NA,  # geo_county
        "Norway",
        "Hepatitis A Virus",
        "BioNumerics Export",
        variant,  # Use variant as sequencing_run for grouping
        notes
      )
    )
    
    success <- success + 1
    if (success %% 100 == 0) cat("  Imported", success, "samples...\n")
    
  }, error = function(e) {
    errors <<- errors + 1
    cat("  Error with", sample_id, ":", e$message, "\n")
  })
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n", rep("=", 50), "\n", sep = "")
cat("IMPORT COMPLETE\n")
cat(rep("=", 50), "\n", sep = "")
cat("Successfully imported:", success, "samples\n")
cat("Errors:", errors, "\n")

# Verify
total <- dbGetQuery(con, "SELECT COUNT(*) as n FROM sequences")$n
cat("\nTotal sequences in database:", total, "\n")

# Show sample data
cat("\nSample entries:\n")
sample_data <- dbGetQuery(con, "
  SELECT s.sample_id, s.sequence_length, m.sampling_date, m.geo_location, m.sequencing_run, m.notes
  FROM sequences s
  LEFT JOIN metadata m ON s.sample_id = m.sample_id
  LIMIT 10
")
print(sample_data)

# Genotype distribution from notes
cat("\nGenotype distribution:\n")
notes_data <- dbGetQuery(con, "SELECT notes FROM metadata WHERE notes IS NOT NULL")
genotypes <- sapply(notes_data$notes, function(n) {
  m <- regmatches(n, regexpr("Genotype: [^;]+", n))
  if (length(m) > 0) sub("Genotype: ", "", m) else NA
})
print(table(genotypes, useNA = "ifany"))

dbDisconnect(con)
cat("\nDone!\n")
