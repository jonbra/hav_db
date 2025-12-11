# =============================================================================
# Extract sequences and metadata from BioNumerics export
# =============================================================================

library(xml2)
library(base64enc)
library(dplyr)
library(tidyr)

# Set paths
export_dir <- "export"
entries_dir <- file.path(export_dir, "_Entries")
experiments_dir <- file.path(export_dir, "_Experiments")
output_fasta <- "data/bionumerics_sequences.fasta"
output_metadata <- "data/bionumerics_metadata.csv"

# Helper function to decode Base64
decode_b64 <- function(x) {

  if (is.na(x) || x == "") return("")
  tryCatch({
    rawToChar(base64decode(x))
  }, error = function(e) x)
}

# -----------------------------------------------------------------------------
# Parse all Entry files to get metadata
# -----------------------------------------------------------------------------
cat("Parsing metadata from _Entries...\n")

entry_files <- list.files(entries_dir, pattern = "\\.xml$", full.names = TRUE)
all_entries <- list()

for (ef in entry_files) {
  cat("  Reading:", basename(ef), "\n")
  doc <- read_xml(ef)
  entries <- xml_find_all(doc, ".//Entry")
  
  for (entry in entries) {
    # Get sample key (ID)
    key_node <- xml_find_first(entry, ".//Key")
    key_b64 <- xml_attr(key_node, "B64")
    sample_id <- if (!is.na(key_b64) && key_b64 == "1") {
      decode_b64(xml_text(key_node))
    } else {
      xml_text(key_node)
    }
    
    # Get experiment path (links to sequence file)
    exp_path_node <- xml_find_first(entry, ".//SEQ/Path")
    exp_path <- ""
    if (!is.na(xml_text(exp_path_node))) {
      path_b64 <- xml_attr(exp_path_node, "B64")
      exp_path <- if (!is.na(path_b64) && path_b64 == "1") {
        decode_b64(xml_text(exp_path_node))
      } else {
        xml_text(exp_path_node)
      }
    }
    
    # Get all fields
    fields <- xml_find_all(entry, ".//Field")
    field_data <- list(sample_id = sample_id, exp_path = exp_path)
    
    for (field in fields) {
      name_node <- xml_find_first(field, ".//Name")
      value_node <- xml_find_first(field, ".//Value")
      
      # Decode field name
      name_b64 <- xml_attr(name_node, "B64")
      field_name <- if (!is.na(name_b64) && name_b64 == "1") {
        decode_b64(xml_text(name_node))
      } else {
        xml_text(name_node)
      }
      
      # Decode field value
      value_b64 <- xml_attr(value_node, "B64")
      field_value <- if (!is.na(value_b64) && value_b64 == "1") {
        decode_b64(xml_text(value_node))
      } else {
        xml_text(value_node)
      }
      
      field_data[[field_name]] <- field_value
    }
    
    all_entries[[length(all_entries) + 1]] <- field_data
  }
}

# Convert to data frame
metadata_df <- bind_rows(all_entries)
cat("Found", nrow(metadata_df), "entries\n")

# -----------------------------------------------------------------------------
# Parse Experiment files to get sequences
# -----------------------------------------------------------------------------
cat("\nParsing sequences from _Experiments...\n")

sequences <- list()

for (i in 1:nrow(metadata_df)) {
  exp_path <- metadata_df$exp_path[i]
  sample_id <- metadata_df$sample_id[i]
  
  if (is.na(exp_path) || exp_path == "") next
  
  # Convert Windows path to current OS path
  exp_file <- file.path(export_dir, gsub("\\\\", "/", exp_path))
  
  if (file.exists(exp_file)) {
    doc <- read_xml(exp_file)
    seq_node <- xml_find_first(doc, ".//SequenceData")
    
    if (!is.na(seq_node)) {
      seq_data <- xml_text(seq_node)
      if (!is.na(seq_data) && nchar(seq_data) > 0) {
        sequences[[sample_id]] <- toupper(seq_data)
      }
    }
  }
  
  if (i %% 50 == 0) cat("  Processed", i, "of", nrow(metadata_df), "\n")
}

cat("Found", length(sequences), "sequences\n")

# -----------------------------------------------------------------------------
# Create output metadata (only for samples with sequences)
# -----------------------------------------------------------------------------
cat("\nCreating output files...\n")

# Filter metadata to only include samples with sequences
metadata_output <- metadata_df %>%
  filter(sample_id %in% names(sequences)) %>%
  select(
    sample_id,
    sampling_date = SAMPLE_DATE,
    sample_year = SAMPLE_YEAR,
    genotype = GENOTYPE,
    outbreak_variant = OUTBREAK_VARIANT,
    origin = ORIGIN,
    source = SOURCE,
    transmission_route = TRANSMISSION_ROUTE,
    patient_initials = PATIENT_INITIALS,
    group = GROUP
  ) %>%
  mutate(
    # Convert date format from DD.MM.YYYY to YYYY-MM-DD
    sampling_date = case_when(
      grepl("^\\d{2}\\.\\d{2}\\.\\d{4}$", sampling_date) ~ {
        parts <- strsplit(sampling_date, "\\.")[[1]]
        paste(parts[3], parts[2], parts[1], sep = "-")
      },
      TRUE ~ sampling_date
    )
  )

# Fix date conversion (need to do row-wise)
metadata_output$sampling_date <- sapply(metadata_output$sampling_date, function(d) {
  if (is.na(d) || d == "") return(NA)
  if (grepl("^\\d{2}\\.\\d{2}\\.\\d{4}$", d)) {
    parts <- strsplit(d, "\\.")[[1]]
    return(paste(parts[3], parts[2], parts[1], sep = "-"))
  }
  return(d)
})

# Create data directory if needed
dir.create("data", showWarnings = FALSE)

# Write metadata CSV
write.csv(metadata_output, output_metadata, row.names = FALSE, na = "")
cat("Wrote metadata to:", output_metadata, "\n")

# -----------------------------------------------------------------------------
# Write FASTA file
# -----------------------------------------------------------------------------
fasta_lines <- c()
for (sample_id in names(sequences)) {
  fasta_lines <- c(fasta_lines, paste0(">", sample_id))
  fasta_lines <- c(fasta_lines, sequences[[sample_id]])
}

writeLines(fasta_lines, output_fasta)
cat("Wrote", length(sequences), "sequences to:", output_fasta, "\n")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat("Total entries in database:", nrow(metadata_df), "\n")
cat("Sequences extracted:", length(sequences), "\n")
cat("Output files:\n")
cat("  - FASTA:", output_fasta, "\n
")
cat("  - Metadata:", output_metadata, "\n")

# Show genotype distribution
cat("\nGenotype distribution:\n")
print(table(metadata_output$genotype, useNA = "ifany"))

cat("\nYear distribution:\n")
print(table(metadata_output$sample_year, useNA = "ifany"))
