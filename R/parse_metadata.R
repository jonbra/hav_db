# =============================================================================
# Metadata Parsing and Cleaning Script
# =============================================================================
# This script parses and cleans metadata from BioNumerics export files.
# It harmonizes country names to English, handles uncertain values,
# and prepares clean data for database import.
#
# Input:  export/export.csv, export/2PA.fa
# Output: data/cleaned_metadata.csv, data/cleaned_sequences.csv
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)

# =============================================================================
# Country Name Mapping (Norwegian/Other -> English)
# =============================================================================

country_mapping <- c(
  # Norwegian to English
  "Norge" = "Norway",
  "Noreg" = "Norway",
  "Sverige" = "Sweden",
  "Danmark" = "Denmark",
  "Finland" = "Finland",
  "Island" = "Iceland",
  "Tyskland" = "Germany",
  "Frankrike" = "France",
  "Spania" = "Spain",
  "Italia" = "Italy",
  "Hellas" = "Greece",
  "Tyrkia" = "Turkey",
  "Polen" = "Poland",
  "Russland" = "Russia",
  "Ukraina" = "Ukraine",
  "Storbritannia" = "United Kingdom",
  "England" = "United Kingdom",
  "Skottland" = "Scotland",
  "Irland" = "Ireland",
  "Nederland" = "Netherlands",
  "Belgia" = "Belgium",
  "Østerrike" = "Austria",
  "Sveits" = "Switzerland",
  "Tsjekkia" = "Czech Republic",
  "Ungarn" = "Hungary",
  "Romania" = "Romania",
  "Bulgaria" = "Bulgaria",
  "Serbia" = "Serbia",
  "Kroatia" = "Croatia",
  "Bosnia" = "Bosnia and Herzegovina",
  "Albania" = "Albania",
  "Kosovo" = "Kosovo",
  "Montenegro" = "Montenegro",
  "Nord-Makedonia" = "North Macedonia",
  "Makedonia" = "North Macedonia",
  "Slovenia" = "Slovenia",
  "Slovakia" = "Slovakia",
  "Estland" = "Estonia",
  "Latvia" = "Latvia",
  "Litauen" = "Lithuania",
  "Hviterussland" = "Belarus",
  
  # Middle East
  "Syria" = "Syria",
  "Libanon" = "Lebanon",
  "Jordan" = "Jordan",
 "Irak" = "Iraq",
  "Iran" = "Iran",
  "Israel" = "Israel",
  "Palestina" = "Palestine",
  "Saudi-Arabia" = "Saudi Arabia",
  "De forente arabiske emirater" = "United Arab Emirates",
  "UAE" = "United Arab Emirates",
  "Kuwait" = "Kuwait",
  "Jemen" = "Yemen",
  "Oman" = "Oman",
  "Qatar" = "Qatar",
  "Bahrain" = "Bahrain",
  
  # Africa
  "Egypt" = "Egypt",
  "Egyp" = "Egypt",
  "Marokko" = "Morocco",
  "Tunisia" = "Tunisia",
  "Algerie" = "Algeria",
  "Libya" = "Libya",
  "Sudan" = "Sudan",
  "Etiopia" = "Ethiopia",
  "Eritrea" = "Eritrea",
  "Somalia" = "Somalia",
  "Kenya" = "Kenya",
  "Uganda" = "Uganda",
  "Tanzania" = "Tanzania",
  "Rwanda" = "Rwanda",
  "Burundi" = "Burundi",
  "Kongo" = "Democratic Republic of the Congo",
  "DR Kongo" = "Democratic Republic of the Congo",
  "Kamerun" = "Cameroon",
  "Nigeria" = "Nigeria",
  "Ghana" = "Ghana",
  "Elfenbenskysten" = "Ivory Coast",
  "Senegal" = "Senegal",
  "Mali" = "Mali",
  "Gambia" = "Gambia",
  "Mosambik" = "Mozambique",
  "Zimbabwe" = "Zimbabwe",
  "Zambia" = "Zambia",
  "Malawi" = "Malawi",
  "Madagaskar" = "Madagascar",
  
  # Asia
  "Kina" = "China",
  "Japan" = "Japan",
  "India" = "India",
  "Pakistan" = "Pakistan",
  "Bangladesh" = "Bangladesh",
  "Sri Lanka" = "Sri Lanka",
  "Nepal" = "Nepal",
  "Afghanistan" = "Afghanistan",
  "Myanmar" = "Myanmar",
  "Thailand" = "Thailand",
  "Vietnam" = "Vietnam",
  "Kambodsja" = "Cambodia",
  "Laos" = "Laos",
  "Malaysia" = "Malaysia",
  "Singapore" = "Singapore",
  "Indonesia" = "Indonesia",
  "Filippinene" = "Philippines",
  "Filipinene" = "Philippines",
  "Taiwan" = "Taiwan",
  "Mongolia" = "Mongolia",
  
  # Americas
  "USA" = "United States",
  "Amerikas forente stater" = "United States",
  "Canada" = "Canada",
  "Kanada" = "Canada",
  "Mexico" = "Mexico",
  "Brasil" = "Brazil",
  "Argentina" = "Argentina",
  "Chile" = "Chile",
  "Peru" = "Peru",
  "Colombia" = "Colombia",
  
  # Oceania
  "Australia" = "Australia",
  "New Zealand" = "New Zealand",
  
  # Common variations and typos
  "UK" = "United Kingdom",
  "Skotland" = "Scotland",
  "Marokko" = "Morocco",
  "Tyskland" = "Germany",
  "Afrika" = "Africa"
)

# =============================================================================
# Helper Functions
# =============================================================================

#' Clean and standardize country/origin field
#' @param origin_str Raw origin string from CSV
#' @return List with 'country', 'location', and 'comment'
clean_origin <- function(origin_str) {
  if (is.na(origin_str) || origin_str == "" || origin_str == "NA") {
    return(list(country = NA_character_, location = NA_character_, comment = NA_character_))
  }
  
  # Trim whitespace
  origin_str <- trimws(origin_str)
  
  # Check for uncertainty markers
  is_uncertain <- grepl("\\?", origin_str)
  
  # Remove question marks and clean
  clean_str <- gsub("\\?", "", origin_str)
  clean_str <- trimws(clean_str)
  
  # Skip empty after cleaning
  if (clean_str == "") {
    return(list(country = NA_character_, location = NA_character_, comment = NA_character_))
  }
  
  # Patterns that indicate Norwegian locations (not countries)
  norwegian_location_patterns <- c(
    "Frognerseteren", "rkontakt", "kontakt", "Oslo", "Bergen", 
    "Trondheim", "Stavanger", "Drammen", "Norsk"
  )
  
  is_norwegian_location <- any(sapply(norwegian_location_patterns, function(pat) {
    grepl(pat, clean_str, ignore.case = TRUE)
  }))
  
  if (is_norwegian_location) {
    # This is a location/note within Norway
    comment <- clean_str
    if (is_uncertain) comment <- paste0(comment, " (uncertain)")
    return(list(country = "Norway", location = clean_str, comment = comment))
  }
  
  # Try to map to English country name
  mapped_country <- NA_character_
  
  # Direct match
  if (clean_str %in% names(country_mapping)) {
    mapped_country <- country_mapping[[clean_str]]
  } else {
    # Case-insensitive match
    idx <- which(tolower(names(country_mapping)) == tolower(clean_str))
    if (length(idx) > 0) {
      mapped_country <- country_mapping[[idx[1]]]
    } else {
      # Check if the value itself is already an English country name
      if (clean_str %in% country_mapping) {
        mapped_country <- clean_str
      } else {
        # Keep original if no mapping found (might be a valid country)
        mapped_country <- clean_str
      }
    }
  }
  
  comment <- if (is_uncertain) "Country uncertain" else NA_character_
  
  return(list(country = mapped_country, location = NA_character_, comment = comment))
}

#' Convert date from DD.MM.YYYY to YYYY-MM-DD
convert_date <- function(date_str) {
  if (is.na(date_str) || date_str == "") return(NA_character_)
  
  # Try DD.MM.YYYY format
  if (grepl("^\\d{2}\\.\\d{2}\\.\\d{4}$", date_str)) {
    parts <- strsplit(date_str, "\\.")[[1]]
    return(paste(parts[3], parts[2], parts[1], sep = "-"))
  }
  
  # Try other common formats
  if (grepl("^\\d{4}-\\d{2}-\\d{2}$", date_str)) {
    return(date_str)  # Already ISO format
  }
  
  return(NA_character_)
}

#' Parse FASTA file and extract sequences with header info
parse_fasta <- function(fasta_file) {
  cat("Reading FASTA file:", fasta_file, "\n")
  
  lines <- readLines(fasta_file, encoding = "latin1", warn = FALSE)
  
  sequences <- list()
  current_id <- NULL
  current_genotype <- NULL
  current_variant <- NULL
  current_seq <- ""
  
  for (line in lines) {
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
      
      current_id <- if (length(parts) >= 2) trimws(parts[2]) else NA
      current_genotype <- if (length(parts) >= 3) trimws(parts[3]) else NA
      current_variant <- if (length(parts) >= 4) trimws(parts[4]) else NA
      current_seq <- ""
      
    } else {
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
  
  cat("  Parsed", length(sequences), "sequences\n")
  return(sequences)
}

#' Parse CSV metadata file
parse_csv_metadata <- function(csv_file) {
  cat("Reading CSV file:", csv_file, "\n")
  
  # Read with latin1 encoding for Norwegian characters
  df <- read.csv(csv_file, sep = ";", stringsAsFactors = FALSE, 
                 fileEncoding = "latin1", na.strings = c("", "NA"))
  
  # Standardize column names
  names(df) <- trimws(names(df))
  
  cat("  Found", nrow(df), "rows\n")
  cat("  Columns:", paste(names(df), collapse = ", "), "\n")
  
  return(df)
}

#' Main function to prepare clean metadata
prepare_metadata <- function(fasta_file = "export/2PA.fa", 
                             csv_file = "export/export.csv",
                             output_dir = "data") {
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("METADATA PARSING AND CLEANING\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Parse FASTA
  sequences <- parse_fasta(fasta_file)
  
  # Parse CSV
  csv_data <- parse_csv_metadata(csv_file)
  csv_data$Key <- as.character(csv_data$Key)
  
  # Process each sequence
  cat("\nProcessing and cleaning metadata...\n")
  
  results <- list()
  
  for (sample_id in names(sequences)) {
    seq_info <- sequences[[sample_id]]
    
    # Find matching CSV row
    csv_row <- csv_data[csv_data$Key == sample_id, ]
    
    # Initialize record
    record <- list(
      sample_id = sample_id,
      sequence = seq_info$sequence,
      sequence_length = nchar(seq_info$sequence),
      sampling_date = NA_character_,
      sample_year = NA_integer_,
      genotype = seq_info$genotype,
      variant = seq_info$variant,
      patient_id = NA_character_,
      geo_location = NA_character_,
      geo_country = NA_character_,
      source = NA_character_,
      transmission_route = NA_character_,
      comment = NA_character_
    )
    
    # Merge CSV data if available
    if (nrow(csv_row) > 0) {
      csv_row <- csv_row[1, ]  # Take first match
      
      # Date
      if (!is.na(csv_row$Sample.date)) {
        record$sampling_date <- convert_date(csv_row$Sample.date)
      }
      
      # Year
      if (!is.na(csv_row$Year)) {
        record$sample_year <- as.integer(csv_row$Year)
      }
      
      # Genotype (prefer CSV over FASTA if available)
      if (!is.na(csv_row$Genotype) && csv_row$Genotype != "") {
        record$genotype <- csv_row$Genotype
      }
      
      # Variant (prefer CSV over FASTA if available)
      if (!is.na(csv_row$OUTBREAK_VARIANT) && csv_row$OUTBREAK_VARIANT != "") {
        record$variant <- csv_row$OUTBREAK_VARIANT
      }
      
      # Patient ID
      if (!is.na(csv_row$Patient.initials)) {
        record$patient_id <- csv_row$Patient.initials
      }
      
      # Origin -> Country and Location
      if (!is.na(csv_row$Origin)) {
        origin_info <- clean_origin(csv_row$Origin)
        record$geo_country <- origin_info$country
        record$geo_location <- origin_info$location
        if (!is.na(origin_info$comment)) {
          record$comment <- origin_info$comment
        }
      }
      
      # Source
      if (!is.na(csv_row$Source)) {
        record$source <- csv_row$Source
      }
      
      # Transmission route
      if (!is.na(csv_row$Transmission.route)) {
        record$transmission_route <- csv_row$Transmission.route
      }
    }
    
    results[[sample_id]] <- record
  }
  
  # Convert to data frames
  cat("Creating output data frames...\n")
  
  # Sequences data frame
  sequences_df <- data.frame(
    sample_id = sapply(results, function(x) x$sample_id),
    sequence = sapply(results, function(x) x$sequence),
    sequence_length = sapply(results, function(x) x$sequence_length),
    stringsAsFactors = FALSE
  )
  
  # Metadata data frame
  metadata_df <- data.frame(
    sample_id = sapply(results, function(x) x$sample_id),
    sampling_date = sapply(results, function(x) x$sampling_date),
    sample_year = sapply(results, function(x) x$sample_year),
    genotype = sapply(results, function(x) x$genotype),
    variant = sapply(results, function(x) x$variant),
    patient_id = sapply(results, function(x) x$patient_id),
    geo_location = sapply(results, function(x) x$geo_location),
    geo_country = sapply(results, function(x) x$geo_country),
    source = sapply(results, function(x) x$source),
    transmission_route = sapply(results, function(x) x$transmission_route),
    comment = sapply(results, function(x) x$comment),
    stringsAsFactors = FALSE
  )
  
  # Replace "NA" strings with actual NA
  metadata_df[metadata_df == "NA"] <- NA
  
  # Save to CSV
  seq_output <- file.path(output_dir, "cleaned_sequences.csv")
  meta_output <- file.path(output_dir, "cleaned_metadata.csv")
  
  write.csv(sequences_df, seq_output, row.names = FALSE, na = "")
  write.csv(metadata_df, meta_output, row.names = FALSE, na = "")
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("OUTPUT SUMMARY\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("Sequences:", nrow(sequences_df), "\n")
  cat("Saved to:", seq_output, "\n\n")
  cat("Metadata:", nrow(metadata_df), "\n")
  cat("Saved to:", meta_output, "\n\n")
  
  # Show statistics
  cat("Genotype distribution:\n")
  print(table(metadata_df$genotype, useNA = "ifany"))
  
  cat("\nCountry distribution:\n")
  print(table(metadata_df$geo_country, useNA = "ifany"))
  
  cat("\nYear distribution:\n")
  print(table(metadata_df$sample_year, useNA = "ifany"))
  
  cat("\nRecords with comments:", sum(!is.na(metadata_df$comment)), "\n")
  
  return(list(sequences = sequences_df, metadata = metadata_df))
}

# =============================================================================
# Run if called directly
# =============================================================================
if (!interactive()) {
  prepare_metadata()
}
