# BLAST wrapper functions for HAV Database
# Provides both real NCBI BLAST+ and fallback R-only nearest-match search
library(Biostrings)

# =============================================================================
# NCBI BLAST+ Functions (requires blast conda package)
# =============================================================================

#' Create a BLAST database from sequences
#' @param fasta_path Path to input FASTA file
#' @param db_path Output path for BLAST database (without extension)
#' @param db_type Database type: "nucl" or "prot"
#' @return Path to created database
make_blast_db <- function(fasta_path, db_path, db_type = "nucl") {
  if (!file.exists(fasta_path)) stop("FASTA file not found: ", fasta_path)
  
  # Get makeblastdb binary from conda environment
  makeblastdb_bin <- tryCatch(
    get_conda_bin("makeblastdb"),
    error = function(e) NULL
  )
  if (is.null(makeblastdb_bin) || !nzchar(makeblastdb_bin)) {
    stop("makeblastdb not found. Install BLAST+: conda env update -f conda/environment.yml")
  }
  
  # Create output directory if needed
  dir.create(dirname(db_path), showWarnings = FALSE, recursive = TRUE)
  
  args <- c("-in", fasta_path, "-dbtype", db_type, "-out", db_path)
  stderr_file <- tempfile("makeblastdb_stderr_")
  status <- system2(makeblastdb_bin, args = args, stdout = TRUE, stderr = stderr_file)
  
  if (!file.exists(paste0(db_path, ".ndb")) && !file.exists(paste0(db_path, ".nin"))) {
    stderr_msg <- tryCatch(paste(readLines(stderr_file), collapse = "\n"), error = function(e) "")
    stop("makeblastdb failed to create database. stderr: ", stderr_msg)
  }
  
  db_path
}

#' Run BLASTn search
#' @param query_fasta Path to query FASTA file
#' @param db_path Path to BLAST database
#' @param max_hits Maximum number of hits per query
#' @param evalue E-value threshold
#' @param num_threads Number of threads
#' @return Data frame with BLAST results (outfmt 6 columns + SNP count)
run_blastn <- function(query_fasta, db_path, max_hits = 10, evalue = 10, num_threads = 2) {
  if (!file.exists(query_fasta)) stop("Query FASTA not found: ", query_fasta)
  
  # Get blastn binary from conda environment
  blastn_bin <- tryCatch(
    get_conda_bin("blastn"),
    error = function(e) NULL
  )
  if (is.null(blastn_bin) || !nzchar(blastn_bin)) {
    stop("blastn not found. Install BLAST+: conda env update -f conda/environment.yml")
  }
  
  output_file <- tempfile("blastn_out_", fileext = ".tsv")
  
  # BLAST outfmt 6 columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
  args <- c(
    "-query", query_fasta,
    "-db", db_path,
    "-outfmt", "6",
    "-max_target_seqs", as.character(max_hits),
    "-evalue", as.character(evalue),
    "-num_threads", as.character(num_threads),
    "-out", output_file
  )
  
  stderr_file <- tempfile("blastn_stderr_")
  status <- system2(blastn_bin, args = args, stdout = TRUE, stderr = stderr_file)
  
  if (status != 0) {
    stderr_msg <- tryCatch(paste(readLines(stderr_file), collapse = "\n"), error = function(e) "")
    stop("blastn failed (status ", status, "): ", stderr_msg)
  }
  
  # Parse tabular output
  if (!file.exists(output_file) || file.info(output_file)$size == 0) {
    return(data.frame(
      query_id = character(),
      subject_id = character(),
      identity_pct = numeric(),
      alignment_length = integer(),
      mismatches = integer(),
      gap_opens = integer(),
      query_start = integer(),
      query_end = integer(),
      subject_start = integer(),
      subject_end = integer(),
      evalue = numeric(),
      bit_score = numeric(),
      snp_count = integer(),
      stringsAsFactors = FALSE
    ))
  }
  
  hits <- read.table(
    output_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("query_id", "subject_id", "identity_pct", "alignment_length",
                  "mismatches", "gap_opens", "query_start", "query_end",
                  "subject_start", "subject_end", "evalue", "bit_score")
  )
  
  # Calculate SNP count from alignment length and identity
  # SNPs ≈ alignment_length × (1 - identity/100)
  hits$snp_count <- round(hits$alignment_length * (1 - hits$identity_pct / 100))
  
  hits
}

#' Run BLAST search against database sequences
#' @param query DNAStringSet, FASTA path, or character vector of sequences
#' @param con Database connection
#' @param max_hits Maximum hits per query
#' @param evalue E-value threshold
#' @param num_threads Number of threads
#' @return List with hits data frame and subject sequences
run_blast_search <- function(query, con, max_hits = 10, evalue = 10, num_threads = 2) {
  # Convert query to DNAStringSet
  qset <- read_query_to_DNAStringSet(query)
  if (length(qset) == 0) {
    return(list(hits = data.frame(), sequences = DNAStringSet()))
  }
  
  # Export all database sequences to temp FASTA
  db_seqs <- get_sequences_with_metadata(con)
  if (nrow(db_seqs) == 0) {
    return(list(hits = data.frame(), sequences = DNAStringSet()))
  }
  
  # Create temp files
  query_fasta <- tempfile("blast_query_", fileext = ".fa")
  subject_fasta <- tempfile("blast_subjects_", fileext = ".fa")
  db_path <- tempfile("blast_db_")
  
  # Write query sequences
  Biostrings::writeXStringSet(qset, filepath = query_fasta, format = "fasta")
  
  # Write subject sequences
  subj_set <- DNAStringSet(toupper(db_seqs$sequence))
  names(subj_set) <- db_seqs$sample_id
  Biostrings::writeXStringSet(subj_set, filepath = subject_fasta, format = "fasta")
  
  # Create BLAST database
  make_blast_db(subject_fasta, db_path)
  
  # Run BLAST
  hits <- run_blastn(query_fasta, db_path, max_hits = max_hits, evalue = evalue, num_threads = num_threads)
  
  # Clean up temp files
  unlink(c(query_fasta, subject_fasta))
  unlink(list.files(dirname(db_path), pattern = basename(db_path), full.names = TRUE))
  
  if (nrow(hits) == 0) {
    return(list(hits = data.frame(), sequences = DNAStringSet()))
  }
  
  # Get hit sequences
  hit_ids <- unique(hits$subject_id)
  hit_seqs <- subj_set[hit_ids]
  
  list(hits = hits, sequences = hit_seqs)
}

# =============================================================================
# R-only Fallback Functions (no external dependencies)
# =============================================================================

read_query_to_DNAStringSet <- function(x) {
  if (inherits(x, "DNAStringSet")) return(x)
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    return(readDNAStringSet(x, format = "fasta"))
  }
  if (is.character(x)) {
    seqs <- DNAStringSet(toupper(x))
    names(seqs) <- paste0("query_", seq_along(seqs))
    return(seqs)
  }
  stop("Unsupported query type; must be DNAStringSet, FASTA path, or character vector")
}

read_subjects_to_DNAStringSet <- function(subjects = NULL, subject_fasta = NULL, con = NULL) {
  if (!is.null(subjects)) {
    if (inherits(subjects, "DNAStringSet")) return(subjects)
    if (is.character(subjects) && length(subjects) == 1 && file.exists(subjects)) {
      return(readDNAStringSet(subjects, format = "fasta"))
    }
    if (is.character(subjects)) {
      s <- DNAStringSet(toupper(subjects))
      if (is.null(names(s))) names(s) <- paste0("subj_", seq_along(s))
      return(s)
    }
  }
  if (!is.null(subject_fasta) && file.exists(subject_fasta)) {
    return(readDNAStringSet(subject_fasta, format = "fasta"))
  }
  if (!is.null(con)) {
    if (!exists("get_sequences_with_metadata")) stop("DB helper `get_sequences_with_metadata` not found")
    df <- get_sequences_with_metadata(con)
    if (nrow(df) == 0) return(DNAStringSet())
    seqs <- DNAStringSet(toupper(df$sequence))
    names(seqs) <- df$sample_id
    return(seqs)
  }
  stop("No valid subjects provided")
}

compute_identity <- function(aln) {
  ap <- as.character(alignedPattern(aln))
  asub <- as.character(alignedSubject(aln))
  a_chars <- strsplit(ap, "")[[1]]
  s_chars <- strsplit(asub, "")[[1]]
  valid_pos <- (a_chars != "-" & s_chars != "-")
  if (sum(valid_pos) == 0) return(0)
  matches <- sum(a_chars[valid_pos] == s_chars[valid_pos])
  matches / sum(valid_pos) * 100
}

#' R-only nearest-match search using pairwiseAlignment (fallback when BLAST unavailable)
#' @param query DNAStringSet, FASTA path, or character vector
#' @param subjects DNAStringSet or NULL
#' @param subject_fasta Path to subject FASTA or NULL
#' @param con Database connection or NULL
#' @param max_hits Maximum hits per query
#' @param min_identity Minimum identity threshold (0-100)
#' @param alignment_type Alignment type for pairwiseAlignment
#' @return List with hits data frame and subject sequences
run_blast_fallback <- function(query,
                      subjects = NULL,
                      subject_fasta = NULL,
                      con = NULL,
                      max_hits = 10,
                      min_identity = 0,
                      alignment_type = "local") {
  qset <- read_query_to_DNAStringSet(query)
  sset <- read_subjects_to_DNAStringSet(subjects = subjects, subject_fasta = subject_fasta, con = con)

  if (length(qset) == 0 || length(sset) == 0) {
    return(list(hits = data.frame(), sequences = DNAStringSet()))
  }

  hits_list <- list()
  for (qi in seq_along(qset)) {
    qname <- names(qset)[qi]
    qseq <- qset[[qi]]
    scores <- numeric(length(sset))
    identities <- numeric(length(sset))
    for (si in seq_along(sset)) {
      sseq <- sset[[si]]
      aln <- tryCatch(pairwiseAlignment(pattern = qseq, subject = sseq, type = alignment_type, scoreOnly = FALSE), error = function(e) NULL)
      if (is.null(aln)) {
        scores[si] <- NA
        identities[si] <- NA
        next
      }
      identities[si] <- compute_identity(aln)
      scores[si] <- score(aln)
    }
    df <- data.frame(
      query_id = rep(qname, length(sset)),
      subject_id = names(sset),
      identity_pct = identities,
      score = scores,
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$identity_pct) & df$identity_pct >= as.numeric(min_identity), , drop = FALSE]
    if (nrow(df) > 0) {
      df <- df[order(-df$identity_pct, -df$score), , drop = FALSE]
      df <- head(df, as.integer(max_hits))
    }
    hits_list[[qname]] <- df
  }

  all_hits <- do.call(rbind, lapply(names(hits_list), function(nm) {
    if (nrow(hits_list[[nm]]) == 0) return(NULL)
    hits_list[[nm]]
  }))
  if (is.null(all_hits) || nrow(all_hits) == 0) {
    return(list(hits = data.frame(), sequences = DNAStringSet()))
  }
  subj_ids <- unique(all_hits$subject_id)
  subj_seqs <- sset[subj_ids]
  return(list(hits = all_hits, sequences = subj_seqs))
}

#' Check if BLAST+ is available
#' @return TRUE if blastn is available, FALSE otherwise
blast_available <- function() {
  tryCatch({
    bin <- get_conda_bin("blastn")
    nzchar(bin) && file.exists(bin)
  }, error = function(e) FALSE)
}
