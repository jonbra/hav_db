# R-only BLAST-like nearest-match helper (uses Biostrings pairwiseAlignment)
# Provides: run_blast(query, subjects/subject_fasta/con, max_hits, min_identity)
library(Biostrings)

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

# run_blast: R-only nearest-match search using pairwiseAlignment
# query: DNAStringSet or FASTA path or character vector
# subjects/subject_fasta/con: subject definitions (DB) or NULL to error
# returns list(hits = data.frame, sequences = DNAStringSet)
run_blast <- function(query,
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
      identity = identities,
      score = scores,
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$identity) & df$identity >= as.numeric(min_identity), , drop = FALSE]
    if (nrow(df) > 0) {
      df <- df[order(-df$identity, -df$score), , drop = FALSE]
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
