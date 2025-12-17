# Alignment utilities using MAFFT
# Pure R functions, no Shiny/reactivity. Deterministic I/O.

write_fasta_from_strings <- function(seqs, ids = NULL, path) {
  if (is.null(ids)) ids <- paste0("seq", seq_along(seqs))
  if (length(ids) != length(seqs)) stop("ids must match seqs length")
  ss <- Biostrings::DNAStringSet(seqs)
  names(ss) <- ids
  Biostrings::writeXStringSet(ss, filepath = path, format = "fasta")
  path
}

run_mafft <- function(input_fasta, output_fasta, threads = 1, opts = "--auto") {
  if (!nzchar(Sys.which("mafft"))) stop("mafft not found on PATH")
  args <- c("--thread", as.character(threads), opts, input_fasta)
  res <- system2("mafft", args = args, stdout = output_fasta, stderr = TRUE)
  if (!file.exists(output_fasta)) stop("mafft failed to produce output file")
  output_fasta
}

align_sequences <- function(seqs, ids = NULL, threads = 1) {
  # seqs: character vector or DNAStringSet
  in_fa <- tempfile(fileext = ".fa")
  out_fa <- tempfile(fileext = ".aligned.fa")
  if (inherits(seqs, "DNAStringSet")) {
    Biostrings::writeXStringSet(seqs, filepath = in_fa, format = "fasta")
  } else if (is.character(seqs)) {
    write_fasta_from_strings(seqs, ids = ids, path = in_fa)
  } else {
    stop("seqs must be character vector or DNAStringSet")
  }
  run_mafft(in_fa, out_fa, threads = threads)
  Biostrings::readDNAStringSet(out_fa, format = "fasta")
}
