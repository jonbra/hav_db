# Alignment utilities using MAFFT
# Pure R functions, no Shiny/reactivity. Deterministic I/O.

# Helper function to get binary path from hav_db conda environment
get_conda_bin <- function(tool) {
  # Try hav_db conda environment first
  conda_base <- Sys.getenv("CONDA_PREFIX_1", Sys.getenv("CONDA_PREFIX", ""))
  if (nzchar(conda_base)) {
    # Check in hav_db environment under the conda base
    hav_db_bin <- file.path(dirname(conda_base), "envs", "hav_db", "bin", tool)
    if (file.exists(hav_db_bin)) return(hav_db_bin)
  }
  
  # Common conda installation paths
  home <- Sys.getenv("HOME")
  conda_paths <- c(
    file.path(home, "miniforge3", "envs", "hav_db", "bin", tool),
    file.path(home, "miniconda3", "envs", "hav_db", "bin", tool),
    file.path(home, "anaconda3", "envs", "hav_db", "bin", tool),
    file.path("/opt", "conda", "envs", "hav_db", "bin", tool)
  )
  
  for (p in conda_paths) {
    if (file.exists(p)) return(p)
  }
  
  # Fall back to PATH
  path_bin <- Sys.which(tool)
  if (nzchar(path_bin)) return(path_bin)
  
  stop(sprintf("%s not found. Ensure the hav_db conda environment is set up: conda env create -f conda/environment.yml", tool))
}

write_fasta_from_strings <- function(seqs, ids = NULL, path) {
  if (is.null(ids)) ids <- paste0("seq", seq_along(seqs))
  if (length(ids) != length(seqs)) stop("ids must match seqs length")
  ss <- Biostrings::DNAStringSet(seqs)
  names(ss) <- ids
  Biostrings::writeXStringSet(ss, filepath = path, format = "fasta")
  path
}

run_mafft <- function(input_fasta, output_fasta, threads = 1, opts = "--auto") {
  bin <- get_conda_bin("mafft")
  args <- c("--thread", as.character(threads), opts, input_fasta)
  res <- system2(bin, args = args, stdout = output_fasta, stderr = TRUE)
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
