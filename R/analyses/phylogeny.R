# Phylogeny builders (NJ, UPGMA, IQ-TREE wrapper)

build_nj_tree_from_alignment <- function(aligned_dna, model = "K80") {
  dnabin <- ape::as.DNAbin(as.matrix(aligned_dna))
  distm <- ape::dist.dna(dnabin, model = model, pairwise.deletion = TRUE)
  tr <- ape::nj(distm)
  ape::ladderize(tr)
}

build_upgma_tree_from_alignment <- function(aligned_dna, model = "K80") {
  dnabin <- ape::as.DNAbin(as.matrix(aligned_dna))
  distm <- ape::dist.dna(dnabin, model = model, pairwise.deletion = TRUE)
  hc <- hclust(as.dist(distm), method = "average")
  tr <- ape::as.phylo(hc)
  ape::ladderize(tr)
}

run_iqtree <- function(alignment_fasta, prefix = tempfile("iq"), threads = 1, extra_args = NULL) {
  # Try common binary names
  bin <- Sys.which("iqtree")
  if (!nzchar(bin)) bin <- Sys.which("iqtree2")
  if (!nzchar(bin)) stop("iqtree not found on PATH")

  # Use --redo by default to overwrite previous runs/checkpoints
  args <- c("--redo", "-s", alignment_fasta, "-nt", as.character(threads), "-pre", prefix)
  if (!is.null(extra_args)) args <- c(args, extra_args)

  outf <- tempfile("iq_out")
  errf <- tempfile("iq_err")
  status <- system2(bin, args = args, stdout = outf, stderr = errf)
  if (status != 0) {
    stderr_msg <- tryCatch(paste(readLines(errf), collapse = "\n"), error = function(e) "(failed to read iqtree stderr)")
    stop(sprintf("iqtree failed (status %d): %s", status, stderr_msg))
  }

  treefile <- paste0(prefix, ".treefile")
  if (!file.exists(treefile)) {
    stop("iqtree did not produce treefile; check iqtree output and permissions")
  }

  ape::read.tree(treefile)
}
