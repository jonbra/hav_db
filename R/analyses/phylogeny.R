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

run_iqtree <- function(alignment_fasta, prefix = tempfile("iq"), threads = 2, 
                       model = "GTR+G+I", bootstrap = 1000, extra_args = NULL) {
  # Use iqtree from hav_db conda environment (IQ-TREE 3.x)
  # The get_conda_bin helper is defined in alignment.R and sourced before this file
  bin <- tryCatch(
    get_conda_bin("iqtree"),
    error = function(e) NULL
  )
  if (is.null(bin) || !nzchar(bin)) {
    stop("iqtree not found. Ensure the hav_db conda environment is set up: conda env create -f conda/environment.yml")
  }

  # Use --redo by default to overwrite previous runs/checkpoints
  # Use -T for threads (IQ-TREE 2.x/3.x syntax, replaces -nt)
  args <- c("--redo", "-s", alignment_fasta, "-T", as.character(threads), "-pre", prefix)
  
  # Add substitution model (default GTR+G+I)
  if (!is.null(model) && nzchar(model)) {
    args <- c(args, "-m", model)
  }
  
  # Add ultrafast bootstrap replicates if bootstrap > 0
  if (!is.null(bootstrap) && !is.na(bootstrap) && bootstrap > 0) {
    args <- c(args, "-B", as.character(as.integer(bootstrap)))
  }
  
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
