# =============================================================================
# Analysis Functions for HAV Database
# =============================================================================
# Multiple sequence alignment, phylogenetic tree construction, and clustering.

library(Biostrings)
library(ape)
library(jsonlite)

# Source alignment helpers first (provides get_conda_bin function)
aln_file <- NULL
if (exists("app_dir") && nzchar(app_dir)) {
  aln_file <- file.path(app_dir, "R", "analyses", "alignment.R")
} else {
  aln_file <- file.path("R", "analyses", "alignment.R")
}
if (!is.null(aln_file) && file.exists(aln_file)) {
  try(source(aln_file), silent = TRUE)
}

# Source phylogeny wrappers (depends on get_conda_bin from alignment.R)
phy_file <- NULL
if (exists("app_dir") && nzchar(app_dir)) {
  phy_file <- file.path(app_dir, "R", "analyses", "phylogeny.R")
} else {
  phy_file <- file.path("R", "analyses", "phylogeny.R")
}
if (!is.null(phy_file) && file.exists(phy_file)) {
  try(source(phy_file), silent = TRUE)
}

# Source BLAST wrapper (provides run_blast_search, blast_available, etc.)
blast_file <- NULL
if (exists("app_dir") && nzchar(app_dir)) {
  blast_file <- file.path(app_dir, "R", "analyses", "blast_wrapper.R")
} else {
  blast_file <- file.path("R", "analyses", "blast_wrapper.R")
}
if (!is.null(blast_file) && file.exists(blast_file)) {
  try(source(blast_file), silent = TRUE)
}

# NOTE: the `msa` Bioconductor package is not used. Alignments are performed
# by calling the external `mafft` binary (installed with conda or package
# manager). This avoids compiling bundled C/C++ code and keeps the runtime
# environment lightweight.

# =============================================================================
# Multiple Sequence Alignment
# =============================================================================

#' Perform multiple sequence alignment
#' @param dna_set DNAStringSet object
#' @param method Alignment method: currently uses external "mafft" via system call
#' @return MsaDNAMultipleAlignment object
run_msa <- function(dna_set, method = "mafft") {
  # dna_set: DNAStringSet or character vector of sequences
  if (length(dna_set) < 2) stop("Need at least 2 sequences for alignment")

  # We only support mafft via system call for now. The `method` argument is
  # kept for API compatibility but ignored when using mafft.
  # Use get_conda_bin to find mafft in hav_db conda environment
  mafft_bin <- tryCatch(
    get_conda_bin("mafft"),
    error = function(e) NULL
  )
  if (is.null(mafft_bin) || !nzchar(mafft_bin)) {
    stop("mafft not found. Ensure the hav_db conda environment is set up: conda env create -f conda/environment.yml")
  }

  # Write input sequences to a temporary FASTA file
  in_fa <- tempfile(fileext = ".fa")
  out_fa <- tempfile(fileext = ".aligned.fa")

  # Accept DNAStringSet or character vector
  if (inherits(dna_set, "DNAStringSet")) {
    Biostrings::writeXStringSet(dna_set, filepath = in_fa, format = "fasta")
  } else if (is.character(dna_set)) {
    seqs <- Biostrings::DNAStringSet(dna_set)
    Biostrings::writeXStringSet(seqs, filepath = in_fa, format = "fasta")
  } else {
    stop("dna_set must be a DNAStringSet or character vector")
  }

  # Run mafft --auto for reasonable defaults. Write stdout to out_fa.
  args <- c("--auto", in_fa)
  res <- system2(mafft_bin, args = args, stdout = out_fa, stderr = tempfile())
  if (res != 0) stop("mafft alignment failed (check installation and input sequences)")

  # Read aligned sequences back into R as a DNAStringSet
  aligned <- Biostrings::readDNAStringSet(out_fa, format = "fasta")

  aligned
}

#' Run MSA on database sequences
#' @param con Database connection
#' @param sample_ids Vector of sample IDs
#' @param method Alignment method
#' @return MsaDNAMultipleAlignment object
run_msa_from_db <- function(con, sample_ids, method = "mafft") {
  dna_set <- db_to_dna_stringset(con, sample_ids)
  run_msa(dna_set, method)
}

#' Convert MSA to DNAbin format (for ape)
#' @param msa_result MsaDNAMultipleAlignment object
#' @return DNAbin object
msa_to_dnabin <- function(msa_result) {
  # Accept either a DNAStringSet (returned by our mafft wrapper) or an
  # MsaDNAMultipleAlignment object from the `msa` package. Convert to
  # a character matrix and then to DNAbin.
  if (inherits(msa_result, "DNAStringSet")) {
    aligned_seqs <- msa_result
  } else if (inherits(msa_result, "MsaDNAMultipleAlignment") || inherits(msa_result, "MsaAAMultipleAlignment")) {
    aligned_seqs <- as(msa_result, "DNAStringSet")
  } else {
    stop("Unsupported msa_result type for msa_to_dnabin")
  }

  seq_matrix <- as.matrix(aligned_seqs)
  ape::as.DNAbin(seq_matrix)
}

#' Export alignment to FASTA
#' @param msa_result MsaDNAMultipleAlignment object
#' @param output_file Path for output file
#' @return Path to output file
export_alignment_fasta <- function(msa_result, output_file) {
  if (inherits(msa_result, "DNAStringSet")) {
    aligned_seqs <- msa_result
  } else if (inherits(msa_result, "MsaDNAMultipleAlignment") || inherits(msa_result, "MsaAAMultipleAlignment")) {
    aligned_seqs <- as(msa_result, "DNAStringSet")
  } else {
    stop("Unsupported msa_result type for export_alignment_fasta")
  }
  # Ensure sequence names are present (IQ-TREE requires non-empty unique names)
  seq_names <- names(aligned_seqs)
  if (is.null(seq_names) || any(nzchar(seq_names) == FALSE)) {
    names(aligned_seqs) <- paste0("seq", seq_len(length(aligned_seqs)))
  }
  Biostrings::writeXStringSet(aligned_seqs, output_file, format = "fasta")
  output_file
}

#' Get alignment summary statistics
#' @param msa_result MsaDNAMultipleAlignment object
#' @return List with alignment statistics
alignment_stats <- function(msa_result) {
  if (inherits(msa_result, "DNAStringSet")) {
    aligned_seqs <- msa_result
  } else if (inherits(msa_result, "MsaDNAMultipleAlignment") || inherits(msa_result, "MsaAAMultipleAlignment")) {
    aligned_seqs <- as(msa_result, "DNAStringSet")
  } else {
    stop("Unsupported msa_result type for alignment_stats")
  }

  aln_matrix <- as.matrix(aligned_seqs)
  
  # Calculate statistics
  n_seqs <- nrow(aln_matrix)
  aln_length <- ncol(aln_matrix)
  
  # Gap statistics
  gap_counts <- apply(aln_matrix, 2, function(x) sum(x == "-"))
  total_gaps <- sum(gap_counts)
  gap_positions <- sum(gap_counts > 0)
  
  # Conservation (positions with no gaps and all same nucleotide)
  conserved <- sum(apply(aln_matrix, 2, function(x) {
    x_no_gap <- x[x != "-"]
    length(unique(x_no_gap)) == 1 && length(x_no_gap) == length(x)
  }))
  
  list(
    n_sequences = n_seqs,
    alignment_length = aln_length,
    total_gaps = total_gaps,
    gap_positions = gap_positions,
    conserved_positions = conserved,
    conservation_pct = round(conserved / aln_length * 100, 1),
    sequence_names = rownames(aln_matrix)
  )
}

# =============================================================================
# Phylogenetic Tree Construction
# =============================================================================

#' Build neighbor-joining tree from alignment
#' @param msa_result MsaDNAMultipleAlignment object
#' @param model Distance model: "raw", "K80", "K81", "F81", "F84", "T92", "TN93", "JC69"
#' @return phylo object
build_nj_tree <- function(msa_result, model = "K80") {
  dna_bin <- msa_to_dnabin(msa_result)
  
  # Calculate distance matrix
  dist_matrix <- ape::dist.dna(dna_bin, model = model, pairwise.deletion = TRUE)
  
  # Handle any NA values
  if (any(is.na(dist_matrix))) {
    warning("Some distances could not be calculated, using pairwise deletion")
    dist_matrix[is.na(dist_matrix)] <- max(dist_matrix, na.rm = TRUE)
  }
  
  # Build NJ tree
  tree <- ape::nj(dist_matrix)
  
  # Root at midpoint
  tree <- ape::ladderize(tree)
  
  tree
}

#' Build UPGMA tree from alignment
#' @param msa_result MsaDNAMultipleAlignment object
#' @param model Distance model
#' @return phylo object
build_upgma_tree <- function(msa_result, model = "K80") {
  dna_bin <- msa_to_dnabin(msa_result)
  dist_matrix <- ape::dist.dna(dna_bin, model = model, pairwise.deletion = TRUE)
  
  if (any(is.na(dist_matrix))) {
    dist_matrix[is.na(dist_matrix)] <- max(dist_matrix, na.rm = TRUE)
  }
  
  # UPGMA clustering
  hc <- hclust(dist_matrix, method = "average")
  tree <- ape::as.phylo(hc)
  
  ape::ladderize(tree)
}

#' Build tree from database sequences
#' @param con Database connection
#' @param sample_ids Vector of sample IDs
#' @param method Tree method: "nj" or "upgma"
#' @param alignment_method MSA method
#' @param dist_model Distance model
#' @return phylo object
build_tree_from_db <- function(con, sample_ids, method = "nj", 
                                alignment_method = "mafft", dist_model = "K80") {
  # Run alignment
  msa_result <- run_msa_from_db(con, sample_ids, alignment_method)
  
  # Build tree
  if (method == "nj") {
    build_nj_tree(msa_result, dist_model)
  } else if (method == "upgma") {
    build_upgma_tree(msa_result, dist_model)
  } else if (method == "iqtree" || method == "iq") {
    # Use iqtree external binary. Export alignment to fasta and call wrapper.
    fa <- tempfile(fileext = ".fa")
    export_alignment_fasta(msa_result, fa)
    # Prefix for iqtree files
    pref <- tempfile("iqtree")
    if (!exists("run_iqtree")) stop("IQ-TREE support not available (run_iqtree missing)")
    run_iqtree(fa, prefix = pref, threads = 1)
    treefile <- paste0(pref, ".treefile")
    if (!file.exists(treefile)) stop("IQ-TREE did not produce a treefile")
    ape::read.tree(treefile)
  } else {
    stop("Invalid method. Choose 'nj', 'upgma', or 'iqtree'")
  }
}

#' Export tree to Newick format
#' @param tree phylo object
#' @param output_file Path for output file
#' @return Path to output file
export_tree_newick <- function(tree, output_file) {
  ape::write.tree(tree, file = output_file)
  output_file
}

#' Get tree statistics
#' @param tree phylo object
#' @return List with tree statistics
tree_stats <- function(tree) {
  list(
    n_tips = length(tree$tip.label),
    n_internal_nodes = tree$Nnode,
    tip_labels = tree$tip.label,
    is_rooted = ape::is.rooted(tree),
    is_binary = ape::is.binary(tree),
    total_branch_length = sum(tree$edge.length)
  )
}

# =============================================================================
# Distance Matrix Functions
# =============================================================================

#' Calculate pairwise distance matrix
#' @param msa_result MsaDNAMultipleAlignment object
#' @param model Distance model
#' @param as_percent Return as percent identity instead of distance
#' @return Distance matrix
calculate_distance_matrix <- function(msa_result, model = "K80", as_percent = FALSE) {
  dna_bin <- msa_to_dnabin(msa_result)
  dist_matrix <- ape::dist.dna(dna_bin, model = model, pairwise.deletion = TRUE)
  
  if (as_percent) {
    # Convert to percent identity
    dist_matrix <- (1 - as.matrix(dist_matrix)) * 100
  }
  
  dist_matrix
}

#' Calculate distance matrix from database sequences
#' @param con Database connection
#' @param sample_ids Vector of sample IDs
#' @param model Distance model
#' @return Distance matrix
distance_matrix_from_db <- function(con, sample_ids, model = "K80") {
  msa_result <- run_msa_from_db(con, sample_ids, "mafft")
  calculate_distance_matrix(msa_result, model)
}

# =============================================================================
# Clustering Functions
# =============================================================================

#' Cluster sequences by distance threshold
#' @param msa_result MsaDNAMultipleAlignment object
#' @param threshold Distance threshold for clustering
#' @param model Distance model
#' @return Data frame with sample_id and cluster assignment
cluster_by_distance <- function(msa_result, threshold = 0.03, model = "K80") {
  dist_matrix <- calculate_distance_matrix(msa_result, model)
  
  # Hierarchical clustering
  hc <- hclust(dist_matrix, method = "complete")
  
  # Cut tree at threshold
  clusters <- cutree(hc, h = threshold)
  
  data.frame(
    sample_id = names(clusters),
    cluster = clusters,
    stringsAsFactors = FALSE
  )
}

#' Cluster sequences by number of clusters
#' @param msa_result MsaDNAMultipleAlignment object
#' @param k Number of clusters
#' @param model Distance model
#' @return Data frame with sample_id and cluster assignment
cluster_by_k <- function(msa_result, k, model = "K80") {
  dist_matrix <- calculate_distance_matrix(msa_result, model)
  
  hc <- hclust(dist_matrix, method = "complete")
  clusters <- cutree(hc, k = k)
  
  data.frame(
    sample_id = names(clusters),
    cluster = clusters,
    stringsAsFactors = FALSE
  )
}

#' Cluster database sequences
#' @param con Database connection
#' @param sample_ids Vector of sample IDs
#' @param threshold Distance threshold (or NULL to use k)
#' @param k Number of clusters (or NULL to use threshold)
#' @return Data frame with cluster assignments
cluster_from_db <- function(con, sample_ids, threshold = 0.03, k = NULL) {
  msa_result <- run_msa_from_db(con, sample_ids, "mafft")
  
  if (!is.null(k)) {
    cluster_by_k(msa_result, k)
  } else {
    cluster_by_distance(msa_result, threshold)
  }
}

# =============================================================================
# Result Serialization (for database storage)
# =============================================================================

#' Serialize analysis result to JSON
#' @param result Analysis result object
#' @param type Type of result: "alignment_stats", "tree", "clusters", "distance"
#' @return JSON string
serialize_result <- function(result, type) {
  if (type == "tree") {
    # Convert phylo to Newick string
    newick <- ape::write.tree(result)
    jsonlite::toJSON(list(
      type = "tree",
      newick = newick,
      stats = tree_stats(result)
    ), auto_unbox = TRUE)
  } else if (type == "alignment_stats") {
    jsonlite::toJSON(list(
      type = "alignment_stats",
      data = result
    ), auto_unbox = TRUE)
  } else if (type == "clusters") {
    jsonlite::toJSON(list(
      type = "clusters",
      data = result
    ), auto_unbox = TRUE)
  } else if (type == "distance") {
    jsonlite::toJSON(list(
      type = "distance",
      data = as.matrix(result)
    ), auto_unbox = TRUE)
  } else {
    jsonlite::toJSON(result, auto_unbox = TRUE)
  }
}

#' Deserialize analysis result from JSON
#' @param json_str JSON string
#' @return Original result object
deserialize_result <- function(json_str) {
  parsed <- jsonlite::fromJSON(json_str)
  
  if (parsed$type == "tree") {
    # Convert Newick back to phylo
    tree <- ape::read.tree(text = parsed$newick)
    return(tree)
  }
  
  parsed$data
}

# =============================================================================
# Analysis Pipeline (combines multiple steps)
# =============================================================================

#' Run complete phylogenetic analysis
#' @param con Database connection
#' @param sample_ids Vector of sample IDs
#' @param alignment_method MSA method
#' @param tree_method Tree building method
#' @param dist_model Distance model
#' @param save_to_db Save results to database
#' @param analysis_name Name for saved analysis
#' @return List with alignment, tree, and statistics
run_phylogenetic_analysis <- function(con, sample_ids, 
                                       alignment_method = "mafft",
                                       tree_method = "nj",
                                       dist_model = "K80",
                                       save_to_db = FALSE,
                                       analysis_name = NULL) {
  
  cat("Running multiple sequence alignment...\n")
  msa_result <- run_msa_from_db(con, sample_ids, alignment_method)
  aln_stats <- alignment_stats(msa_result)
  
  cat("Building phylogenetic tree...\n")
  if (tree_method == "nj") {
    tree <- build_nj_tree(msa_result, dist_model)
  } else if (tree_method == "upgma") {
    tree <- build_upgma_tree(msa_result, dist_model)
  } else if (tree_method == "iqtree" || tree_method == "iq") {
    # Export alignment and run IQ-TREE
    fa <- tempfile(fileext = ".fa")
    export_alignment_fasta(msa_result, fa)
    pref <- tempfile("iqtree")
    if (!exists("run_iqtree")) stop("IQ-TREE wrapper not available. Ensure R/analyses/phylogeny.R is present and sourced.")
    tree <- run_iqtree(fa, prefix = pref, threads = 1)
  } else {
    stop("Unsupported tree_method; choose 'nj', 'upgma', or 'iqtree'")
  }
  tree_info <- tree_stats(tree)
  
  cat("Calculating distance matrix...\n")
  dist_matrix <- calculate_distance_matrix(msa_result, dist_model)
  
  results <- list(
    alignment = msa_result,
    alignment_stats = aln_stats,
    tree = tree,
    tree_stats = tree_info,
    distance_matrix = dist_matrix,
    sample_ids = sample_ids,
    parameters = list(
      alignment_method = alignment_method,
      tree_method = tree_method,
      dist_model = dist_model
    )
  )
  
  # Optionally save to database
  if (save_to_db && !is.null(analysis_name)) {
    cat("Saving results to database...\n")
    result_json <- serialize_result(tree, "tree")
    params_json <- jsonlite::toJSON(results$parameters, auto_unbox = TRUE)
    
    save_analysis_result(con, 
      analysis_name = analysis_name,
      analysis_type = "phylogeny",
      sample_ids = sample_ids,
      result_data = result_json,
      parameters = params_json
    )
  }
  
  cat("Analysis complete!\n")
  results
}
