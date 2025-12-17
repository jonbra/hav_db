# Distance and SNP utilities

calculate_distance_matrix_from_alignment <- function(aligned_dna, model = "K80", as_percent = FALSE) {
  # aligned_dna: DNAStringSet
  if (!inherits(aligned_dna, "DNAStringSet")) stop("aligned_dna must be DNAStringSet")
  mat <- as.matrix(aligned_dna)
  dnabin <- ape::as.DNAbin(mat)
  distm <- ape::dist.dna(dnabin, model = model, pairwise.deletion = TRUE)
  if (as_percent) {
    as.matrix((1 - as.matrix(distm)) * 100)
  } else {
    distm
  }
}

compute_snp_distance <- function(aligned_dna) {
  # simple SNP count (pairwise) ignoring gaps
  if (!inherits(aligned_dna, "DNAStringSet")) stop("aligned_dna must be DNAStringSet")
  mat <- as.matrix(aligned_dna)
  n <- nrow(mat)
  names <- rownames(mat)
  out <- matrix(0, n, n, dimnames = list(names, names))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    a <- mat[i, ]
    b <- mat[j, ]
    valid <- !(a == "-" | b == "-")
    out[i, j] <- sum(a[valid] != b[valid])
  }
  out
}
