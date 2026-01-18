manova_var_explained <- function(D, group) {
  # D: distance matrix (n x n) or dist object, containing UNSQUARED distances
  # group: length-n grouping vector/factor
  
  n <- nrow(D)
  stopifnot(n == ncol(D))
  
  # basic checks (optional but helpful)
  if (anyNA(D)) stop("D has NA values.")
  if (max(abs(D - t(D))) > 1e-8) stop("D must be symmetric.")
  if (max(abs(diag(D))) > 1e-8) stop("diag(D) should be ~0 for a distance matrix.")
  
  # convert one-hot to group labels if needed
  if (length(group) != n) stop("group must have length n.")
  group <- factor(group)
  
  # PERMANOVA uses squared distances
  D2 <- D^2
  
  # total SS
  SST <- sum(D2) / (2 * n)
  
  # within-group SS
  idx_list <- split(seq_len(n), group)
  SSW <- sum(vapply(idx_list, function(ii) {
    nk <- length(ii)
    sum(D2[ii, ii, drop = FALSE]) / (2 * nk)
  }, numeric(1)))
  
  # explained (between / model) SS
  SSA <- SST - SSW
  
  # percent explained
  R2 <- SSA / SST
  
  list(
    SST = SST,
    SSW = SSW,
    SSA = SSA,
    R2 = R2
  )
}

W2_gauss_1d <- function(mu1, var1, mu2, var2) {
  s1 <- sqrt(var1); s2 <- sqrt(var2)
  sqrt((mu1 - mu2)^2 + (s1 - s2)^2)
}

