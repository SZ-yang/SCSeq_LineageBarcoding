.remove_zero_variance <- function(data, 
                                  response_col, 
                                  group_col,
                                  tol = 1e-6){
  stopifnot(is.data.frame(data))
  
  # Convert grouping column to factor
  data[[group_col]] <- as.factor(data[[group_col]])
  groups <- levels(data[[group_col]])
  k <- length(groups)
  
  var_i  <- numeric(k)
  
  for (i in seq_along(groups)) {
    subset_data <- data[data[[group_col]] == groups[i], response_col, drop = TRUE]
    var_i[i]    <- var(subset_data)
  }
  
  names(var_i) <- groups
  
  idx <- which(var_i < tol)
  if(length(idx) > 0){
    rm_idx <- which(data[[group_col]] %in% names(var_i)[idx])
    data <- data[-rm_idx,]
    
    data[[group_col]] <- droplevels(data[[group_col]])
  }
  
  return(data)
}

.welch_anova <- function(data, response_col, group_col) {
  # Convert grouping column to factor
  data[[group_col]] <- as.factor(data[[group_col]])
  groups <- levels(data[[group_col]])
  k <- length(groups)
  
  if (k < 2) {
    stop("You need at least two groups for ANOVA.")
  }
  
  # Group-wise statistics
  n_i    <- numeric(k)
  mean_i <- numeric(k)
  var_i  <- numeric(k)
  
  for (i in seq_along(groups)) {
    subset_data <- data[data[[group_col]] == groups[i], response_col, drop = TRUE]
    n_i[i]      <- length(subset_data)
    mean_i[i]   <- mean(subset_data)
    var_i[i]    <- var(subset_data)
    
    if (n_i[i] < 2) {
      stop("Each group must have at least 2 observations.")
    }
    if (var_i[i] == 0) {
      stop("Group ", groups[i], " has zero variance; Welch's ANOVA is not well-defined.")
    }
  }
  
  # Welch weights
  w_i <- n_i / var_i
  W <- sum(w_i)  # sum of all weights
  
  # Weighted grand mean
  y_bar_w <- sum(w_i * mean_i) / W
  
  # Q: weighted between-group sum of squares
  Q <- sum(w_i * (mean_i - y_bar_w)^2)
  
  # Numerator df
  df_num <- k - 1
  
  # G: for Welch's correction factor
  G <- sum((w_i / W)^2 / (n_i - 1))
  
  # Welch-type denominator factor
  denom_factor <- 1 + (2 * (k - 1) / (k^2 - 1)) * G
  
  # Welch's F*
  F_star <- (Q / df_num) / denom_factor
  
  # Approx. denominator df
  df_den <- (k - 1)^2 / (3 * G)
  
  # p-value
  p_val <- 1 - pf(F_star, df1 = df_num, df2 = df_den)
  
  # ------------------------------------------------------------------------
  # Compute Weighted Total Sum of Squares (TSS_w)
  # TSS_w = Q + sum_i (n_i - 1)
  # This ensures Q <= TSS_w and yields an R^2 in [0, 1].
  # ------------------------------------------------------------------------
  n_total <- sum(n_i)
  TSS_w <- Q + (n_total - k)
  R2_welch_weighted <- Q / TSS_w
  
  # Return a structured result (like an 'htest' object)
  results <- list(
    statistic    = F_star,
    parameter    = c(df_num = df_num, df_den = df_den),
    p.value      = p_val,
    R2_welch     = R2_welch_weighted
  )
  
  return(results)
}