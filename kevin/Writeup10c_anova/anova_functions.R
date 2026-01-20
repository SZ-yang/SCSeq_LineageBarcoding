.anova_percentage <- function(x, y){
  stopifnot(!any(is.na(y)))
  
  y <- as.numeric(y)
  total_std <- sum((y - mean(y))^2)
  
  if(is.factor(x)){
    across_lineage_std <- sum(sapply(levels(x), function(level){
      idx <- which(x == level)
      mean_val <- mean(y[idx])
      length(idx) * (mean_val - mean(y))^2 
    }))
    
    explained_variance <- across_lineage_std/total_std
    
    return(explained_variance)
  } else{
    stopifnot(is.numeric(x))
    
    tmp <- data.frame(y = y, x = x)
    lm_res <- stats::lm(y ~ ., data = tmp)
    summary_res <- summary(lm_res)
    
    return(summary_res$r.squared)
  } 
}


