source("kevin/Writeup7_dylan/welch_anova.R")

n <- 300
# Example data
set.seed(123)
df <- data.frame(
  measurement = c(rnorm(n, mean = 10, sd = 2),
                  rnorm(n, mean = 10.5, sd = 3),
                  rnorm(n, mean = 10, sd = 4)),
  group = rep(c("A", "B", "C"), each = n)
)

# Perform Welch’s ANOVA
res <- welch_anova(data       = df,
                   response_col = "measurement",
                   group_col    = "group")

# Check results
res

# Extract key values
res$statistic     # The Welch's F* statistic
res$parameter     # Numerator & denominator degrees of freedom
res$p.value       # p-value
res$R2_welch      # % of variance explained by the group


#########

data       = df
response_col = "measurement"
group_col    = "group"
# Build the formula for oneway.test
formula <- as.formula(paste(response_col, "~", group_col))

# Perform Welch's ANOVA (oneway.test with var.equal = FALSE)
test_result <- oneway.test(formula, data = data, var.equal = FALSE)

test_result$parameter
test_result$statistic
test_result$p.value
