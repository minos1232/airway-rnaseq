# 1. Load the full count matrix and log2-transform all genes
counts_mat <- read.csv("~/bioinformatics/airway-rnaseq/results/gene_count_matrix.csv", row.names = 1)
log_counts <- log2(counts_mat + 1) # Remove high-leverage influential points (compresses values to ~0-17)

# 2. Extract two untreated replicates (X = SRR1039508, Y = SRR1039512)
x <- log_counts$SRR1039508
y <- log_counts$SRR1039512

# 3. Calculate Correlation (r) and R-squared (R^2)
r_val <- cor(x, y)
r_sq  <- r_val^2
cat("Correlation (r):", round(r_val, 5), "\n") # Unlike print(), cat() doesn't print [1]
cat("R-squared (R^2):", round(r_sq, 5), "\n\n")

# 4. Fit Least-Squares Regression Line: y_hat = m*x + b
fit <- lm(SRR1039512 ~ SRR1039508, data = log_counts) # lm() fits least-squares line y_hat = mx + b
print(summary(fit)$coefficients) # Formula syntax: Y ~ X

# 5. Check the first 3 residuals (e = y - y_hat)
cat("\nFirst 3 gene residuals:\n")
print(round(head(residuals(fit), 3), 5))