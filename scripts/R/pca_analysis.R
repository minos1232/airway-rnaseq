# =========================================================================
# Principal Component Analysis (PCA) of Airway RNA-seq Samples
# =========================================================================
# Filters out zero-variance genes, applies log2(count + 1) transformation,
# performs scaled PCA across all 8 samples, and visualizes sample clustering
# by dexamethasone treatment and donor identity.

suppressPackageStartupMessages(library(ggplot2))

results_dir <- "~/bioinformatics/airway-rnaseq/results"

# 1. Load Count Matrix & Sample Metadata
my_counts <- read.csv(file.path(results_dir, "gene_count_matrix.csv"), row.names = 1)
metadata  <- read.csv(file.path(results_dir, "metadata.csv"))

# 2. Filter Out Zero-Variance (Undetected) Genes & Log2-Transform
gene_variances     <- apply(my_counts, 1, var)
my_counts_filtered <- my_counts[gene_variances > 0, ]
log_counts         <- log2(my_counts_filtered + 1)

# 3. Run Scaled PCA (samples as rows, genes as columns)
pca_result <- prcomp(t(log_counts), scale. = TRUE)
cat("PCA Variance Summary:\n")
print(summary(pca_result)$importance[, 1:4])

# 4. Merge PCA Coordinates with Sample Metadata
pca_scores        <- as.data.frame(pca_result$x)
pca_scores$SRR    <- rownames(pca_scores)
pca_scores_merged <- merge(pca_scores, metadata, by = "SRR")

# 5. Plot & Save PCA Colored by Treatment and Donor
treatment_plot <- ggplot(pca_scores_merged, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA of Airway RNA-seq Samples (by Treatment)",
    x = "PC1",
    y = "PC2"
  )
ggsave(file.path(results_dir, "pca_treatment_plot.png"), plot = treatment_plot, width = 8, height = 6)

donor_plot <- ggplot(pca_scores_merged, aes(x = PC1, y = PC2, color = donor)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA of Airway RNA-seq Samples (by Donor)",
    x = "PC1",
    y = "PC2"
  )
ggsave(file.path(results_dir, "pca_donor_plot.png"), plot = donor_plot, width = 8, height = 6)
