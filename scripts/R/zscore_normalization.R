suppressPackageStartupMessages(library(tidyverse))

# Loading DESeq2 results table (read_csv from readr is faster than base read.csv)
res_df <- read_csv(
  "~/bioinformatics/airway-rnaseq/results/deseq2_results.csv",
  show_col_types = FALSE
)

# Filter for significant genes & add a direction column
sig_genes <- res_df %>%
  # Remove low-count genes (NA padj), keep FDR < 0.01, and |log2FC| > 2
  filter(!is.na(padj), padj < 0.01, abs(log2FoldChange) > 2) %>%
  # Create new column called regulation
  mutate(regulation = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
  # Sort rows by padj (ascending by default)
  arrange(padj)

# Grab top 20 gene IDs from the gene_id column using the $ operator
top20_ids <- head(sig_genes$gene_id, 20)

# 1. Load the count matrix using column 1 (gene_id) as the row labels
counts_mat <- read.csv(
  "~/bioinformatics/airway-rnaseq/results/gene_count_matrix.csv",
  row.names = 1
)

# 2. Slice the matrix to keep ONLY the 20 rows matching top20_ids
top20_counts <- counts_mat[top20_ids, ]

# 3. Log2-transform counts: log2(count + 1)
top20_log <- log2(top20_counts + 1)

# 4. Calculate row-wise Z-scores
# scale() calculates mean and SD vertically (column-wise), so we transpose
# before scaling and transpose back afterward to restore genes as rows
z_scores_r <- t(scale(t(top20_log)))

# Check that every gene's Z-score row mean is 0 and row SD is 1
print(round(rowMeans(z_scores_r), 4)) # Round to 4 decimal places
print(apply(z_scores_r, 1, sd))       # Base R lacks rowSds(), so apply sd() across rows (MARGIN = 1)

# Inspect the first 3 genes and first 4 samples
print(round(z_scores_r[1:3, 1:4], 7))