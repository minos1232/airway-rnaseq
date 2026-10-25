import os
import numpy as np
import pandas as pd

# Expand the '~' path just like R does automatically
results_dir = os.path.expanduser("~/bioinformatics/airway-rnaseq/results")

# ---------------------------------------------------------
# 1. Load & Filter DESeq2 Results (Equivalent to dplyr)
# ---------------------------------------------------------
res_df = pd.read_csv(f"{results_dir}/deseq2_results.csv")

# In Pandas, we filter using boolean masks combined with '&' (AND)
sig_genes = res_df[
    (res_df["padj"].notna())
    & (res_df["padj"] < 0.01)
    & (res_df["log2FoldChange"].abs() > 2)
].copy() #to make independant table in memory

# Equivalent to mutate(regulation = ifelse(...))
sig_genes["regulation"] = np.where(
    sig_genes["log2FoldChange"] > 0, "Up", "Down"
)

# Equivalent to arrange(padj)
sig_genes = sig_genes.sort_values("padj")

# Grab top 20 gene IDs as a Python list
top20_ids = sig_genes["gene_id"].head(20).tolist()

# ---------------------------------------------------------
# 2. Load Counts & Compute Row-Wise Z-Scores (Unit 4 Math)
# ---------------------------------------------------------
# index_col=0 is the exact equivalent of R's row.names = 1
counts_mat = pd.read_csv(f"{results_dir}/gene_count_matrix.csv", index_col=0)

# .loc[rows, columns] is the equivalent of R's counts_mat[top20_ids, ]
top20_counts = counts_mat.loc[top20_ids, :]

# Log2 transform: log2(count + 1)
top20_log = np.log2(top20_counts + 1)

# Calculate row means (axis=1 means across columns/samples for each row)
row_means = top20_log.mean(axis=1)

# Calculate row standard deviations (ddof=1 matches R's sample SD formula)
row_stds = top20_log.std(axis=1, ddof=1) #ddof=1 is the default for R's sd() function, which uses n-1 in the denominator

# Apply Unit 4 formula: Z = (x - mean) / std along rows (axis=0 aligns row index)
z_scores_py = top20_log.sub(row_means, axis=0).div(row_stds, axis=0)

# Compare the first 3 genes and first 4 samples to your R output!
print(z_scores_py.iloc[0:3, 0:4].round(7))