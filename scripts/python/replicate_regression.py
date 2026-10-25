"""
Biological Replicate Reproducibility via Linear Regression
------------------------------------------------------------------
Evaluates expression consistency across ~60,000 genes between two untreated
biological replicates (SRR1039508 vs. SRR1039512) using Pearson correlation (r),
coefficient of determination (R^2), least-squares linear regression, and residuals.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import linregress

# -------------------------------------------------------------------------
# 1. Load Count Matrix & Log2-Transform
# -------------------------------------------------------------------------
results_dir = os.path.expanduser("~/bioinformatics/airway-rnaseq/results")
counts_mat = pd.read_csv(f"{results_dir}/gene_count_matrix.csv", index_col=0)

# Log2(count + 1) prevents extreme housekeeping genes from acting as high-leverage points
log_counts = np.log2(counts_mat + 1)

# Select two untreated control replicates
x = log_counts["SRR1039508"]
y = log_counts["SRR1039512"]

# -------------------------------------------------------------------------
# 2. Fit Least-Squares Regression & Compute Correlation
# -------------------------------------------------------------------------
res = linregress(x, y)

print(f"Correlation (r): {res.rvalue:.5f}")
print(f"R-squared (R^2): {res.rvalue**2:.5f}")
print(f"Slope (m):       {res.slope:.5f}")
print(f"Intercept (b):   {res.intercept:.5f}")

# -------------------------------------------------------------------------
# 3. Compute Predicted Values (y_hat) and Residuals (e = y - y_hat)
# -------------------------------------------------------------------------
y_hat = res.slope * x + res.intercept
residuals = y - y_hat

print("\nFirst 3 gene residuals:")
print(residuals.head(3).round(5))