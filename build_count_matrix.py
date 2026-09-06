import pandas as pd

tx2gene = pd.read_csv("reference/tx2gene_from_fasta.tsv", sep="\t", header=None, names=["transcript_id", "gene_id"])

samples = ["SRR1039508", "SRR1039509", "SRR1039512", "SRR1039513",
           "SRR1039516", "SRR1039517", "SRR1039520", "SRR1039521"]

all_gene_counts = {}

for srr in samples:
    sample1 = pd.read_csv(f"quant/{srr}/quant.sf", sep="\t") # load this sample's quant.sf
    merged = sample1.merge(tx2gene, left_on="Name", right_on="transcript_id")    # merge with tx2gene
    gene_counts = merged.groupby("gene_id")["NumReads"].sum() # groupby gene_id, sum NumReads
    all_gene_counts[srr] = gene_counts # store result in all_gene_counts[srr]
    pass

count_matrix = pd.DataFrame(all_gene_counts)
count_matrix = count_matrix.round(0).astype(int) # round to nearest integer and convert to int
print(count_matrix.head())
print(count_matrix.shape)
count_matrix.to_csv("results/gene_count_matrix.csv")

data = {
    'SRR': ['SRR1039508', 'SRR1039509', 'SRR1039512', 'SRR1039513',
            'SRR1039516', 'SRR1039517', 'SRR1039520', 'SRR1039521'],
    'donor': ['N61311', 'N61311', 'N052611', 'N052611',
              'N080611', 'N080611', 'N061011', 'N061011'],
    'treatment': ['untreated', 'dexamethasone', 'untreated', 'dexamethasone',
                  'untreated', 'dexamethasone', 'untreated', 'dexamethasone'],
}
metadata = pd.DataFrame(data).set_index("SRR")
print(metadata)
metadata.to_csv("results/metadata.csv")