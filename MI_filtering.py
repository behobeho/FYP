import numpy as np
import pandas as pd
import seaborn as sns

regime = 0
run = 1

counts = pd.read_csv(f"C:/Users/bella/OneDrive - Imperial College London/FYP/Code/BENCHMARKING/{regime}_example_data_noisy{run}.csv")
mi= pd.read_csv(f"C:/Users/bella/OneDrive - Imperial College London/FYP/Code/BENCHMARKING_mi_chi2/mi_{regime}_{run}.csv", index_col=0)
chi2 = pd.read_csv(f"C:/Users/bella/OneDrive - Imperial College London/FYP/Code/BENCHMARKING_mi_chi2/chi2_{regime}_{run}.csv", index_col=0)


n = mi.shape[0]

rows = []
genes = np.arange(1, 1200, 1)


# 1. Extract upper triangular entries (excluding diagonal)
for i in range(n):
    for j in range(i + 1, n):
        rows.append((genes[i], genes[j], mi.iloc[i, j], chi2.iloc[i, j]))

df = pd.DataFrame(rows, columns=["gene1", "gene2", "MI", "Chi2"])

# 2. Sort by descending MI
df = df.sort_values("MI", ascending=False).reset_index(drop=True)
print(df)

# 3. Greedy selection of top MI-connected genes until 100 unique genes
selected_genes = set()
selected_list = []

for g1, g2, mi, chi2 in df.itertuples(index=False):
    if len(selected_genes) >= 100:
        break

    if g1 not in selected_genes:
        selected_genes.add(g1)
        selected_list.append(g1)

    if len(selected_genes) >= 100:
        break

    if g2 not in selected_genes:
        selected_genes.add(g2)
        selected_list.append(g2)

# If you want it as a list instead of mixed order insertions:
selected_genes = list(selected_genes)

print(f"Selected {len(selected_genes)} genes")

counts = counts.drop(columns=["Unnamed: 0"])
counts = counts.T
print(counts)

filtered_counts = counts.iloc[selected_genes]
print(filtered_counts)

