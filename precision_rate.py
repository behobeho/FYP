import pandas as pd

# Load adjacency matrices
grn = pd.read_csv(r"C:\Users\bella\OneDrive - Imperial College London\FYP\Code\GRN_sim_adjacency_matrix.csv", index_col=0)
simulated = pd.read_csv(r"C:\Users\bella\0.001_sim_adjacency_matrix.csv", index_col=0)

print(grn.shape)
print(simulated.shape)

print(len(grn.index))
print(len(simulated.index))

# Replace row/column labels
simulated.columns = grn.columns
simulated.index = grn.columns

# Save corrected matrix
simulated.to_csv("sim_adjacency_matrix_fixed.csv")

print(simulated.head())
# Ensure matrices are aligned
grn = grn.loc[simulated.index, simulated.columns]

# Initialise counters
TP = 0
FP = 0

# Iterate through every cell
for row in grn.index:
    for col in grn.columns:

        grn_value = grn.loc[row, col]
        sim_value = simulated.loc[row, col]

        # True Positive
        if grn_value == 1 and sim_value == 1:
            TP += 1
            print(row, col)

        # False Positive
        elif grn_value == 0 and sim_value == 1:
            FP += 1

# Calculate precision
if (TP + FP) > 0:
    precision = TP / (TP + FP)
else:
    precision = 0

# Print results
print("True Positives:", TP)
print("False Positives:", FP)
print("Precision:", precision)