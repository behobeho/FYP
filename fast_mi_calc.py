import numpy as np
import pandas as pd
from fast_mutual_information import *
import time

runs = np.arange(0, 14, 1)
regs = [80, 60, 40, 20, 0]

for run in runs:
    for reg in regs:
    
        start = time.perf_counter()
        
        alpha = np.load(f"/home/beho/nb_{reg}_{run}/alpha.npy")
        r = np.load(f"/home/beho/nb_{reg}_{run}/r.npy")
        mu = np.load(f"/home/beho/nb_{reg}_{run}/mu.npy")
        data = pd.read_csv(f"/home/beho/HVG_filtered/{reg}_HVG_{run}.csv")
        counts = data.to_numpy(dtype=np.int32)
        counts = counts.T
        print("alpha:", len(alpha), "r", len(r), "mu", (len(mu)))
        print(counts)


        mu = np.asarray(mu, dtype=np.float64).reshape(-1, 1)
        r = np.asarray(r, dtype=np.float64).reshape(-1, 1)
        alpha = np.asarray(alpha, dtype=np.float64).reshape(-1, 1)

        mi_matrix, chi2_matrix = mi_negative_binomial_zi(counts, mu, r, alpha, min_pop=50)
        np.savetxt(f"BENCHMARKING_mi_chi2/mi_{reg}_{run}.csv", mi_matrix, delimiter=',')
        np.savetxt(f"BENCHMARKING_mi_chi2/chi2_{reg}_{run}.csv", chi2_matrix, delimiter=',')
        
        elapsed = time.perf_counter() - start

        print(
            f"Regime={reg:2d}, Run={run:2d}, "

            f"Time={elapsed:.2f} s"
        )




