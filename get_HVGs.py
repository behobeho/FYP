from __future__ import annotations
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

data = "C:/Users/bella/OneDrive - Imperial College London/FYP/Code/concat_all_interpro.csv"
num_genes = 5000

def get_HVG(data, num_genes):

    #read and clean up CSV
    df = pd.read_csv(data, index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce')
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna(axis=0)

    #create anndata
    adata = sc.AnnData(df.T)
    adata.X = adata.X.astype(np.float32)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    #calculate 
    sc.pp.highly_variable_genes(adata, n_top_genes=num_genes)
    HVG = adata.var[['highly_variable']].reset_index()
    HVG.columns = ['accession', 'highly_variable']

    #append HVGs to original dataframe and create new df that only contains top HVGs
    concatHVG = df.merge(HVG, on="accession", how="outer")
    hvg_df = concatHVG[concatHVG.highly_variable==True].reset_index()
    hvg_df.drop(hvg_df.iloc[0, :])

    return hvg_df



print(get_HVG(data, num_genes))
get_HVG(data, num_genes).to_csv("interpro_top5000.csv")