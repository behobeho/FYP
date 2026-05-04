#!/usr/bin/env python3

import argparse
import os
import pandas as pd
"""
def process(input_tsv):

    # 1. Read the TSV into a DataFrame (no header assumed)
    df = pd.read_csv(input_tsv, sep='\t', header=0)

    
    # 2. Extract the 5th column (index 4) containing "accession"
    sig_acc = df.iloc[:, 11]

    # 3. Count each accession and build a new DataFrame
    counts_df = (
        sig_acc
        .value_counts()
        .reset_index()
        .rename(columns={'index': 'interpro_id', 4: 'count'})
    ).sort_values('count', ascending = False)

    # 4. Construct output filename: replace .tsv with _counts.csv
    base, _ = os.path.splitext(input_tsv)
    output_csv = f"{base}_counts.csv"


    #5. Save to CSV
    counts_df.to_csv(output_csv, index=False)

    print(f"Saved counts to {output_csv}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Count occurrences of the 'Signature accession' column in an InterProScan TSV file"
    )
    parser.add_argument(
        'input_tsv',
        help='Path to the InterProScan output file (TSV format)'
    )
    args = parser.parse_args()

    process(args.input_tsv)
    
"""
def concat_data(counts_dir): #give function a directory of .tsv files and it concatenates them all to a single csv

    files = sorted([
        os.path.join(counts_dir, f) 
        for f in os.listdir(counts_dir) 
        if f.endswith(".csv")
    ])

    if not files:
        raise ValueError(f"no .tsv files found in {counts_dir}")
    
    df_master = None

    for file in files:

        df_new = pd.read_csv(file, header=None, names=["count", "accession", "description"])
        df_new = df_new.drop("description", axis=1)
        
        sample = os.path.basename(file).replace("_MERGED_FASTQ_InterPro.csv", "")

        df_new = df_new.rename(columns={"count": sample})

        print(df_new.columns)

        if df_master is None:

            df_master = df_new

        else:
            df_master = df_master.merge(df_new[["accession", sample]], on="accession", how="outer")

    df_master = df_master.fillna(0)

    return df_master


counts = "C:/Users/bella/OneDrive - Imperial College London/FYP/Code/interpro"
concat_file = concat_data(counts)

concat_file.to_csv("concat_all_interpro.csv", index=False)
