import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import gzip
from scbig import run_scbig, setup_seed
import torch

def main(path_to_data, path_to_result):
    setup_seed(100)
    torch.cuda.is_available = lambda: False
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    # Read and prepare the data 
    dataset = pd.read_csv(path_to_data, sep="\t", index_col=0)
    dataset = dataset.T

    # Create an AnnData object from the DataFrame
    adata = sc.AnnData(dataset.values, obs=pd.DataFrame(index=dataset.index), var=pd.DataFrame(index=dataset.columns))
    
    # Preprocess data
    adata.raw = adata.copy()
    sc.pp.calculate_qc_metrics(adata, inplace=True, log1p=False, use_raw=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    # Calculate the library size factor
    adata.obs['cs_factor'] = adata.obs.total_counts / np.median(adata.obs.total_counts)
    # Log Normalization
    sc.pp.log1p(adata)
    # Calculate the gene size factor
    adata.var['gs_factor'] = np.max(adata.X, axis=0, keepdims=True).reshape(-1)

    # Train the scBiG model
    adata, record = run_scbig(adata, impute=True, return_all=True)
    
    # Log-normalization
    adata.X = adata.obsm['imputed']
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Save adjusted count matrix
    imputed_dataset = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    imputed_dataset = imputed_dataset.T
    with gzip.open(path_to_result, 'wt') as f:
        imputed_dataset.to_csv(f, sep="\t", index=True, header=True)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Perform scBiG imputation on a normalized gene counts matrix.")
    parser.add_argument("input_path", type=str, help="Path to the input tab-delimited gene counts matrix file.")
    parser.add_argument("output_path", type=str, help="Path to save the imputed counts matrix as a tab-delimited file.")

    args = parser.parse_args()

    # Run the main function with the specified input and output paths
    main(args.input_path, args.output_path)
