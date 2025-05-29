import argparse
import pandas as pd
import scanpy as sc
import scvi
import gzip

def main(path_to_data, path_to_result):
    # Read and prepare the data 
    dataset = pd.read_csv(path_to_data, sep="\t", index_col=0)
    dataset = dataset.T

    # Create an AnnData object from the DataFrame
    adata = sc.AnnData(dataset.values, obs=pd.DataFrame(index=dataset.index), var=pd.DataFrame(index=dataset.columns))

    # Set up the scVI model
    scvi.model.SCVI.setup_anndata(adata)

    # Initialize the scVI model
    model = scvi.model.SCVI(adata)

    # Train the scVI model
    model.train(max_epochs=400)

    # Perform imputation/denoising
    result = model.get_normalized_expression()

    # Save adjusted count matrix
    imputed_dataset = pd.DataFrame(result, index=adata.obs_names, columns=adata.var_names)
    with gzip.open(path_to_result, 'wt') as f:
        imputed_dataset.to_csv(f, sep="\t", index=True, header=True)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Perform scVI imputation on a gene counts matrix.")
    parser.add_argument("input_path", type=str, help="Path to the input tab-delimited gene counts matrix file.")
    parser.add_argument("output_path", type=str, help="Path to save the imputed counts matrix as a tab-delimited file.")

    args = parser.parse_args()

    # Run the main function with the specified input and output paths
    main(args.input_path, args.output_path)
