library(reticulate)
# Conda env for scVI
# reticulate::conda_install("scvi-env", c("scvi-tools", "scanpy"), python_version = 3.10, pip = TRUE)
use_condaenv("scvi-env")

# Command line arguments
args <- commandArgs(TRUE)
path_to_data <- args[1] # Raw un-normalized counts
path_to_result <- args[2] # Imputed counts

# path_to_data <- "results/2_k562_analysis/gene_expression_filtered/sample_002.raw_counts.txt.gz"
# path_to_result <- "results/2_k562_analysis/gene_expression_imputed/sample_004.scvi.txt.gz"

# Import python packages with reticulate
pd <- import("pandas",  convert=FALSE)
np <- import("numpy",  convert=FALSE)
sc <- import("scanpy", convert=FALSE)
scvi <- import("scvi",  convert=FALSE)
gzip <- import("gzip",  convert=FALSE)

print("packages loaded ok")

# Read and prepare the data 
dataset = pd$read_csv(path_to_data, sep="\t", index_col=0)
dataset = dataset$T

print("dataset read ok")

# Create an AnnData object from the DataFrame
adata = sc$AnnData(dataset$values, obs=pd$DataFrame(index=dataset$index), var=pd$DataFrame(index=dataset$columns))
adata
print("anndata created ok")

# Set up the scVI model
scvi$model$SCVI$setup_anndata(adata)

print("model set up ok")


# Initialize the scVI model
model = scvi$model$SCVI(adata)

print("model initializeds ok")

# Train the scVI model
model$train()

# Perform imputation/denoising
result = model$get_normalized_expression(adata, library_size=1e4)

# Log-transform adjusted count matrix
imputed_dataset = np$log1p(result$T)

# Save adjusted count matrix
with(gzip$open(path_to_result, "wt") %as% f, {
  imputed_dataset$to_csv(f, sep="\t", index=TRUE, header=TRUE)
})


