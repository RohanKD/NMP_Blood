#install.packages("Seurat")
#install.packages("Matrix")
library(Seurat)
library(Matrix)

# 
# # Load the barcodes
# barcodes_file <- "C:\\Users\\rohan\\OneDrive\\Documents\\gene_analysis\\barcodes.tsv" # Replace with the actual file path
# barcodes <- read.table(barcodes_file, header = FALSE, col.names = "BARCODE",sep = "\t")
# 
# # Load the features
# feature_col_names <- c("GENE", "SYMBOL", "GE")
# features_file <- "C:\\Users\\rohan\\OneDrive\\Documents\\gene_analysis\\features.tsv" # Replace with the actual file path
# features <- read.table(features_file, header = FALSE, col.names = feature_col_names,sep = "\t")
# 
# # Load the gene expression matrix in Matrix Market format
# matrix_file <- "C:\\Users\\rohan\\OneDrive\\Documents\\gene_analysis\\matrix.mtx" # Replace with the actual file path
# gene_expression <- readMM(matrix_file)
# # create seurat object
# seurat_object <- CreateSeuratObject(counts = gene_expression, project = "tailbudSEQ")
# 
# # add in features and barcodes
# seurat_object$barcodes <- barcodes$barcode
# seurat_object$features <- features$gene

# extra optional data normalization/analysis
#seurat_object <- NormalizeData(seurat_object)
#seurat_object <- FindVariableFeatures(seurat_object)
#seurat_object <- ScaleData(seurat_object)
#seurat_object <- RunPCA(seurat_object)
#seurat_object <- FindNeighbors(seurat_object)
#seurat_object <- FindClusters(seurat_object)

# load in h5ad file and convert to seurat object using 10x library
# seurat_object <- Read10X_h5("C:\\Users\\rohan\\OneDrive\\Documents\\gene_analysis\\tailbudSEQ.h5ad")
#install.packages("anndata")
library(anndata)
library(Seurat)

# Replace "path/to/your_file.h5ad" with the actual file path to your h5ad file
h5ad_file <- "C:\\Users\\rohan\\OneDrive\\Documents\\gene_analysis\\tailbudSEQ.h5ad"
library(reticulate)
library(anndata
        )
# python_executable <- Sys.which("python")
# print(python_executable)
use_python("C:\\Users\\rohan\\anaconda3\\python.exe")  # Replace with the path to your Python executable

# Import the Python 'anndata' package
anndata <- import("anndata")

# Replace "path/to/your_file.h5ad" with the actual file path to your h5ad file

# Read the h5ad file using read_h5ad function
h5ad <- anndata$read_h5ad(h5ad_file)

# Access the gene expression matrix and other data
gene_expression <- h5ad$X
features <- h5ad$var
barcodes <- rownames(h5ad$X)

# Create Seurat object using the gene expression matrix
library(Seurat)
seurat_object <- CreateSeuratObject(counts = gene_expression)

# Set features and barcodes metadata
seurat_object$features <- features
seurat_object$barcodes <- barcodes

# Continue with further Seurat processing (e.g., normalization, dimensionality reduction, clustering)
