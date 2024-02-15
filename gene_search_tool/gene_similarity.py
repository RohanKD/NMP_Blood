import anndata
from tqdm import tqdm


def compare_expr_values(reference, comparison):
    return (reference - comparison).power(2).sum()


def normalize_coeff(data):
    return [(x - min(data)) / (max(data) - min(data)) for x in data]


import pandas as pd
import scanpy as sc

if __name__ == "__main__":
    # load data into dataframe
    adata = sc.read_mtx("C:\\Users\\rohan\OneDrive\Documents\gene_analysis\\full_tailbud\matrix.mtx")
    adata.obs = pd.read_csv("C:\\Users\\rohan\OneDrive\Documents\gene_analysis\\full_tailbud\\barcodes.tsv",
                            index_col=0, sep="\t")
    adata.var = pd.read_csv("C:\\Users\\rohan\OneDrive\Documents\gene_analysis\\full_tailbud\features.tsv", index_col=0,
                            sep="\t")
    while True:
        # get user input for the gene to search for
        gene = str(input("search for a gene: "))

        # search using nested for loop
        try:
            reference_values = adata.X[:, adata.var_names.tolist().index(gene)]
            break
        except:
            print("gene not found in database")

    rankings = []

    for gene in tqdm(range(3200)):
        genes = adata.X[:, gene]
        # compare current gene array to reference_values of the searched gene
        sum = compare_expr_values(reference_values, genes)

        rankings.append(sum)

    ranked = sorted(rankings)
    for i in ranked:
        with open("genes.txt", "a") as f:
            f.write(str(adata.var_names[rankings.index(i)]) + " r = " + str(i) + "\n")

    # iterate through all genes for all cells O(300k)
    # compute differences for each cell for each gene and test for top 10 genes by least square difference
    # compute R through comparisons using benchmarks from the least similar and normalizing values to be between 0 and 1
    # output all statistics to a text file
    # add to SCnet (Zebrafish gene search)
    # pseudocode-ish
