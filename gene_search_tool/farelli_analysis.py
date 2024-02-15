import anndata
import pandas as pd
# import scanpy as sc
# import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
#
#
# def fix_arr(data):
#     index = 0
#     for i in data:
#         try:
#             float(i)
#         except:
#             data[index] = 0
#         index += 1
#
#
# def mtx_load():
#     mtx_prefix = "C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\\farell_1\Sur2023_counts.mtx"
#     adata = anndata.read_mtx(mtx_prefix)
#     genes = pd.read_csv('C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\\farell_1\GSE223922_Sur2023_counts_rows_genes.txt',
#                         header=None, sep='\t')
#     adata.var["genes"] = genes[0].values
#     adata.var_names = adata.var["genes"]
#     adata.var_names_make_unique(join="-")
#     adata.var.index.name = None
#     cells = pd.read_csv('C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\\farell_1\GSE223922_Sur2023_counts_cols_cells.txt',
#                         header=None, sep='\t')
#     adata.obs["barcode"] = cells[0].values
#     adata.obs_names = adata.obs['barcode']
#     adata.obs_names_make_unique(join="-")
#     adata.obs.index.name = None
#
#     adata = adata[:, adata.var_names.notnull()]
#     adata.var_names = [i.upper() for i in list(adata.var_names)]
#     return adata
#
if __name__ == "__main__":
    # super slow code do not do this
        # df = pd.read_csv("C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\\URD_Dropseq_Expression_Log2TPM.txt", sep='\t')
    # print(df.head())
    # # Convert the pandas DataFrame to an AnnData object
    # adata = anndata.AnnData(df)
    #
    # # Display the AnnData object
    # print(adata)
    # loading zebrafish embryo sample from daniocell
    mtx_prefix = "C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\\farell_1\matrix.mtx"
    adata = anndata.read_mtx(mtx_prefix).T
    genes = pd.read_csv(
        'C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\\farell_1\genes.tsv',
        header=None, sep='\t')
    adata.var["genes"] = genes[0].values
    adata.var_names = adata.var["genes"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None
    cells = pd.read_csv(
        'C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\\farell_1\\barcodes.tsv',
        header=None, sep='\t')
    adata.obs["barcode"] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    adata = adata[:, adata.var_names.notnull()]
    adata.var_names = [i.upper() for i in list(adata.var_names)]


#     # find cell ontology classes and filter, also find <10 hpf cells to compare sox2, sox3 expression with.
#
#
#
#     # data contains expression matrix, genes, and barcodes (individual cell names)
#     # var_names = genes
#     # obs_names = barcodes = cells
#     # adata.X = raw gene expr matrix
#     # adata.X[cell:gene]
#     # print("start")
#     # adata = anndata.read_h5ad("C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\\farell_1\zf_atlas_10hpf_v4_release.h5ad")
#     # index = 0
#     # # brachyury(tbxta) at index 22498 "high expr" is above 2.5
#     # # sox2 at index 30972
#     # # shape is (1251, 32060)
#     # # expr_values = []
#     # # expr_brachyury = []
#     # # expr_sox2 = []
#     # NMP_index = [34, 69, 149, 540, 758, 1037, 1205]
#     # non_NMP_index = [1, 2, 3, 4, 5, 6, 7]
#     # switch = 0
#     #
#     # for i in range(1251):
#     #     if adata.X[i, 30972] > 1.50:
#     #
#     #         if adata.X[i, 22498] > 1:
#     #             print(adata.obs_names[i])
#     #             NMP_index.append(i)
#     #             switch = 1
#     #     if switch == 0:
#     #         non_NMP_index.append(i)
#     #     switch = 0
#
#     # # checking if there's "high" expr for a gene in the NMP population compared to non NMP population
#     # # NMP population identified as containing relativey high sox2 and tbxta expression levels
#     # # correct here
#     #
#     # correlated_genes = []
#     # z_scores = []
#     # # test adata csr_matrix casting capabilities
#     #
#     # for gene in range(32060):
#     #     NMP_expr_levels = []
#     #     non_NMP_expr_levels = []
#     #     for cell in NMP_index:
#     #         # NMP_expr_levels.append
#     #         var = (float(adata.X[cell, gene]))
#     #         NMP_expr_levels.append(var)
#     #     for cell in non_NMP_index[0:20]:
#     #         # non_NMP_expr_levels.append
#     #         var = (float(adata.X[cell, gene]))
#     #         non_NMP_expr_levels.append(var)
#     #     # conduct two sample t-test for determining whether there is a statistically significant diff between the two cell populations
#     #     # for now a less robust test involving mean difference and std difference
#     #     fix_arr(NMP_expr_levels)
#     #     fix_arr(non_NMP_expr_levels)
#     #     # print('NMP_expr_levels')
#     #     # print(NMP_expr_levels)
#     #     mean_diff = abs(np.mean(NMP_expr_levels) - np.mean(non_NMP_expr_levels))
#     #     # print("mean_diff "+ str(mean_diff))
#     #     std_avg = (np.std(NMP_expr_levels) + np.std(non_NMP_expr_levels)) / 2
#     #     # print("std_avg "+str(std_avg))
#     #     z_score = mean_diff / std_avg
#     #     if z_score > 0.5:
#     #         correlated_genes.append(adata.var_names[gene])
#     #         # print(adata.var_names[gene])
#     #         z_scores.append(z_score)
#     #
#     # # print(correlated_genes)
#     # # print(z_scores)
#     # # print(len(z_scores))
#     # search = sorted(z_scores, reverse = True)
#     # highest_genes = []
#     # for score in search[0:10]:
#     #     highest_genes.append(correlated_genes[z_scores.index(score)])
#     # print(highest_genes)
