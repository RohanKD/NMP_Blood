# import anndata
# import numpy as np
# import random
#
# adata = anndata.read_h5ad("C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\zf_atlas_16hpf_v4_release.h5ad")
# non_NMP_adata = adata[adata.obs_names[0]]
# non_NMP_index = [random.randint(0, 6000) for _ in range(40)]
# if len(non_NMP_index) > 0:
#     NMP_adata = adata[adata.obs_names[non_NMP_index[0]]]
#     for i in non_NMP_index[1:]:
#         NMP_adata = anndata.AnnData.concatenate(*[NMP_adata, adata[adata.obs_names[i]]],
#                                                 join="inner")
#
#
# def loc(gene_symbol):
#     return adata.var_names.tolist().index(gene_symbol)
#
#
# if __name__ == "__main__":
#     similarity = []
#     with open('genes.txt', "r") as f:
#         for line in f:
#             gene = line.split()[0]
#             similarity.append(np.sum(non_NMP_adata.X[:, loc(gene)]))
#     print("length")
#     print(len(similarity))
#
#     loc_dict = {}
#     for i, val in enumerate(similarity):
#         loc_dict[val] = non_NMP_adata.var_names[i]
#     similarity = sorted(similarity, reverse=True)
#     output = []
#     unique_genes = set()
#     for i in range(len(similarity)):
#         gene = loc_dict[similarity[i]]
#         if gene not in unique_genes:
#             output.append(str(gene) + " - " + str(similarity[i]))
#             unique_genes.add(gene)
#
#     with open('genes_.txt', 'w') as f:
#         for value in output:
#             f.write(str(value) + '\n')
import math
import random

import anndata
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
from tqdm import tqdm


def create_gene(adata, non_NMP_adata):
    # non_NMP_adata = adata[adata.obs_names[0]]
    # non_NMP_index = [random.randint(0, 6000) for _ in range(40)]
    # if len(non_NMP_index) > 0:
    #     nonNMP_adata = adata[adata.obs_names[non_NMP_index[0]]]
    #     for i in non_NMP_index[1:]:
    #         NMP_adata = anndata.AnnData.concatenate(*[NMP_adata, adata[adata.obs_names[i]]],
    #                                                 join="inner")
    # shape = (41, 32060)
    # iterate through NMP_adata
    similarity = []
    gene_names = adata.var_names.tolist()
    corresponding = []
    for i in gene_names:
        corresponding.append(i)
    for i in (tqdm(range(27771))):
        similarity.append(np.sum(adata.X[non_NMP_adata, i]))
    pairs = [(name, score) for name, score in zip(corresponding, similarity)]
    pairs = sorted(pairs, key=lambda tup: tup[1], reverse=True)
    loc = {}
    for pair in pairs:
        loc[pair[0]] = pair[1]
    # with open('genes.txt', 'w') as f:
    #
    #     f.write(str(pair[0]) + " - " + str(pair[1]) + "\n")
    return loc
