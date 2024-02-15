import anndata
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
from tqdm import tqdm


def create_genes(adata, NMP_adata):
    # # finding possible NMPs in the data
    # NMP_index = []
    # # prunes data to only include tail bud cells
    # # adata = adata[adata.obs['zebrafish_anatomy_ontology_class'] == "tail bud"]
    # # graph values to find determining point for high sox2 and tbxta
    # for i in range(6297):
    #     if adata.X[i, 30972] > .90:
    #         if adata.X[i, 22498] > .9:
    #             NMP_index.append(i)
    #
    # if len(NMP_index) > 0:
    #     NMP_adata = adata[adata.obs_names[NMP_index[0]]]
    #     for i in NMP_index[1:]:
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
        similarity.append(np.sum(adata.X[NMP_adata, i]))
    pairs = [(name, score) for name, score in zip(corresponding, similarity)]
    pairs = sorted(pairs, key=lambda tup: tup[1], reverse=True)
    loc = {}
    for pair in pairs:
        loc[pair[0]] = pair[1]
    # with open('genes.txt', 'w') as f:
    #
    #     f.write(str(pair[0]) + " - " + str(pair[1]) + "\n")
    return loc
