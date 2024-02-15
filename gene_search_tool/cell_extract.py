import anndata
import numpy as np
import random
from tqdm import tqdm
from scipy.sparse import csr_matrix

random.seed(2023)


def generate_random_list(original_list, range):
    random_list = []
    while len(random_list) != len(original_list):
        random_number = random.randint(1, range)

        random_list.append(random_number)
    return random_list


# this whole file is trash. gotta work on fixing it

def extract(adata, genes, levels):
    adata_X = csr_matrix(adata.X)
    # for inputted genes find cells with higest expression of all 3 using gene expr distribtion for each gene
    # adata = anndata.read_h5ad("C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\zf_atlas_16hpf_v4_release.h5ad")
    index_ind = []
    gene_ind = []
    for i in genes:
        gene_ind.append(adata.var_names.tolist().index(i))
    needed = len(genes)
    length = adata.obs_names.shape[0]
    for cell in tqdm(range(length)):
        count = 0
        inc = -1
        for gene in gene_ind:
            inc += 1

            if (adata_X[cell, gene] > 1):
                count += 1

        # generalize this to any such n genes
        if count == needed:
            index_ind.append(cell)

    # distr = adata.X[:, adata.var_names.tolist().index(gene)]
    # top_indexes = get_top_indexes(distr)
    # index_ind.append(top_indexes)

    print("Length of index_ind:", len(index_ind))

    ext_indices = generate_random_list(index_ind, len(adata.obs_names))

    # common_cell_adata = adata[adata.obs_names[index_ind[0]]]
    # for i in index_ind[1:]:
    #     common_cell_adata = anndata.AnnData.concatenate(*[common_cell_adata, adata[adata.obs_names[i]]],
    #                                                     join="inner")
    #
    # non_common_adata = adata[adata.obs_names[ext_indices[0]]]
    #
    # for i in ext_indices[1:]:
    #     non_common_adata = anndata.AnnData.concatenate(*[non_common_adata, adata[adata.obs_names[i]]],
    #                                                    join="inner")

    return index_ind[0:500], ext_indices[0:500]
