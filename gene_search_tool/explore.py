import anndata
from tqdm import tqdm

if __name__ == "__main__":

    adata = anndata.read_h5ad('C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\zf_atlas_full_v4_release.h5ad')
    # options = ['10hpf', '12hpf', '14hpf', '16hpf', '19hpf', '1dpf', '2dpf', '3dpf', '5dpf', '10dpf']
    options = ['adaxial cell', 'common myeloid progenitor', 'ectodermal cell', 'floor plate', 'hatching gland cell', 'head mesenchyme', 'heart']
    data = adata[adata.obs['zebrafish_anatomy_ontology_class']== "notochord"]
    for i in tqdm(options):
        sample2_adata = adata[adata.obs['zebrafish_anatomy_ontology_class'] == i]  ## with first 1000 genes
        data = anndata.AnnData.concatenate(*[data, sample2_adata],
                                                      join="inner")  ## will automatically find the common genes to concatenate cells


    print(data.X.shape)
