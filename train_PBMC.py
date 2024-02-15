import os
import anndata
import pandas as pd
import numpy as np
import tensorflow as tf
## problem of different versions of Anndata:
## it would read data into bytes (higher version) or string (lower version)..
from sklearn.preprocessing import OneHotEncoder

from utils import _utils

_ZhengPBMC_celltype_categories = {
            'B cells': [b'b_cells', 'b_cells'],
            'CD4 T cells': [b'cd4_t_helper', b'naive_t', b'memory_t', b'regulatory_t',
                'cd4_t_helper', 'naive_t', 'memory_t', 'regulatory_t'],
            'CD8 T cells': [b'cytotoxic_t', b'naive_cytotoxic',
                'cytotoxic_t', 'naive_cytotoxic'],
            'NK cells': [b'cd56_nk', 'cd56_nk'],
            'Monocytes': [b'cd14_monocytes', 'cd14_monocytes'],
}

def main():
    '''Load PBMC Zheng FACS-sorted data
    '''
    adata = anndata.read_h5ad('C:\\Users\\rohan\PycharmProjects\SciFaitEmory\humanPBMC_scRNAseq\PBMC_Zheng_FACS\FACS_adata.h5ad')

    if isinstance(adata.obs_names[0], (bytes)):
        adata.obs_names = adata.obs_names.str.decode('utf-8')
    if isinstance(adata.var_names[0], (bytes)):
        adata.var_names = adata.var_names.str.decode('utf-8')
    print(adata.obs)
    print(adata)
    print(adata.var)

    data_dir = 'C:\\Users\\rohan\PycharmProjects\SciFaitEmory\humanPBMC_scRNAseq'
    plate_protocols = ["CEL-Seq2", "Smart-seq2"]
    metadata_df = pd.read_csv(data_dir + os.sep + "PBMC_protocols\metadata.txt",
                              header=0, sep="\t")
    ## plate-based data
    read_adata = anndata.read_mtx(data_dir + os.sep + "PBMC_protocols\counts.read.txt").T
    read_cells = pd.read_csv(data_dir + os.sep + "PBMC_protocols\cells.read.new.txt",
                             header=None)
    read_genes = pd.read_csv(data_dir + os.sep + "PBMC_protocols\genes.read.txt",
                             header=None)
    read_adata.var['gene_symbols'] = [x.split('_')[1] for x in read_genes[0].values]
    read_adata.var_names = read_adata.var['gene_symbols']
    read_adata.var_names_make_unique(join="-")  # make unique
    read_adata.var_names.name = None
    read_adata.obs['barcode'] = read_cells[0].values
    read_adata.obs_names = read_adata.obs['barcode']
    read_adata.obs_names_make_unique(join="-")  ## make unique
    read_adata.obs_names.name = None
    plate_metadata = metadata_df[metadata_df['Method'].isin(plate_protocols)]
    common_cells = set(plate_metadata['NAME']).intersection(set(read_adata.obs_names))
    common_cells = list(common_cells)
    read_adata = read_adata[common_cells]  # 1052 cells
    obs_df = read_adata.obs.merge(plate_metadata, how='left',
                                  left_index=True, right_on='NAME')
    obs_df.index = obs_df['barcode'].values
    read_adata.obs = obs_df

    ## umi-based data
    umi_adata = anndata.read_mtx(data_dir + os.sep + "PBMC_protocols\counts.umi.txt").T
    umi_cells = pd.read_csv(data_dir + os.sep + "PBMC_protocols\cells.umi.new.txt",
                            header=None)
    umi_genes = pd.read_csv(data_dir + os.sep + "PBMC_protocols\genes.umi.txt",
                            header=None)
    umi_adata.var['gene_symbols'] = [x.split('_')[1] for x in umi_genes[0].values]
    umi_adata.var_names = umi_adata.var['gene_symbols']
    umi_adata.var_names_make_unique(join="-")  # make unique
    umi_adata.var_names.name = None
    umi_adata.obs['barcode'] = umi_cells[0].values
    umi_adata.obs_names = umi_adata.obs['barcode']
    umi_adata.obs_names_make_unique(join="-")  ## make unique
    umi_adata.obs_names.name = None
    droplet_metadata = metadata_df[~metadata_df['Method'].isin(plate_protocols)]
    common_cells = set(droplet_metadata['NAME']).intersection(set(umi_adata.obs_names))
    common_cells = list(common_cells)
    umi_adata = umi_adata[common_cells]  # 29969 cells
    obs_df = umi_adata.obs.merge(droplet_metadata, how='left',
                                 left_index=True, right_on='NAME')
    obs_df.index = obs_df['barcode'].values
    umi_adata.obs = obs_df

    ## concatenate adata together
    adata = read_adata.concatenate(umi_adata, batch_key="protocol_type",
                                   batch_categories=['plate', 'droplet'])
    adata.obs.rename(columns={'CellType': 'cell.type'}, inplace=True)

    adata_obs = adata.obs
    adata_obs['Method'].replace(['10x Chromium (v2)', '10x Chromium (v2) A',
                                 '10x Chromium (v2) B'], '10x-v2', inplace=True)
    adata_obs['Method'].replace(['10x Chromium (v3)'], '10x-v3', inplace=True)
    adata.obs = adata_obs

    test_adata = adata
    train_adata = adata
    common_genes = set(train_adata.var_names).intersection(set(test_adata.var_names))
    train_adata = train_adata[:, list(common_genes)]
    test_adata = test_adata[:, list(common_genes)]
    train_adata = _utils._process_adata(train_adata, process_type='train')
    train_adata = _utils._select_feature(train_adata,
                                                  fs_method='F-test',
                                                  num_features=1000)  ## use F-test to select 1000 informative genes
    train_adata = _utils._scale_data(train_adata)  ## center-scale
    # _utils._visualize_data(combined_train_adata, output_dir=".",
    #                        prefix="traindata_vis")  ## visualize cell types with selected features on a low dimension (you might need to change some parameters to let them show all the cell labels)
    # train an MLP model on it
    MLP_DIMS = _utils.MLP_DIMS  ## get MLP structure from _utils.py
    x_train = _utils._extract_adata(train_adata)
    enc = OneHotEncoder(handle_unknown='ignore')
    y_train = enc.fit_transform(train_adata.obs[[_utils.Celltype_COLUMN]]).toarray()
    mlp = _utils._init_MLP(x_train, y_train, dims=MLP_DIMS,
                           seed=_utils.RANDOM_SEED)
    mlp.compile()
    mlp.fit(x_train, y_train)
    mlp.model.save('./trained_MLP')  ## save the model so that you can load and play with it
    encoders = dict()
    for idx, cat in enumerate(enc.categories_[0]):
        encoders[idx] = cat
    print("encoders")
    print(encoders)
    # set(adata.obs['Sample'])
    ## preprocess the test data and predict cell types
    combined_test_adata = _utils._process_adata(test_adata, process_type='test')
    combined_test_adata = combined_test_adata[:,
                          list(
                              train_adata.var_names)]  ## extract out the features selected in the training dataset
    test_data_mat = _utils._extract_adata(combined_test_adata)
    test_data_mat = (test_data_mat - np.array(train_adata.var['mean'])) / np.array(
        train_adata.var['std'])
    y_pred = tf.nn.softmax(mlp.model.predict(test_data_mat)).numpy()
    pred_celltypes = _utils._prob_to_label(y_pred, encoders)
    combined_test_adata.obs[_utils.PredCelltype_COLUMN] = pred_celltypes

    ## let us evaluate the performance --> luckily you will have the accuracy over 99%
    from sklearn.metrics import accuracy_score, adjusted_rand_score, f1_score

    print("Overall Accuracy:",
          accuracy_score(combined_test_adata.obs[_utils.Celltype_COLUMN], combined_test_adata.obs[
              _utils.PredCelltype_COLUMN]))
    print("ARI:",
          adjusted_rand_score(combined_test_adata.obs[_utils.Celltype_COLUMN], combined_test_adata.obs[
              _utils.PredCelltype_COLUMN]))
    print("Macro F1:",
          f1_score(combined_test_adata.obs[_utils.Celltype_COLUMN], combined_test_adata.obs[_utils.PredCelltype_COLUMN],
                   average='macro'))
    ## a lot more evaluation metrics can be found on sklearn.metrics and you can explore with them
    # save feature file
    model_save_dir = "C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\models"
    train_adata.var.loc[:, ['mean', 'std']].to_csv(
        "C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\models\\features.txt")
    print(train_adata.var)
    ## save enc information
    with open(model_save_dir + os.sep + "onehot_encoder.txt", 'w') as f:
        for idx, cat in enumerate(enc.categories_[0]):
            f.write('%d:%s\n' % (idx, cat))


if __name__=="__main__":
    main()