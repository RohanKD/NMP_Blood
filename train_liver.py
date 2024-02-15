import os
import anndata
import pandas as pd
import numpy as np
import tensorflow as tf
## problem of different versions of Anndata:
## it would read data into bytes (higher version) or string (lower version)..
from sklearn.preprocessing import OneHotEncoder

from utils import _utils


def main():
    # '''Load PBMC Zheng FACS-sorted data
    # '''
    # adata = anndata.read_mtx('C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\sham1_matrix.mtx\matrix.mtx')
    # genes = pd.read_csv("C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\sham1_matrix.mtx\\features.tsv")
    # adata.var["genes"] = genes[0].values
    # adata.var_names = adata.var["genes"]
    # adata.var_names_make_unique(join="-")
    # adata.var.index.name = None
    # cells = pd.read_csv("C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\sham1_matrix.mtx\\sham1_barcodes.tsv")
    # adata.obs["barcode"] = cells[0].values
    # adata.obs_names = adata.obs['barcode']
    # adata.obs_names_make_unique(join="-")
    # adata.obs.index.name = None
    #
    # adata = adata[:, adata.var_names.notnull()]
    # adata.var_names = [i.upper() for i in list(adata.var_names)]



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