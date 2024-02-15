import anndata
import numpy as np
import scipy

import matplotlib.pyplot as plt
if __name__ == "__main__":
    adata = anndata.read_h5ad("C:\\Users\\rohan\PycharmProjects\SciFaitEmory\data\zf_atlas_10hpf_v4_release.h5ad")



    vals = []
    # for i in range(32060):

    plt.hist(adata.X[:, adata.var_names.tolist().index("tppp3")].toarray().flatten())  ## plot the distribution of first gene
    plt.show()
    plt.hist(np.log(adata.X[:, 0] + 0.1))  ## make a log-transformed change and check what has changed
    plt.show()
