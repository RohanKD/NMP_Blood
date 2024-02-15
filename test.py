from bioinfokit import analys, visuz
import pandas as pd
if __name__=="__main__":
    df = pd.read_csv("C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\models\\features.txt")
    df = df.set_index(df.columns[0])
    df= df.head(100)
    visuz.gene_exp.hmap(df=df, rowclus=False, colclus=False, dim=(3, 6), tickfont=(6, 4))

import anndata
from utils import _utils

if __name__=="__main__":
    target_adata = anndata.read_h5ad("../data/scnet_mouse/Mouse_pFC.h5ad")  # test
    _utils._visualize_data(target_adata, output_dir="C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\models",prefix="traindata_vis")
"""
    featureCounts -a Danio_rerio.GRCz11.110.gtf.gz 
              -o gene_counts.txt 
              -T 4 
              -p 
              -B 
              -C 
              -s 0 
              *.bam
              
              """