### meta.py #######################################################################################
# purpose: getting the meta data for cells:

### PREAMBLE ######################################################################################
# load in libraries
import pandas as pd
import scanpy as sc
import sys
import re
import glob
import os
import anndata as ad
from tqdm import tqdm
from datetime import date
import numpy as np

# import date:
today = date.today()
today = today.strftime("%Y-%m-%d")

# load in the h5ad file:
normalized = sc.read_h5ad('/u/project/cluo/heffel/BICAN/JOINT/hold_adata_04282025.h5ad')

# get the normalized umap:
meta_data = normalized.obs
umap_df = pd.DataFrame(
    normalized.obsm["X_umap"],
    columns=["UMAP_1", "UMAP_2"],
    index=normalized.obs_names
)

combined_df = pd.concat([meta_data, umap_df], axis=1)

# output this combined_df:
combined_df.to_csv(f'/u/home/l/lixinzhe/project-cluo/data/{today}-combined-meta-QCed.csv', index = True)