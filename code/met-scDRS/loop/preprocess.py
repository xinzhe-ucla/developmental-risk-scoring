### preprocess.py #################################################################################
# preprocess:
import os
import pandas as pd
import anndata
from anndata import AnnData
from datetime import date
import mygene
import numpy as np
import scanpy as sc

today = date.today()
print(today)

###########################################################################################
######                        Adapted to new Heng generated loops                    ######
###########################################################################################
# load in the methylation ratio
met_ratio = '/u/home/h/hex002/project-cluo/BICAN/loop_methylation/merged/all_groups_merged.h5ad'
loop_ratio = sc.read_h5ad(met_ratio)

# check data:
assert loop_ratio.var_names.is_unique, 'Gene names not unique'
assert np.sum(np.isnan(loop_ratio.X)) == 0, 'NAs in the methylation matrix'

# No preprocessing will be done because Heng has already preprocessed it well!

###########################################################################################
######                                  Get covariates                               ######
###########################################################################################
# load in the meta data:
meta = pd.read_csv('/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz', index_col = 0)
meta_reordered = meta.loc[loop_ratio.obs_names,:]

# write out the covaraite:
covariate = pd.DataFrame({
    "cell": meta_reordered.index,
    "const": 1,
    "global": meta_reordered['mCG/CG_global']
})

assert (covariate.cell == loop_ratio.obs_names).all(), 'not all covariate rownames match with loop ratio'

# write out the covariate table:
covariate.to_csv('/u/project/cluo/lixinzhe/data/BICAN3/Heng_all_groups_merged.cov', sep = '\t', index = False)
