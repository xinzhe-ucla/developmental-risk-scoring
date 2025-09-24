### updated_merge.py ##############################################################################
# load in packages:
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
import mygene
import polars as pl
import gc
import importlib.resources as pkg_resources

# import date:
today = date.today()
today = today.strftime("%Y-%m-%d")

# import function:
# identify and rename genes that are duplicated:
def make_unique(names):
    from collections import defaultdict
    counter = defaultdict(int)
    unique_names = []
    
    for name in names:
        count = counter[name]
        if count == 0:
            unique_names.append(name)
        else:
            unique_names.append(f"{name}.{count}")
        counter[name] += 1
    
    return unique_names

def ENSG_to_HGNC(adata, convert = False):
    """
    Change ENSG symbols to HGNC symbols
    CAUTION: collapse all genes that has the same converted hgnc id by summation
    
    Parameters
    ----------
    adata: ann.data object
        adata that contains ENSG as variable names
    convert: bool, optional
        should a new adata be returned with hgnc symbols. Collapse by summation to keep gene names unique
    """
    file_path = pkg_resources.files('spatial_invader.data').joinpath('2025-09-06-bioMart-ensg-to-hgnc-conversion-table.txt')
    conversion_file = pd.read_csv(file_path, sep = '\t')
    conversion_map = dict(zip(conversion_file['Gene stable ID'], conversion_file['Gene name']))
    
    # remove the suffix of the ensg name and map to hgnc:
    adata.var['ensg_stripped'] = adata.var_names.str.replace(r"\..*$", "", regex = True)
    adata.var['hgnc'] = adata.var['ensg_stripped'].map(conversion_map)
    adata.var.loc[adata.var['hgnc'].isna(), 'hgnc'] = adata.var_names[adata.var['hgnc'].isna()]
    
    if convert:
        # group by stripped names and sum columns
        df = pd.DataFrame(adata.X, columns=adata.var['hgnc'], index=adata.obs_names)
        df = df.groupby(df.columns, axis=1).mean()
        
        # rebuild AnnData
        adata = ad.AnnData(df, obs=adata.obs.copy())
        return adata
    else:
        return adata

# specify data:
sample_to_read = "/u/project/cluo/heffel/BICAN3/DATA/mcg_Exc-merge_unnormalized_genebody_blacklist_allgenes_88932x49073.h5ad"

# specify meta data path:
meta = pd.read_csv("/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz", index_col = 0)

# read in the data, convert to hgnc and output
adata = sc.read_h5ad(sample_to_read)
adata = ENSG_to_HGNC(adata, convert = True)
adata.write_h5ad("/u/project/cluo/lixinzhe/data/BICAN3/mcg_Exc-merge_unnormalized_genebody_blacklist_hgnc_88932x49073.h5ad")
