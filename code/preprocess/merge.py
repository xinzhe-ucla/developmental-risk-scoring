### merge.py ######################################################################################
# purpose: merge the developmental samples into a common h5ad object

### PREAMBLE ######################################################################################
# load in libraries:
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


# specify data:
sample_to_read = [
    '/u/project/cluo/heffel/BICAN/yzcl54',
    '/u/project/cluo/heffel/BICAN/yzcl73',
    '/u/project/cluo/heffel/BICAN/yzcl74',
    '/u/project/cluo/heffel/BICAN/yzcl75',
    '/u/project/cluo/heffel/BICAN/yzcl76',
    '/u/project/cluo/heffel/BICAN/yzcl78',
    '/u/project/cluo/heffel/BICAN/DLPFC/210111',
    '/u/project/cluo/heffel/BICAN/DLPFC/210224',
    '/u/project/cluo/heffel/BICAN/DLPFC/210316',
    '/u/project/cluo/heffel/BICAN/DLPFC/210505',
    '/u/project/cluo/heffel/BICAN/DLPFC/210528',
    '/u/project/cluo/heffel/BICAN/source_fastq/Lee_2019',
    '/u/project/cluo/heffel/BICAN/source_fastq/2019_7mo',
    '/u/project/cluo/heffel/BICAN/HPC/yzcl20',
    '/u/project/cluo/heffel/BICAN/HPC/yzcl23',
    '/u/project/cluo/heffel/BICAN/yzcl32',
    '/u/project/cluo/heffel/BICAN/yzcl22',
]

adatas = []
loaded_file = []

# load in the samples h5ad (obs::QC)
for sample in tqdm(sample_to_read):
    sample_id = re.sub(r'.*/', '', sample)
    pattern = os.path.join(sample, f'mcg*unnormalized_genebody.h5ad')
    matched_files = glob.glob(pattern)
    adata = sc.read_h5ad(matched_files[0])
    adata.obs['batch'] = f'sample_{sample}'
    adatas.append(adata)
    loaded_file.append(matched_files)

# load in the sample cell metadata:
meta_data = {}
for sample in tqdm(sample_to_read):
    sample_id = re.sub(r'.*/', '', sample)
    pattern = os.path.join(sample, f'metadata_*.tsv')
    matched_files = glob.glob(pattern)
    meta = pd.read_csv(matched_files[0], sep = '\t')
    meta['sample'] = sample_id
    meta_data[sample] = meta

# combine the samples
adata_combined = ad.concat(adatas, axis=0, label='sample', keys=[f'sample_{re.sub(r".*/", "", s)}' for s in sample_to_read]) # about 24G in memory
meta_combined = pd.concat(meta_data, axis = 0)

# seems like there are some samples that contains NAs, identify them and remove them:
nan_rows = np.isnan(adata_combined.X).any(axis=1)
cells_with_nans = adata_combined.obs[nan_rows]
print(f"{nan_rows.sum()} cells contain NaNs") # 9 cells contains nan

# filter them off:
adata_combined = adata_combined[~nan_rows].copy()

# change the gene names:
mg = mygene.MyGeneInfo()
ensg_ids = [gene_id.split('.')[0] for gene_id in adata_combined.var_names] # remove the suffixes
query_result = mg.querymany(ensg_ids, scopes="ensembl.gene", fields="symbol", species="human")

# create a dictionary on the ensg id to hgnc symbol
ensg_to_symbol = {
    item['query']: item.get('symbol', item['query'])  # Use ENSG if no symbol found
    for item in query_result if not item.get('notfound', False)
}

# Replace var_names with HGNC symbols where available
adata_combined.var_names = [ensg_to_symbol.get(gene_id, gene_id) for gene_id in ensg_ids]
adata_combined.var_names = make_unique(adata_combined.var_names.tolist())

# check if all genes are unique:
assert adata_combined.var_names.is_unique, "Gene names are still not unique!"

# output one that doesn't have any preprocessing:
adata_combined.write(f'/u/home/l/lixinzhe/project-cluo/data/{today}-combined-mcg-raw.h5ad')

### MET_SCDRS_PREPROCESS ##########################################################################
# Compute the gene variances:
print('filtering low variance genes: ')
gene_variances = pd.Series(adata_combined.X.var(axis=0), index=adata_combined.var_names)
percentile_5th = gene_variances.quantile(0.05)
variance_mask = gene_variances >= percentile_5th

# filter based on the variance mask:
adata_combined = adata_combined[:, variance_mask]

# finally flip the adata_combined:
adata_combined.X = 1 - adata_combined.X
adata_meta = adata_combined.obs

### OUTPUT ########################################################################################
adata_combined.write(f'/u/home/l/lixinzhe/project-cluo/data/{today}-combined-mcg-processed.h5ad')
adata_meta.to_csv(f'/u/home/l/lixinzhe/project-cluo/data/{today}-combined-meta.csv', index = True)
