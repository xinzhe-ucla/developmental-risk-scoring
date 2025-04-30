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

# import date:
today = date.today()
today = today.strftime("%Y-%m-%d")

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

# for each of the sample, load in the h5ad:
adata_combined = ad.concat(adatas, axis=0, label='sample', keys=[f'sample_{re.sub(r".*/", "", s)}' for s in sample_to_read]) # about 24G in memory
meta_combined = pd.concat(meta_data, axis = 0)
adata_meta = adata_combined.obs

### OUTPUT ########################################################################################
adata_combined.write(f'/u/home/l/lixinzhe/project-cluo/data/{today}-combined-mcg.h5ad')
adata_meta.to_csv(f'/u/home/l/lixinzhe/project-cluo/data/{today}-combined-meta.csv', index = True)
