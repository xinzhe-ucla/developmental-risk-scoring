import pandas as pd
from datetime import date
today = date.today()
import re
import os
import subprocess
from tqdm import tqdm

# load in the gtex:
gtex_files = [
    "/u/home/l/lixinzhe/project-geschwind/data/GTEX_eQTL/V7/Brain_Caudate_basal_ganglia.allpairs.txt.gz",
    "/u/home/l/lixinzhe/project-geschwind/data/GTEX_eQTL/V7/Brain_Nucleus_accumbens_basal_ganglia.allpairs.txt.gz",
    "/u/home/l/lixinzhe/project-geschwind/data/GTEX_eQTL/V7/Brain_Putamen_basal_ganglia.allpairs.txt.gz",
    '/u/home/l/lixinzhe/project-geschwind/data/GTEX_eQTL/V7/Brain_Frontal_Cortex_BA9.allpairs.txt.gz'
]

gtex_col = {}
for gtex_file in tqdm(gtex_files):
    tissue = re.sub('.*/', '', gtex_file)
    tissue = re.sub('.allpairs.txt.gz', '', tissue)
    
    gtex = pd.read_csv(
        gtex_file,
        sep = '\t'
        )

    gtex = gtex[gtex.pval_nominal < 0.05].copy()
    gtex['chromosome'] = [re.sub('_.*', '',f) for f in gtex.variant_id]
    gtex['stripped_id'] = gtex["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)

    # split the variant id:
    cols = gtex["variant_id"].str.split("_", expand=True)
    gtex["chromosome"] = cols[0]
    gtex["position"] = cols[1].astype(int)
    gtex["ref"] = cols[2]
    gtex["alt"] = cols[3]
    gtex_col[tissue] = gtex
    output_dir = '/u/home/l/lixinzhe/project-geschwind/data/GTEX_eQTL/V7_sig_only/'
    gtex.to_csv(f'{output_dir}{tissue}.marginal_significant.txt', sep = '\t', index = False)
