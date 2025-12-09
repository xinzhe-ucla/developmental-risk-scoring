### preprocess.py #################################################################################
# preprocess:
import os
import pandas as pd
import anndata
from anndata import AnnData
from datetime import date
import mygene

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

today = date.today()
print(today)

### PREAMBLE ######################################################################################
df=[]
for f in os.listdir('/u/project/jflint/heffel/Heng/hic_gene_score_h5ad/'):
    if 'h5' in f:
        tdf=pd.read_hdf('/u/project/jflint/heffel/Heng/hic_gene_score_h5ad/'+f, key='X')
        df.append(tdf)
        df=pd.concat(df)
        gdata=AnnData(df)

# save the gdata:
gdata.write(f'/u/project/cluo/lixinzhe/data/{today}_heng_loop_combined.h5ad')

# change the gene names:
mg = mygene.MyGeneInfo()
ensg_ids = [gene_id.split('.')[0] for gene_id in gdata.var_names] # remove the suffixes
query_result = mg.querymany(ensg_ids, scopes="ensembl.gene", fields="symbol", species="human")

# create a dictionary on the ensg id to hgnc symbol
ensg_to_symbol = {
    item['query']: item.get('symbol', item['query'])  # Use ENSG if no symbol found
    for item in query_result if not item.get('notfound', False)
}

# Replace var_names with HGNC symbols where available
gdata.var_names = [ensg_to_symbol.get(gene_id, gene_id) for gene_id in ensg_ids]
gdata.var_names = make_unique(gdata.var_names.tolist())

# check if all genes are unique:
assert gdata.var_names.is_unique, "Gene names are still not unique!"

# output one that doesn't have any preprocessing:
gdata.write(f'/u/home/l/lixinzhe/project-cluo/data/{today}-combined-3C-gene-score.h5ad')
