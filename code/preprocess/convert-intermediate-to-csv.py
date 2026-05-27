import pickle
import numpy as np

# load in intermediate:
with open('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged.h5ad', 'rb') as f:
    intermediate = pickle.load(f)

# save the matrix:
np.save("/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged_X_only", intermediate.X)
intermediate.obs.to_csv("/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged_cell_meta.csv")
intermediate.var.to_csv("/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged_gene_meta.csv")
