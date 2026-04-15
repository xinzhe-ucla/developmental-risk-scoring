###########################################################################################
######          Visualize the percent significance over different runs               ######
###########################################################################################
library(ggplot2)
library(dplyr)

# load in data:
mean_var_logit = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-logit-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)
mean_var_arcsine = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-arcsine-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)
mean_var_library = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-library-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)

mean_var_length_logit = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-logit-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)
mean_var_length_arcsine = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-arcsine-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)
mean_var_length_library = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-library-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)

# load in the cell type:
meta = read.table('/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/metadata_10292025_subset.csv.gz', sep = ',', header=TRUE, row.names = 1)

# once this is loaded in, look at the proportion of significant cells in 