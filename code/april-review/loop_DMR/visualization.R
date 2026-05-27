### visualize the ldsc result as a heatmap of age by cell type

### PREAMBLE ######################################################################################
# load in the packages for us to investigate the plot:
library(ggplot2)
library(dplyr)

# load in the data:
ldsc_dir = "/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/pairwise_with_fetal_brain/"
ldsc_results = list.files(ldsc_dir, pattern = '.results')

# for each of the results, read in the
all_dev_time = all_cell_types = {}
for (h_est in ldsc_results){
    #get meta data:
    dev_time = gsub('_dms2.*', '', gsub('methylpy_pair_', '', h_est))
    cell_type = 'MSN_DRD1-eccentric-CASZ1'
    
    # aggregate to get cell types:
    all_cell_types = cell_type
    all_dev_time = c(all_dev_time, dev_time)
}

# load in the h2 estimate:
table_columns = c('time_point', 'Prop._h2', 'Prop._h2_std_error', 'Enrichment', 'Enrichment_std_error', 'Enrichment_p', 'neg_log_p')
all_cell_types = unique(all_cell_types)
all_dev_time = unique(all_dev_time)
h2_df = data.frame(matrix(NA, nrow = length(all_dev_time), ncol = length(table_columns)))
colnames(h2_df) = table_columns
rownames(h2_df) = all_dev_time

# load in the matrix to fill out the h_est
for (h_est in ldsc_results){
    loaded = read.table(paste0(ldsc_dir, h_est), sep = '\t', header = TRUE)
    rownames(loaded) = loaded$Category
    #get meta data:
    dev_time = gsub('_dms2.*', '', gsub('methylpy_pair_', '', h_est))
    cell_type = 'MSN_DRD1-eccentric-CASZ1'
    
    # record:
    h2_df[dev_time, 'time_point'] = dev_time
    h2_df[dev_time, 'Prop._h2'] = loaded['L2_2', 'Prop._h2']
    h2_df[dev_time, 'Prop._h2_std_error'] = loaded['L2_2', 'Prop._h2_std_error']
    h2_df[dev_time, 'Enrichment'] = loaded['L2_2', 'Enrichment']
    h2_df[dev_time, 'Enrichment_std_error'] = loaded['L2_2', 'Enrichment_std_error']
    h2_df[dev_time, 'Enrichment_p'] = loaded['L2_2', 'Enrichment_p']
    h2_df[dev_time, 'neg_log_p'] = -log(loaded['L2_2', 'Enrichment_p'])
}
write.table(h2_df, sep = ',', file = paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-pairwise-with-fetal-brain-DNAase-h2-table.csv'))