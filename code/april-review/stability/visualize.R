###########################################################################################
######          Visualize the percent significance over different runs               ######
###########################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)


# load in data:
mean_var_logit = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-logit-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)
mean_var_arcsine = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-arcsine-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)
mean_var_library = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-library-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)

mean_var_length_logit = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-logit-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)
mean_var_length_arcsine = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-arcsine-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)
mean_var_length_library = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-library-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)

# load in the cell type:
# meta = read.table('/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/metadata_10292025_subset.csv.gz', sep = ',', header=TRUE, row.names = 1)
meta = read.table('/u/project/cluo/heffel/BICAN3/REVISION/metadata_passQC_05212026.tsv.gz', sep = '\t', header=TRUE, row.names = 1)

# label the '' as 'unknown'
meta$newL1 = meta$L1
meta$newL2 = meta$L2
meta$adjusted_L3 = meta$L3
meta$fine2_age_groups = meta$age_group

# edit the unknowns:
meta$newL1[meta$newL1 == ''] = 'unknown'
meta$newL1[meta$newL1 == '?'] = 'unknown'
meta$newL2[meta$newL2 == ''] = 'unknown'
meta$newL2[meta$newL2 == '?'] = 'unknown'
meta$adjusted_L3[meta$adjusted_L3 == ''] = 'unknown'
meta$adjusted_L3[meta$adjusted_L3 == '?'] = 'unknown'
meta$adjusted_L3 = paste0(meta$newL1, '_', meta$adjusted_L3)

# once this is loaded in, look at the proportion of significant cells in 
stopifnot(rownames(mean_var_logit)== rownames(mean_var_length_logit))
new_name = gsub('-0-0-0', '', rownames(mean_var_logit))
new_name = gsub('-1-0$', '', new_name)
new_name = gsub('-1$', '', new_name)
new_name = gsub('-1-0-0$', '', new_name)
rownames(mean_var_logit) = rownames(mean_var_arcsine) = rownames(mean_var_library) = new_name
rownames(mean_var_length_logit) = rownames(mean_var_length_arcsine) = rownames(mean_var_length_library) = new_name

# perform fdr correction:
mean_var_logit$FDR = p.adjust(mean_var_logit$pval, method = 'fdr')
mean_var_arcsine$FDR = p.adjust(mean_var_arcsine$pval, method = 'fdr')
mean_var_library$FDR = p.adjust(mean_var_library$pval, method = 'fdr')

mean_var_length_logit$FDR = p.adjust(mean_var_length_logit$pval, method = 'fdr')
mean_var_length_arcsine$FDR = p.adjust(mean_var_length_arcsine$pval, method = 'fdr')
mean_var_length_library$FDR = p.adjust(mean_var_length_library$pval, method = 'fdr')


###########################################################################################
######       Check proportion of significance at DRD cell types                      ######
###########################################################################################
cell_types = c(
    'Inh_DRD1-BACH2',
    'Inh_DRD1-eccentric-CASZ1',
    'Inh_DRD1-EPHA4',
    'Inh_DRD2-BACH2',
    'Inh_DRD2-eccentric-CASZ1',
    'Inh_DRD2-EPHA4'
    )

# for each cell type in cell types interest, query the number of significant cells in cell types:
result_df = data.frame(matrix(NA, ncol = 6, nrow = length(cell_types)))
colnames(result_df) = cell_types
rownames(result_df) = c(
    'mean_var_logit', 
    'mean_var_arcsine',
    'mean_var_library',
    'mean_var_length_logit',
    'mean_var_length_arcsine',
    'mean_var_length_library'
    )


for (cell_type in cell_types){
    # get cell in cell type:
    cell_in_cell_type = rownames(meta)[meta$adjusted_L3 == cell_type]
    cell_in_cell_type = intersect(new_name, cell_in_cell_type)
    
    # compute proportoin that are significant:
    result_df['mean_var_logit', cell_type] = sum(mean_var_logit[cell_in_cell_type, 'FDR'] < 0.1) / length(cell_in_cell_type)
    result_df['mean_var_arcsine', cell_type] = sum(mean_var_arcsine[cell_in_cell_type, 'FDR'] < 0.1) / length(cell_in_cell_type)
    result_df['mean_var_library', cell_type] = sum(mean_var_library[cell_in_cell_type, 'FDR'] < 0.1) / length(cell_in_cell_type)
    result_df['mean_var_length_logit', cell_type] = sum(mean_var_length_logit[cell_in_cell_type, 'FDR'] < 0.1) / length(cell_in_cell_type)
    result_df['mean_var_length_arcsine', cell_type] = sum(mean_var_length_arcsine[cell_in_cell_type, 'FDR'] < 0.1) / length(cell_in_cell_type)
    result_df['mean_var_length_library', cell_type] = sum(mean_var_length_library[cell_in_cell_type, 'FDR'] < 0.1) / length(cell_in_cell_type)
}

###########################################################################################
######                                        visualize                              ######
###########################################################################################

df_long <- result_df %>%
  rownames_to_column("method") %>%
  pivot_longer(
    cols = -method,
    names_to = "celltype",
    values_to = "value"
  )
  
df_long = df_long %>% data.frame()
df_mean_var_length = df_long[grep('mean_var_length', df_long$method), ]
df_mean_var = df_long[!grepl("mean_var_length", df_long$method), ]

# plot bar plot:
gplot = ggplot(df_mean_var, aes(x = celltype, y = value, fill = method)) +
  geom_col(position = "dodge") +
  labs(
    x = "Cell type",
    y = "percent significant",
    fill = "Method"
  ) +
  scale_fill_manual(values = c(
    "mean_var_logit" = "#bebada",
    "mean_var_arcsine" = "#fb8072",
    "mean_var_library" = "#80b1d3"
  )) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# use the measured width and height for drawing:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-mean-var-robustness-comparison.pdf')
pdf(
    file = output.path,
    width = 5,
    height = 5
    );
print(gplot)
dev.off();


# plot bar plot:
gplot = ggplot(df_mean_var_length, aes(x = celltype, y = value, fill = method)) +
  geom_col(position = "dodge") +
  labs(
    x = "Cell type",
    y = "percent significant",
    fill = "Method"
  ) +
  scale_fill_manual(values = c(
    "mean_var_length_logit" = "#bebada",
    "mean_var_length_arcsine" = "#fb8072",
    "mean_var_length_library" = "#80b1d3"
  )) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# use the measured width and height for drawing:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-mean-var-length-robustness-comparison.pdf')
pdf(
    file = output.path,
    width = 5,
    height = 5
    );
print(gplot)
dev.off();
