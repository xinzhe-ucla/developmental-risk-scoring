###########################################################################################
######          Visualize the percent significance over different runs               ######
###########################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# load in data:
mean_var_length_arcsine = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-arcsine-inv_std/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)

# load in the cell type:
meta = read.table('/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/metadata_10292025_subset.csv.gz', sep = ',', header=TRUE, row.names = 1)

# once this is loaded in, look at the proportion of significant cells in 
new_name = gsub('-0-0-0', '', rownames(mean_var_length_arcsine))
new_name = gsub('-1-0$', '', new_name)
new_name = gsub('-1$', '', new_name)
new_name = gsub('-1-0-0$', '', new_name)
rownames(mean_var_length_arcsine) = new_name

mean_var_length_arcsine$FDR = p.adjust(mean_var_length_arcsine$pval, method = 'fdr')

# load in previous result too:
mean_var_length_arcsine_previous = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)
new_name = gsub('-0-0-0', '', rownames(mean_var_length_arcsine_previous))
new_name = gsub('-1-0$', '', new_name)
new_name = gsub('-1$', '', new_name)
new_name = gsub('-1-0-0$', '', new_name)
rownames(mean_var_length_arcsine_previous) = new_name
mean_var_length_arcsine_previous$FDR = p.adjust(mean_var_length_arcsine_previous$pval, method = 'fdr')

# check the score:
print(sum(mean_var_length_arcsine$FDR < 0.1))
print(sum(mean_var_length_arcsine_previous$FDR < 0.1))

meta1 = meta[rownames(mean_var_length_arcsine),]
meta1$score = mean_var_length_arcsine[, 'zscore']
meta1$fdr = mean_var_length_arcsine[, 'FDR']
cell_type_result = meta1 %>% group_by(adjusted_L3) %>% summarise(mean_score = mean(score)) %>% as.data.frame()
age_result = meta1 %>% group_by(fine2_age_groups) %>% summarise(mean_score = mean(score)) %>% as.data.frame()

# check the score for old one:
meta2 = meta[rownames(mean_var_length_arcsine_previous),]
meta2$score = mean_var_length_arcsine_previous[, 'zscore']
meta2$fdr = mean_var_length_arcsine_previous[, 'FDR']

cell_type_result2 = meta2 %>% group_by(adjusted_L3) %>% summarise(mean_score = mean(score)) %>% as.data.frame()
age_result2 = meta2 %>% group_by(fine2_age_groups) %>% summarise(mean_score = mean(score)) %>% as.data.frame()

# check correlation;
cor(age_result$mean_score, age_result2$mean_score) # 0.9999995
cor(cell_type_result2$mean_score, cell_type_result$mean_score) # 0.9999927

### Check proportion also 
cell_type_result = meta1 %>%
  group_by(adjusted_L3) %>%
  summarise(
    prop_pass_fdr_0.1 = mean(fdr < 0.1, na.rm = TRUE)
  ) %>% data.frame()

cell_type_result2 = meta2 %>%
  group_by(adjusted_L3) %>%
  summarise(
    prop_pass_fdr_0.1 = mean(fdr < 0.1, na.rm = TRUE)
  ) %>% data.frame()
cor(cell_type_result2$prop_pass_fdr_0, cell_type_result$prop_pass_fdr_0) # 0.9997652
