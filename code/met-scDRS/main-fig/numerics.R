### Getting the numerics for publications:

###########################################################################################
######                       DDR1/2 significant proportion                           ######
###########################################################################################
# get the % of significant cells in DDR1 and DDR2:
significant_ratio = read.table(file = '/u/home/l/lixinzhe/project-geschwind/plot/2025-12-22-developmental-cell-type-brain-traits-proportion-adult-only.csv', sep = ',', check.names = F)
disease = 'PASS_Schizophrenia_Pardinas2018'
cell_type = c("Inh_DRD2-EPHA4", "Inh_DRD1-EPHA4", "Inh_DRD1-BACH2", "Inh_DRD2-BACH2", "Inh_DRD2-eccentric-CASZ1" ,"Inh_DRD1-eccentric-CASZ1")

significant_ratio[disease, cell_type]

## get the p value for all age group, LGE derived MS are more strongly enriched vs MGE and CGE derived inhibitory neurons.



###########################################################################################
######                                LGE vs CGE and MGE                             ######
###########################################################################################
# load in data
library(data.table)
library(tidyverse)
library(clusterProfiler);
library(ComplexHeatmap)
require(circlize);

meta = read.table('/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/metadata_10292025_subset.csv.gz', sep = ',', header=TRUE, row.names = 1)
scDRS.directory = '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/'

score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

mdd = data.frame(
    fread(
        file =  "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/PASS_MDD_Howard2019.score.gz",
        sep = '\t',
        header = TRUE,
        data.table = FALSE
        ),
        row.names = 1
    )
# format the meta data:
new_name = gsub('-0-0-0', '', rownames(mdd))
new_name = gsub('-1-0$', '', new_name)
new_name = gsub('-1$', '', new_name)
new_name = gsub('-1-0-0$', '', new_name)

# get the common cells:
common_cells = intersect(new_name, rownames(meta))
meta = meta[common_cells, ]

# label the '' as 'unknown'
meta$newL1[meta$newL1 == ''] = 'unknown'
meta$newL1[meta$newL1 == '?'] = 'unknown'
meta$newL2[meta$newL2 == ''] = 'unknown'
meta$newL2[meta$newL2 == '?'] = 'unknown'
meta$adjusted_L3[meta$adjusted_L3 == ''] = 'unknown'
meta$adjusted_L3[meta$adjusted_L3 == '?'] = 'unknown'
meta$adjusted_L3 = paste0(meta$newL1, '_', meta$adjusted_L3)

# read into the empty list:
for (result in score.files) {
    met_scdrs_result <- data.frame(
        fread(
        file = result,
        sep = '\t',
        header = TRUE,
        data.table = FALSE
        ),
        row.names = 1)
    rownames(met_scdrs_result) = new_name
    met_scdrs_result = met_scdrs_result[common_cells, ]
    met_scdrs_result$fdr = p.adjust(met_scdrs_result$pval, method = 'fdr')
    risk.score[[result]] = met_scdrs_result
    }

# simplify names:
list.names <- gsub(scDRS.directory, '', score.files);
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

# get the trait information
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);
brain_traits = trait.info$Trait_Identifier[trait.info$Category == 'brain']

### get p value for scz trait ###
score = risk.score[['PASS_Schizophrenia_Pardinas2018']]
stopifnot(rownames(score) == rownames(meta))

MSN_sig = sum(score[meta$newL2 == 'MSN', 'fdr'] < 0.01)
MSN_total = sum(meta$newL2 == 'MSN')
MSN_not_sig = MSN_total - MSN_sig

MGE_sig = sum(score[meta$newL2 == 'MGE', 'fdr'] < 0.01)
MGE_total = sum(meta$newL2 == 'MGE')
MGE_not_sig = MGE_total - MGE_sig

CGE_sig = sum(score[meta$newL2 == 'CGE', 'fdr'] < 0.01)
CGE_total = sum(meta$newL2 == 'CGE')
CGE_not_sig = CGE_total - CGE_sig

# make the contingency table:
contingency_table = data.frame(matrix(NA, 2, 2))
colnames(contingency_table) = c('MSN', 'CGE_and_MGE')
rownames(contingency_table) = c('num_sig', 'num_insig')

# perform chisquare test:
contingency_table['num_sig', 'MSN'] = MSN_sig
contingency_table['num_insig', 'MSN'] = MSN_not_sig
contingency_table['num_sig', 'CGE_and_MGE'] = CGE_sig + MGE_sig
contingency_table['num_insig', 'CGE_and_MGE'] = CGE_not_sig + MGE_not_sig

print(chisq.test(contingency_table))

### get p value for BP trait ###
score = risk.score[['PASS_BIP_Mullins2021']]
stopifnot(rownames(score) == rownames(meta))

MSN_sig = sum(score[meta$newL2 == 'MSN', 'fdr'] < 0.01)
MSN_total = sum(meta$newL2 == 'MSN')
MSN_not_sig = MSN_total - MSN_sig

MGE_sig = sum(score[meta$newL2 == 'MGE', 'fdr'] < 0.01)
MGE_total = sum(meta$newL2 == 'MGE')
MGE_not_sig = MGE_total - MGE_sig

CGE_sig = sum(score[meta$newL2 == 'CGE', 'fdr'] < 0.01)
CGE_total = sum(meta$newL2 == 'CGE')
CGE_not_sig = CGE_total - CGE_sig

# make the contingency table:
contingency_table = data.frame(matrix(NA, 2, 2))
colnames(contingency_table) = c('MSN', 'CGE_and_MGE')
rownames(contingency_table) = c('num_sig', 'num_insig')

# perform chisquare test:
contingency_table['num_sig', 'MSN'] = MSN_sig
contingency_table['num_insig', 'MSN'] = MSN_not_sig
contingency_table['num_sig', 'CGE_and_MGE'] = CGE_sig + MGE_sig
contingency_table['num_insig', 'CGE_and_MGE'] = CGE_not_sig + MGE_not_sig

print(chisq.test(contingency_table))

