### load in data ##################################################################################
# load in data
### look at only adults

## plot out the disease x all cell types in adult only to find DRD tye enriched results
# load in the data:
library(data.table)
library(tidyverse)
library(clusterProfiler);
library(ComplexHeatmap)
require(circlize);

# load in umap:
bican_umap = read.csv('/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv', row.names = 1)

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

### look at the adult only data:
adult_cells <- rownames(meta)[meta$fine2_age_groups == 'adult']

# get the proportion of cells that are significant in each of the disease x cell type
meta_adult = meta[adult_cells, ]

bican_umap = bican_umap[rownames(meta), ]
meta = cbind(meta, bican_umap)

###########################################################################################
######                                    make figure D                              ######
###########################################################################################
disease = 'PASS_Schizophrenia_Pardinas2018'
risk_score = risk.score[[disease]][, 'zscore']
fdr = risk.score[[disease]][, 'fdr']
meta$risk_score = risk_score
meta$fdr = fdr

# for schizophrenia at adult, look at the box plot for these cell types:
pvalb_cell = rownames(meta)[meta$adjusted_L3 == 'Inh_PVALB']
lamp5_cell = rownames(meta)[meta$adjusted_L3 == 'Inh_LAMP5']
l4_rorb = rownames(meta)[meta$adjusted_L3 == 'Exc_L4-RORB']

# also look at the DRD1-eccentric-CASZ1
drd1_eccentric_cell = rownames(meta)[meta$adjusted_L3 == 'Inh_DRD1-eccentric-CASZ1']

plot_df = meta[c(pvalb_cell, lamp5_cell, l4_rorb, drd1_eccentric_cell), ]
plot_df$cell_type = factor(plot_df$adjusted_L3, levels = c('Inh_DRD1-eccentric-CASZ1', 'Exc_L4-RORB', 'Inh_PVALB', 'Inh_LAMP5'))
plot_df$age_group = factor(plot_df$fine2_age_groups, levels = c('2T', '3T', '1m', '4-7m', 'adult'))

gplot <- ggplot(plot_df, aes(x = age_group, y = risk_score, fill = cell_type)) +
    geom_boxplot() +
    scale_fill_manual(
        name = "L3",  # legend title
        values = c(
        'Exc_L4-RORB' = '#fdc086',
        'Inh_LAMP5' = '#fc8d62',
        'Inh_PVALB' = '#8da0cb',
        'Inh_DRD1-eccentric-CASZ1' = '#66c2a5'
        ),
        labels = c(
        "Exc_L4-RORB" = "L4-RORB",
        "Inh_LAMP5"  = "LAMP5",
        "Inh_PVALB"   = "PVALB",
        "Inh_DRD1-eccentric-CASZ1" = "DRD1\neccentric"
        )
    ) +
    xlab('Cell type') +
    ylab('met-scDRS') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size = 20))

output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-SCZ-L3-met-scdrs-by-age-box-plot.png')
png(
    filename = output.path,
    width = 7,
    height = 5,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

###########################################################################################
######                                    make figure H                              ######
###########################################################################################
disease = 'PASS_ADHD_Demontis2018'
risk_score = risk.score[[disease]][, 'zscore']
fdr = risk.score[[disease]][, 'fdr']
meta$risk_score = risk_score
meta$fdr = fdr

# for schizophrenia at adult, look at the box plot for these cell types:
pvalb_cell = rownames(meta)[meta$adjusted_L3 == 'Inh_PVALB']
lamp5_cell = rownames(meta)[meta$adjusted_L3 == 'Inh_LAMP5']
l4_rorb = rownames(meta)[meta$adjusted_L3 == 'Exc_L4-RORB']

# also look at the DRD1-eccentric-CASZ1
drd1_eccentric_cell = rownames(meta)[meta$adjusted_L3 == 'Inh_DRD1-eccentric-CASZ1']

plot_df = meta[c(pvalb_cell, lamp5_cell, l4_rorb, drd1_eccentric_cell), ]
plot_df$cell_type = factor(plot_df$adjusted_L3, levels = c('Inh_DRD1-eccentric-CASZ1', 'Exc_L4-RORB', 'Inh_PVALB', 'Inh_LAMP5'))
plot_df$age_group = factor(plot_df$fine2_age_groups, levels = c('2T', '3T', '1m', '4-7m', 'adult'))

gplot <- ggplot(plot_df, aes(x = age_group, y = risk_score, fill = cell_type)) +
    geom_boxplot() +
    scale_fill_manual(
        name = "L3",  # legend title
        values = c(
        'Exc_L4-RORB' = '#fdc086',
        'Inh_LAMP5' = '#fc8d62',
        'Inh_PVALB' = '#8da0cb',
        'Inh_DRD1-eccentric-CASZ1' = '#66c2a5'
        ),
        labels = c(
        "Exc_L4-RORB" = "L4-RORB",
        "Inh_LAMP5"  = "LAMP5",
        "Inh_PVALB"   = "PVALB",
        "Inh_DRD1-eccentric-CASZ1" = "DRD1\neccentric"
        )
    ) +
    xlab('Cell type') +
    ylab('met-scDRS') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size = 20))

output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-ADHD-L3-met-scdrs-by-age-box-plot.png')
png(
    filename = output.path,
    width = 7,
    height = 5,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

