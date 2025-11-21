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
######                                    FIGURE B                                   ######
###########################################################################################
## bind the umap into meta as well:
meta$score = risk.score[['PASS_Schizophrenia_Pardinas2018']]$zscore
meta$fdr = risk.score[['PASS_Schizophrenia_Pardinas2018']]$fdr

# draw umap:
p.cutoff = 0.1
significant.cell <- rownames(meta)[meta$fdr < p.cutoff];
insignificant.cell <- setdiff(rownames(meta), significant.cell);

# next, we will plot out the umap:
plot.df <- meta;

gplot <- ggplot(plot.df, aes(x = X_umap, y = Y_umap)) +
    geom_point(data = plot.df[insignificant.cell, ], colour = 'grey', size = 1) +
    geom_point(data = plot.df[significant.cell, ], aes(colour = score), size = 1) +
    scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
    theme_classic() +
    ggtitle('SCZ Met-scDRS') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 35))+
    xlim(-15, 15) +
    ylim(-15, 18)

# draw out the plot:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-all-SCZ-scDRS-score-umap.png')
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

# draw umap:
p.cutoff = 0.1
significant.cell <- rownames(meta)[meta$fdr < p.cutoff];
insignificant.cell <- setdiff(rownames(meta), significant.cell);

###########################################################################################
######                                    FIGURE C                              ######
###########################################################################################
# next, we plot out umap based only on 1 m 
month_1 = meta$fine2_age_groups == '1m'
plot.df <- meta[month_1, ];
significant.cell <- rownames(plot.df)[plot.df$fdr < p.cutoff];
insignificant.cell <- setdiff(rownames(plot.df), significant.cell);


gplot <- ggplot(plot.df, aes(x = X_umap, y = Y_umap)) +
    geom_point(data = plot.df[insignificant.cell, ], colour = 'grey', size = 1) +
    geom_point(data = plot.df[significant.cell, ], aes(colour = score), size = 1) +
    scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
    theme_classic() +
    ggtitle('1-month SCZ Met-scDRS') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 35))+
    xlim(-15, 15) +
    ylim(-15, 18)

# draw out the plot:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-1m-SCZ-scDRS-score-umap.png')
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

###########################################################################################
######                                    FIGURE D                                   ######
###########################################################################################
## bind the umap into meta as well:
meta$score = risk.score[['PASS_ADHD_Demontis2018']]$zscore
meta$fdr = risk.score[['PASS_ADHD_Demontis2018']]$fdr

# draw umap:
p.cutoff = 0.1
significant.cell <- rownames(meta)[meta$fdr < p.cutoff];
insignificant.cell <- setdiff(rownames(meta), significant.cell);

# next, we will plot out the umap:
plot.df <- meta;

gplot <- ggplot(plot.df, aes(x = X_umap, y = Y_umap)) +
    geom_point(data = plot.df[insignificant.cell, ], colour = 'grey', size = 1) +
    geom_point(data = plot.df[significant.cell, ], aes(colour = score), size = 1) +
    scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
    theme_classic() +
    ggtitle('ADHD Met-scDRS') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 35))+
    xlim(-15, 15) +
    ylim(-15, 18)

# draw out the plot:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-all-ADHD-scDRS-score-umap.png')
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

# draw umap:
p.cutoff = 0.1
significant.cell <- rownames(meta)[meta$fdr < p.cutoff];
insignificant.cell <- setdiff(rownames(meta), significant.cell);

###########################################################################################
######                                    FIGURE E                              ######
###########################################################################################
# next, we plot out umap based only on 1 m 
month_1 = meta$fine2_age_groups == '1m'
plot.df <- meta[month_1, ];
significant.cell <- rownames(plot.df)[plot.df$fdr < p.cutoff];
insignificant.cell <- setdiff(rownames(plot.df), significant.cell);


gplot <- ggplot(plot.df, aes(x = X_umap, y = Y_umap)) +
    geom_point(data = plot.df[insignificant.cell, ], colour = 'grey', size = 1) +
    geom_point(data = plot.df[significant.cell, ], aes(colour = score), size = 1) +
    scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
    theme_classic() +
    ggtitle('1-month ADHD Met-scDRS') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 35))+
    xlim(-15, 15) +
    ylim(-15, 18)

# draw out the plot:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-1m-ADHD-scDRS-score-umap.png')
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();
