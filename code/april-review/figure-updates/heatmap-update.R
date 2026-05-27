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
imaging_trait = list.names[grep('UKB_IDP', list.names)]
trait.info.expand = data.frame(
    c(trait.info$Trait_Identifier, imaging_trait),
    c(trait.info$Category, rep('brain', length(imaging_trait)))
)
colnames(trait.info.expand) = colnames(trait.info)
trait.info = trait.info.expand

# get the brain trait:
brain_traits = trait.info$Trait_Identifier[trait.info$Category == 'brain']

### look at the adult only data:
adult_cells <- rownames(meta)[meta$fine2_age_groups == 'adult']

# get the proportion of cells that are significant in each of the disease x cell type
meta_adult = meta[adult_cells, ]

significance_matrix = data.frame(matrix(NA, nrow = length(risk.score), ncol = length(unique(meta_adult$adjusted_L3))))
colnames(significance_matrix) = unique(meta_adult$adjusted_L3)
rownames(significance_matrix) = names(risk.score)

for (disease in names(risk.score)){
    adult_only_risk_score = risk.score[[disease]][adult_cells, 'zscore']
    adult_only_fdr = risk.score[[disease]][adult_cells, 'fdr']
    
    # bind the meta to the risk score
    meta_adult$risk_score = adult_only_risk_score
    meta_adult$fdr = adult_only_fdr
    
    # calculate proportion:
    proprotion_sig = meta_adult %>% group_by(adjusted_L3) %>% 
        summarize(
            n = n(),
            n_sig = sum(fdr < 0.1, na.rm = TRUE),
            proportion = n_sig / n,
            .groups = "drop"
        )
    
    # coerce into named vector:
    sig_by_celltype_in_disease = c(proprotion_sig$proportion)
    names(sig_by_celltype_in_disease) = proprotion_sig$adjusted_L3
    sig_by_celltype_in_disease = sig_by_celltype_in_disease[colnames(significance_matrix)]
    
    # put the result into a significant proportion matrix:
    significance_matrix[disease, ] = sig_by_celltype_in_disease
}
adult_only = significance_matrix

### look at all the cells:
stopifnot(rownames(meta) == rownames(risk.score[[1]]))

significance_matrix = data.frame(matrix(NA, nrow = length(risk.score), ncol = length(unique(meta$adjusted_L3))))
colnames(significance_matrix) = unique(meta$adjusted_L3)
rownames(significance_matrix) = names(risk.score)

for (disease in names(risk.score)){
    met_z_score = risk.score[[disease]][, 'zscore']
    met_fdr = risk.score[[disease]][, 'fdr']
    
    # bind the meta to the risk score
    meta$risk_score = met_z_score
    meta$fdr = met_fdr
    
    # calculate proportion:
    proprotion_sig = meta %>% group_by(adjusted_L3) %>% 
        summarize(
            n = n(),
            n_sig = sum(fdr < 0.1, na.rm = TRUE),
            proportion = n_sig / n,
            .groups = "drop"
        )
    
    # coerce into named vector:
    sig_by_celltype_in_disease = c(proprotion_sig$proportion)
    names(sig_by_celltype_in_disease) = proprotion_sig$adjusted_L3
    sig_by_celltype_in_disease = sig_by_celltype_in_disease[colnames(significance_matrix)]
    
    # put the result into a significant proportion matrix:
    significance_matrix[disease, ] = sig_by_celltype_in_disease
}
all_sig_mat = significance_matrix



### proposed figure ###############################################################################
# figure A: figure of heatmap on proportion of significant cells in traits vs cell types:

### make the matrix ###
significance.matrix = adult_only

# find out the set of brain traits:
trait.info$Trait_Identifier <- gsub('PASS_', '', trait.info$Trait_Identifier)
trait.info$Trait_Identifier <- gsub('UKB_460K.', '', trait.info$Trait_Identifier)
trait.info$Trait_Identifier <- gsub('cov_', '', trait.info$Trait_Identifier)
trait.info$Trait_Identifier <- gsub('repro_', '', trait.info$Trait_Identifier)

rownames(significance.matrix) <- gsub('PASS_', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('UKB_460K.', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('cov_', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('repro_', '', rownames(significance.matrix))

trait.class <- trait.info$Category[match(rownames(significance.matrix), trait.info$Trait_Identifier)];

# select traits to plot:
publication.traits = c(
    'ADHD_Demontis2018',
    'BIP_Mullins2021',
    'MDD_Howard2019',
    'UKB_IDP0013',
    'UKB_IDP0127',
    'Schizophrenia_Pardinas2018',
    'EDU_YEARS',
    'Type_1_Diabetes',
    'biochemistry_Cholesterol',
    'Type_2_Diabetes',
    'blood_EOSINOPHIL_COUNT',
    'cancer_BREAST'
    )
#publication.traits <- rownames(significance.matrix)[trait.class == 'brain'];

# grab out the cell types:
cell.types = unique(meta_adult[, 'adjusted_L3'])
excitatory = sort(cell.types[grep('^Exc', cell.types)])
inhibitory = sort(cell.types[grep('^Inh', cell.types)])
others = setdiff(cell.types, c(excitatory, inhibitory))

cell.type.order = c(
    excitatory,
    inhibitory,
    others
    )

column.split = c(
    rep('Excitatory', length(excitatory)),
    rep('Inhibitory', length(inhibitory)),
    rep('Others', length(others))
    )

# make color function
col.fun <- colorRamp2(
    c(
        0,
        1
        ),
    c('white', '#de2d26')
    );
heatmap.legend.param <- list(
    at = c(
        0,
        1
        )
    );

# create heatmap:
plot <- Heatmap(
    as.matrix(significance.matrix)[publication.traits, cell.type.order],
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    #row_order = publication.traits,
    cluster_rows = TRUE,
    #column_order = cell.type.order,
    cluster_columns = TRUE,
    width = unit(10 * length(cell.type.order),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    # row_split = row.split,
    column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-cell-type-selected-traits-proportion-adult-only.pdf')
pdf(
    file = output.path,
    width = heatmap.width,
    height = heatmap.height
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 10, 70), "mm"));
dev.off();

### make the matrix ###
write.table(
    significance.matrix,
    file = paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-cell-type-selected-traits-proportion-adult-only.csv'),
    sep = ',',
    col.names = TRUE,
    row.names = TRUE,
    quote = FALSE
    )

###########################################################################################
######                            all cells across time                              ######
###########################################################################################
significance.matrix = all_sig_mat

# find out the set of brain traits:
trait.info$Trait_Identifier <- gsub('PASS_', '', trait.info$Trait_Identifier)
trait.info$Trait_Identifier <- gsub('UKB_460K.', '', trait.info$Trait_Identifier)
trait.info$Trait_Identifier <- gsub('cov_', '', trait.info$Trait_Identifier)
trait.info$Trait_Identifier <- gsub('repro_', '', trait.info$Trait_Identifier)

rownames(significance.matrix) <- gsub('PASS_', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('UKB_460K.', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('cov_', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('repro_', '', rownames(significance.matrix))

trait.class <- trait.info$Category[match(rownames(significance.matrix), trait.info$Trait_Identifier)];

# select traits to plot:
#publication.traits <- rownames(significance.matrix)[trait.class == 'brain'];

# grab out the cell types:
cell.types = unique(meta_adult[, 'adjusted_L3'])
excitatory = sort(cell.types[grep('^Exc', cell.types)])
inhibitory = sort(cell.types[grep('^Inh', cell.types)])
others = setdiff(cell.types, c(excitatory, inhibitory))

cell.type.order = c(
    excitatory,
    inhibitory,
    others
    )

column.split = c(
    rep('Excitatory', length(excitatory)),
    rep('Inhibitory', length(inhibitory)),
    rep('Others', length(others))
    )

# make color function
col.fun <- colorRamp2(
    c(
        0,
        1
        ),
    c('white', '#de2d26')
    );
heatmap.legend.param <- list(
    at = c(
        0,
        1
        )
    );

# create heatmap:
plot <- Heatmap(
    as.matrix(significance.matrix)[publication.traits, cell.type.order],
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    #row_order = publication.traits,
    cluster_rows = TRUE,
    #column_order = cell.type.order,
    cluster_columns = TRUE,
    width = unit(10 * length(cell.type.order),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    # row_split = row.split,
    column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-cell-type-selected-traits-proportion-all-timepoint.pdf')
pdf(
    file = output.path,
    width = heatmap.width,
    height = heatmap.height
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 10, 70), "mm"));
dev.off();

### make the matrix ###
write.table(
    significance.matrix,
    file = paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-cell-type-selected-traits-proportion-all-timepoint.csv'),
    sep = ',',
    col.names = TRUE,
    row.names = TRUE,
    quote = FALSE
    )




###########################################################################################
######                            just the volume trait                              ######
###########################################################################################
significance.matrix = all_sig_mat

# find out the set of brain traits:
trait.info$Trait_Identifier <- gsub('PASS_', '', trait.info$Trait_Identifier)
trait.info$Trait_Identifier <- gsub('UKB_460K.', '', trait.info$Trait_Identifier)
trait.info$Trait_Identifier <- gsub('cov_', '', trait.info$Trait_Identifier)
trait.info$Trait_Identifier <- gsub('repro_', '', trait.info$Trait_Identifier)

rownames(significance.matrix) <- gsub('PASS_', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('UKB_460K.', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('cov_', '', rownames(significance.matrix))
rownames(significance.matrix) <- gsub('repro_', '', rownames(significance.matrix))

trait.class <- trait.info$Category[match(rownames(significance.matrix), trait.info$Trait_Identifier)];

# select traits to plot:
#publication.traits <- rownames(significance.matrix)[trait.class == 'brain'];

# grab out the cell types:
cell.types = unique(meta_adult[, 'adjusted_L3'])
excitatory = sort(cell.types[grep('^Exc', cell.types)])
inhibitory = sort(cell.types[grep('^Inh', cell.types)])
others = setdiff(cell.types, c(excitatory, inhibitory))

publication.traits = c(
    "UKB_IDP0013",
    "UKB_IDP0014",
    "UKB_IDP0015",
    "UKB_IDP0016",
    "UKB_IDP0124",
    "UKB_IDP0125",
    "UKB_IDP0126",
    "UKB_IDP0127"
    )

cell.type.order = c(
    excitatory,
    inhibitory,
    others
    )

column.split = c(
    rep('Excitatory', length(excitatory)),
    rep('Inhibitory', length(inhibitory)),
    rep('Others', length(others))
    )

# make color function
col.fun <- colorRamp2(
    c(
        0,
        0.1
        ),
    c('white', '#de2d26')
    );
heatmap.legend.param <- list(
    at = c(
        0,
        0.1
        )
    );

# create heatmap:
plot <- Heatmap(
    as.matrix(significance.matrix)[publication.traits, cell.type.order],
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    #row_order = publication.traits,
    cluster_rows = TRUE,
    #column_order = cell.type.order,
    cluster_columns = TRUE,
    width = unit(10 * length(cell.type.order),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    # row_split = row.split,
    column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-cell-type-heatmap-proportion-all-timepoint-brain-volume.pdf')
pdf(
    file = output.path,
    width = heatmap.width,
    height = heatmap.height
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 10, 70), "mm"));
dev.off();







###########################################################################################
######                              Look at the box plot for IDP                     ######
###########################################################################################
for (disease in imaging_trait){
    risk_score = risk.score[[disease]][, 'zscore']
    fdr = risk.score[[disease]][, 'fdr']
    meta$risk_score = risk_score
    meta$fdr = fdr

    # get these cell types:
    opc = rownames(meta)[meta$adjusted_L3 == 'Glial_OPC']
    odc = rownames(meta)[meta$adjusted_L3 == 'Glial_ODC']
    astro = rownames(meta)[meta$adjusted_L3 == 'Glial_Astro']
    exc = rownames(meta)[meta$newL1 == 'Exc']
    drd = rownames(meta)[grep('Inh_DRD', meta$adjusted_L3)]
    
    plot_df = meta[c(astro, odc, opc, exc, drd), ]
    plot_df$cell_type = plot_df$adjusted_L3
    plot_df[drd, 'cell_type'] = 'Inh_DRD'
    plot_df[exc, 'cell_type'] = plot_df[exc, 'newL1']
    plot_df$cell_type = factor(plot_df$cell_type, levels = c('Glial_OPC', 'Glial_ODC', 'Glial_Astro',  'Inh_DRD', 'Exc'))
    plot_df$age_group = factor(plot_df$fine2_age_groups, levels = c('2T', '3T', '1m', '4-7m', 'adult'))

    gplot <- ggplot(plot_df, aes(x = age_group, y = risk_score, fill = cell_type)) +
        geom_boxplot() +
        scale_fill_manual(
            name = "cell type",  # legend title
            values = c(
            'Exc' = '#fc8d62', # darker orange
            'Inh_DRD' = '#8da0cb', #blue
            'Glial_OPC' = '#b2e2e2', # green
            'Glial_ODC' = '#66c2a4', # greener
            'Glial_Astro' = '#2ca25f'  # even more green
            ),
            labels = c(
            "Glial_OPC" = "OPC",
            "Glial_ODC"  = "ODC",
            "Glial_Astro"   = "Astro",
            "Exc" = "Exc",
            "Inh_DRD" = "DRD"
            )
        ) +
        xlab('Cell type') +
        ylab('met-scDRS') +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme(text = element_text(size = 20))

    output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-', disease, '-L3-met-scdrs-by-age-box-plot.pdf')
    pdf(
        file = output.path,
        width = 7,
        height = 5
        );
    print(gplot)
    dev.off();
}
