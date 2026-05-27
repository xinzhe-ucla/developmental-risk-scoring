### proportion-brain-heatmap.R ####################################################################
# purpose: given the proportional matrix, draw out the brain only heatmaps:

### PREAMBLE ######################################################################################
# ### brain-trait-heatmap.R #########################################################################
# purpose: draw a heatmap that looks at the proportion of significant cells per cell type

### PREAMBLE ######################################################################################
# load packages:
require(ggplot2);
require(ComplexHeatmap);
require(circlize);

# load in data:
significance.matrix <- read.table(
    file = '/u/home/l/lixinzhe/project-geschwind/plot/2025-10-08-revised-BICAN-mcg-l3-proportion.csv',
    sep = ',',
    row.names = 1,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
    );

# load in the trait info:
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

system.date <- Sys.Date()

# label the '' as 'unknown'
meta.data.path = '/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz'
meta = read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    )
meta$newL1[meta$newL1 == ''] = 'unknown'
meta$newL1[meta$newL1 == '?'] = 'unknown'
meta$newL2[meta$newL2 == ''] = 'unknown'
meta$newL2[meta$newL2 == '?'] = 'unknown'
meta$newL3[meta$newL3 == ''] = 'unknown'
meta$newL3[meta$newL3 == '?'] = 'unknown'

# for L2 and L3, add the word Exc and Inh on top:
meta$newL2 = paste0(meta$newL1, '_', meta$newL2)
meta$newL3 = paste0(meta$newL1, '_', meta$newL3)

### VISUALIZATION #################################################################################
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
publication.traits <- rownames(significance.matrix)[trait.class == 'brain'];

# grab out the cell types:
cell.types = unique(meta[, 'newL3'])
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
    cluster_rows = FALSE,
    #column_order = cell.type.order,
    cluster_columns = FALSE,
    width = unit(10 * length(cell.type.order),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    # row_split = row.split,
    column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-developmental-cell-type-brain-traits-proportion.png')
png(
    filename = output.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();
