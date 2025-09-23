### cell_type_plot.R ##############################################################################
# purpose: plot the cell type umap plot:

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(circlize);
library(Seurat);

# define paths:
date <- Sys.Date()
meta.data.path <- "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz"
umap_path = "/u/project/cluo/heffel/BICAN3/3C/bican_2025_3Cumap.csv.gz"

umap = read.csv(
    file = umap_path,
    row.names = 1,
    header = TRUE
    )

# load in the data:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

meta$newL1[meta$newL1 == ''] = 'unknown'
meta$newL1[meta$newL1 == '?'] = 'unknown'
meta$newL2[meta$newL2 == ''] = 'unknown'
meta$newL2[meta$newL2 == '?'] = 'unknown'
meta$newL3[meta$newL3 == ''] = 'unknown'
meta$newL3[meta$newL3 == '?'] = 'unknown'

# for L2 and L3, add the word Exc and Inh on top:
meta$newL2 = paste0(meta$newL1, '_', meta$newL2)
meta$newL3 = paste0(meta$newL1, '_', meta$newL3)

meta$UMAP_1 = umap[rownames(meta), 'X_umap']
meta$UMAP_2 = umap[rownames(meta), 'Y_umap']

### VISUALIZE - L1 #######################################################################
# create plot df:
plot.df <- data.frame(
    UMAP_1 = meta$UMAP_1,
    UMAP_2 = meta$UMAP_2,
    cell_type = as.factor(meta$newL1)
    );
rownames(plot.df) = rownames(meta);

# Create plot:
gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = cell_type, alpha = 0.5)) +
    geom_point() +
    theme_classic() +
    ggtitle('Developmental UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 20)) +
    theme(legend.position="none")
gplot.label <- LabelClusters(plot = gplot, id = 'cell_type', col = 'black', size = 5)

plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/", date, '-QCed-developmental-cell-type-L1-umap.png')
png(
    filename = plot.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot.label)
dev.off();

### L2 ############################################################################################
# create plot df:
plot.df <- data.frame(
    UMAP_1 = meta$UMAP_1,
    UMAP_2 = meta$UMAP_2,
    cell_type = as.factor(meta$newL2)
    );
rownames(plot.df) = rownames(meta);

# Create plot:
gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = cell_type, alpha = 0.5)) +
    geom_point() +
    theme_classic() +
    ggtitle('Developmental UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 20)) +
    theme(legend.position="none")
gplot.label <- LabelClusters(plot = gplot, id = 'cell_type', col = 'black', size = 5)

plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/", date, '-QCed-developmental-cell-type-L2-umap.png')
png(
    filename = plot.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot.label)
dev.off();

### L3 ############################################################################################
# create plot df:
plot.df <- data.frame(
    UMAP_1 = meta$UMAP_1,
    UMAP_2 = meta$UMAP_2,
    cell_type = as.factor(meta$newL3)
    );
rownames(plot.df) = rownames(meta);

# Create plot:
gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = cell_type, alpha = 0.5)) +
    geom_point() +
    theme_classic() +
    ggtitle('Developmental UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 20)) +
    theme(legend.position="none")
gplot.label <- LabelClusters(
    plot = gplot,
    id = 'cell_type',
    col = 'black',
    size = 5,
    repel = TRUE,
    max.overlaps = Inf
    )

plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/", date, '-QCed-developmental-cell-type-L3-umap.png')
png(
    filename = plot.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot.label)
dev.off();

### FINE AGE ######################################################################################
plot.df <- data.frame(
    UMAP_1 = meta$UMAP_1,
    UMAP_2 = meta$UMAP_2,
    cell_type = as.factor(meta$fine_age_groups)
    );
rownames(plot.df) = rownames(meta);

# Create plot:
gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = cell_type, alpha = 0.5)) +
    geom_point() +
    theme_classic() +
    ggtitle('Developmental UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 20))

plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/", date, '-QCed-developmental-fine-age-umap.png')
png(
    filename = plot.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();