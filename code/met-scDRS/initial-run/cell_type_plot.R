### cell_type_plot.R ##############################################################################
# purpose: plot the cell type umap plot:

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(circlize);
library(Seurat);

# define paths:
date <- Sys.Date()
meta.data.path <- "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv"

# load in the data:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# label the '' as 'unknown'
meta$L1[meta$L1 == ''] = 'unknown'
meta$L2[meta$L2 == ''] = 'unknown'
meta$L3[meta$L3 == ''] = 'unknown'

### VISUALIZE - L1 #######################################################################
# create plot df:
plot.df <- data.frame(
    UMAP_1 = meta$UMAP_1,
    UMAP_2 = meta$UMAP_2,
    cell_type = as.factor(meta$L1)
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
    cell_type = as.factor(meta$L2)
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
    cell_type = as.factor(meta$L3)
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
