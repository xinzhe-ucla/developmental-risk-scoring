### visualize the ldsc result as a heatmap of age by cell type

### PREAMBLE ######################################################################################
require(ComplexHeatmap)
require(circlize);

# parse options:
# gather the different regions
require(docopt)
require(tidyverse);

'Usage:
    visualization.R [--ldsc_dir <ldsc_dir> --output <output>]

Options:
    --ldsc_dir directory to heritability estimates .results file
    --output path to draw the heatmap

]' -> doc

# collect user input: 
opts <- docopt(doc)
ldsc_dir <- opts$ldsc_dir;
output.path <- opts$output

# read in the ldsc
ldsc_results = list.files(ldsc_dir, pattern='.results')

# for each of the results, read in the
all_dev_time = all_cell_types = {}
for (h_est in ldsc_results){
    #get meta data:
    dev_time = gsub('_.*', '', h_est)
    cell_type = gsub(paste0(dev_time, '_'), '', h_est)
    cell_type = gsub('.hypo_dmr_overlap.hg19.dmr.results', '', cell_type)
    
    # aggregate to get cell types:
    all_cell_types = c(all_cell_types, cell_type)
    all_dev_time = c(all_dev_time, dev_time)
}

# create a place holder matrix:
all_cell_types = unique(all_cell_types)
all_dev_time = unique(all_dev_time)
h2_df = data.frame(matrix(NA, nrow = length(all_dev_time), ncol = length(all_cell_types)))
colnames(h2_df) = all_cell_types
rownames(h2_df) = all_dev_time

# load in the matrix to fill out the h_est
for (h_est in ldsc_results){
    loaded = read.table(paste0(ldsc_dir, h_est), sep = '\t', header = TRUE)
    
    #get meta data:
    dev_time = gsub('_.*', '', h_est)
    cell_type = gsub(paste0(dev_time, '_'), '', h_est)
    cell_type = gsub('.hypo_dmr_overlap.hg19.dmr.results', '', cell_type)
    
    # record:
    h2_df[dev_time, cell_type] = loaded[1, 'Enrichment']
}

###########################################################################################
######                                    Make Heatmap                               ######
###########################################################################################
# make the column of the heatmap:
excitatory = sort(all_cell_types[grep('^Exc', all_cell_types)])
inhibitory = sort(all_cell_types[grep('^Inh', all_cell_types)])
others = setdiff(all_cell_types, c(excitatory, inhibitory))

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

max_color_scale = ceiling(max(h2_df, na.rm = T))
if (max_color_scale > 100){
    col.fun <- colorRamp2(
        c(
            0,
            50,
            ceiling(max(h2_df, na.rm = T))
            ),
        c('white', '#de2d26', '#7f0000')
        )
    heatmap.legend.param <- list(
        at = c(
            0,
            50,
            ceiling(max(h2_df, na.rm = T))
            )
        );
} else{
    col.fun <- colorRamp2(
        c(
            0,
            ceiling(max(h2_df, na.rm = T))
            ),
        c('white', '#de2d26')
        );

    heatmap.legend.param <- list(
        at = c(
            0,
            ceiling(max(h2_df, na.rm = T))
            )
        );
}


# create heatmap:
plot <- Heatmap(
    as.matrix(h2_df[c('2T', '3T', '1m', '4-7m', 'adult'), cell.type.order]),
    name = 'Enrichment',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = c('2T', '3T', '1m', '4-7m', 'adult'),
    # cluster_rows = TRUE,
    column_order = cell.type.order,
    # cluster_columns = TRUE,
    width = unit(10 * ncol(h2_df),"mm"),
    height = unit(10 * nrow(h2_df),"mm"),
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
pdf(
    file = output.path,
    width = heatmap.width,
    height = heatmap.height
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 10, 70), "mm"));
dev.off();

