### spatial-age-in-traits #########################################################################
library(data.table)
library(tidyverse)

meta = read.table('/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz', sep = ',', header=TRUE, row.names = 1)
scDRS.directory = '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/'

score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

mdd = data.frame(
    fread(
        file =  "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised//PASS_MDD_Howard2019.score.gz",
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
meta$newL3[meta$newL3 == ''] = 'unknown'
meta$newL3[meta$newL3 == '?'] = 'unknown'

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

# look at the dataframe:
disease = 'PASS_Schizophrenia_Pardinas2018'
frame = cbind(risk.score[[disease]], meta)
quick = frame %>% group_by(newL2, fine_age_groups) %>% summarise('average_z' = mean(zscore)) %>% as.data.frame()    

# for each of the cell type in each of the disease, figure out big variance
result_df = data.frame(matrix(NA, nrow = length(unique(meta$newL3)), ncol = length(risk.score)))
rownames(result_df) = unique(meta$newL3)
colnames(result_df) = names(risk.score)

for (disease in names(risk.score)){
    frame = cbind(risk.score[[disease]], meta)
    quick = frame %>% group_by(newL3, fine_age_groups) %>% summarise('average_z' = mean(zscore)) %>% as.data.frame()    
    
    for (cell_type in unique(meta$newL3)){
        # get the variance:
        result_df[cell_type, disease] = var(quick[quick$newL3 == cell_type, 'average_z'])
    }
}

# from the matrix we can see L1-3-MLLT1 have high variance for scz:
 # 'L4-RORB' could work for PASS_Schizophrenia_Pardinas2018
 # result_df[, 'PASS_Schizophrenia_Pardinas2018', drop = F]

for (cell_type in unique(meta$newL3)){
    cells_of_interest = rownames(meta)[meta$newL3 == cell_type]
    plot_df = meta[cells_of_interest, c('newL3', 'fine_age_groups')]
    plot_df$met_scdrs = risk.score[['PASS_Schizophrenia_Pardinas2018']][cells_of_interest ,'zscore']
    print(plot_df %>% group_by(fine_age_groups) %>% summarize(mean(met_scdrs)))

}

# plot it out:
cells_of_interest = rownames(meta)[meta$newL3 == 'L4-PLCH1']
plot_df = meta[cells_of_interest, c('newL3', 'fine_age_groups')]
plot_df$met_scdrs = risk.score[['PASS_Schizophrenia_Pardinas2018']][cells_of_interest ,'zscore']
print(plot_df %>% group_by(fine_age_groups) %>% summarize(mean(met_scdrs)))

gplot <- ggplot(plot_df, aes(x = fine_age_groups, y = met_scdrs, fill = fine_age_groups)) +
    geom_boxplot() +
    scale_fill_manual(
        values = c(
            '1m' = '#fdc086',
            '4m' = '#fc8d62',
            '7m' = '#8da0cb',
            'adult' = '#66c2a5'
            # 'middle temporal gyrus' = '#e78ac3'
            )
        ) +
    xlab('Age groups') +
    ylab('met-scDRS') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size = 20))

output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-scz-fine-age-box-plot.png')
png(
    filename = output.path,
    width = 8,
    height = 7,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();


# cell type by fine age plot for SCZ
disease = 'PASS_Schizophrenia_Pardinas2018'
frame = cbind(risk.score[[disease]], meta)
sig_summary <- frame %>%
  group_by(newL3, fine_age_groups) %>%
  summarise(
    n_total = n(),
    n_sig   = sum(fdr < 0.1, na.rm = TRUE),
    prop_sig = n_sig / n_total,
    .groups = "drop"
  )

# wide matrix: rows = newL3, cols = fine_age_groups, values = prop_sig
mat_df <- sig_summary %>%
  select(newL3, fine_age_groups, prop_sig) %>%
  pivot_wider(names_from = fine_age_groups, values_from = prop_sig)

# to matrix
rn <- mat_df$newL3
mat <- as.matrix(mat_df[,-1])
rownames(mat) <- rn

library(ComplexHeatmap)
library(circlize)

# color scale 0 → 1
col_fun <- colorRamp2(c(0, 0.5, 1), c("#F7FBFF", "#6BAED6", "#08306B"))

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

plot <- Heatmap(
    mat,
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    #row_order = publication.traits,
    cluster_rows = FALSE,
    #column_order = cell.type.order,
    cluster_columns = FALSE,
    width = unit(10 * ncol(mat),"mm"),
    height = unit(10 * nrow(mat),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    # row_split = row.split,
    #column_split = column.split,
    heatmap_legend_param = heatmap.legend.param
    );

output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-scz-fine-age-by-celltype-plot.png')
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
png(
    filename = output.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );

draw(plot)
dev.off();
