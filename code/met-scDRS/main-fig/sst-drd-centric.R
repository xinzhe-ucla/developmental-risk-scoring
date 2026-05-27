### main figure that is sst and drd neurons centric over diseases 

### PREAMBLE ######################################################################################
library(data.table)
library(tidyverse)
library(clusterProfiler);
library(ComplexHeatmap)
require(circlize);

meta = read.table('/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz', sep = ',', header=TRUE, row.names = 1)
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

trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

brain_traits = trait.info$Trait_Identifier[trait.info$Category == 'brain']

###########################################################################################
######                               format risk scores                              ######
###########################################################################################
# stack risk scores:
subset_meta = meta[rownames(risk.score[[1]]), c('newL1', 'newL2', 'newL3', 'fine_age_groups')]
stacked_df_collection = NULL
for (result in names(risk.score)){
    stacked_df = cbind(risk.score[[result]][, c('zscore', 'fdr')], subset_meta)
    stacked_df$disease = result
    stacked_df_collection[[result]] = stacked_df
}
stacked = Reduce('rbind', stacked_df_collection)

# get the SST and the DRD cells:
brain_stacked = stacked[stacked$disease %in% brain_traits, ]

### SST NXPH2 ###

round_up_to_first_signif <- function(x) {
  s <- floor(log10(abs(x)))   # position of first significant digit
  out <- ceiling(x / 10^s) * 10^s
  out
}

plot_traits_and_age_sig <- function(stacked_frame, cell_type = 'SST-NXPH2', output.path){
    cell_type_stacked = stacked_frame[stacked_frame$newL3 == cell_type,]
    df_sum <- cell_type_stacked %>%
        group_by(fine_age_groups, disease) %>%
        summarize(
            n = n(),
            n_sig = sum(fdr < 0.1, na.rm = TRUE),
            proportion = n_sig / n,
            .groups = "drop"
        )
    significance.matrix <- df_sum %>%
        select(disease, fine_age_groups, proportion) %>%
        tidyr::pivot_wider(names_from = fine_age_groups, values_from = proportion) %>%
        tibble::column_to_rownames("disease") %>%
        as.matrix()
    
    # make color function
    col.fun <- colorRamp2(
        c(
            0,
            round_up_to_first_signif(max(significance.matrix))
            ),
        c('white', '#de2d26')
        );
    heatmap.legend.param <- list(
        at = c(
            0,
            round_up_to_first_signif(max(significance.matrix))
            )
        );
    
    # create heatmap:
    plot <- Heatmap(
        as.matrix(significance.matrix),
        name = 'Sig. cells',
        col = col.fun,
        rect_gp = gpar(col = "black", lwd = 2),
        #row_order = publication.traits,
        cluster_rows = FALSE,
        #column_order = cell.type.order,
        cluster_columns = FALSE,
        width = unit(10 * ncol(significance.matrix),"mm"),
        height = unit(10 * nrow(significance.matrix),"mm"),
        column_names_gp = grid::gpar(fontsize = 15),
        row_names_gp = grid::gpar(fontsize = 15),
        # row_split = row.split,
        # column_split = column.split,
        heatmap_legend_param = heatmap.legend.param
        );
    plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
    # measure the size of the heatmap:
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
    draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
    dev.off();
}

plot_traits_and_age_sig(
    brain_stacked,
    cell_type = 'SST-NXPH2',
    paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-SST_NXPH2-time-point-brain-traits-proportion.png')
)

plot_traits_and_age_sig(
    brain_stacked,
    cell_type = 'SST-THRDE',
    paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-SST_THRDE-time-point-brain-traits-proportion.png')
)
  
plot_traits_and_age_sig(
    brain_stacked,
    cell_type = 'DRD1',
    paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-DRD1-time-point-brain-traits-proportion.png')
)

plot_traits_and_age_sig(
    brain_stacked,
    cell_type = 'DRD2',
    paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-DRD2-time-point-brain-traits-proportion.png')
)

###########################################################################################
######                           met-scDRS with age on scz                           ######
###########################################################################################
# DRD
lineage_step1 = rownames(meta)[meta$fine_age_groups == '2T' & meta$newL3 == 'DRD1']
lineage_step2 = rownames(meta)[meta$fine_age_groups == '3T' & meta$newL3 == 'DRD1']
lineage_step3 = rownames(meta)[meta$fine_age_groups == '1m' & meta$newL3 == 'DRD1']
lineage_step4 = rownames(meta)[meta$fine_age_groups == '4m' & meta$newL3 == 'DRD1']
lineage_step5 = rownames(meta)[meta$fine_age_groups == '7m' & meta$newL3 == 'DRD1']
lineage_step6 = rownames(meta)[meta$fine_age_groups == 'adult' & meta$newL3 == 'DRD1']

# for each step:
disease = 'PASS_Schizophrenia_Pardinas2018'
score = risk.score[[disease]]
total_lineage = list(lineage_step1, lineage_step2, lineage_step3, lineage_step4, lineage_step5, lineage_step6)
names(total_lineage) = paste0('step', seq(1, length(total_lineage)))

plot_dfs = NULL
for (step in names(total_lineage)){
    if (length(total_lineage[[step]]) > 0){
        plot_dfs[[step]] = data.frame(
            'score' = score[total_lineage[[step]], ],
            'lineage_step' = step
            )
    }
}
plot_df = Reduce('rbind', plot_dfs)

gplot <- ggplot(plot_df, aes(x = lineage_step, y = score.zscore, fill = lineage_step)) +
    geom_boxplot() +
    # scale_fill_manual(
    #     values = c(
            
    #         '1m' = '#fdc086',
    #         '4m' = '#fc8d62',
    #         '7m' = '#8da0cb',
    #         'adult' = '#66c2a5'
    #         # 'middle temporal gyrus' = '#e78ac3'
    #         )
    #     ) +
    scale_x_discrete(labels = c(
        'step1' = '2T DRD1',
        'step2' = '3T DRD1',
        'step3' = '1M DRD1',
        'step4' = '4M DRD1',
        'step5' = '7M DRD1',
        'step6' = 'Adult DRD1'
        ))+
    xlab('Lineage Step') +
    ylab('met-scDRS') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size = 20)) +
    theme(legend.position = "none")

output.path <- paste0('/u/home/l/lixinzhe/project-cluo/plot/', Sys.Date(), '-SCZ-DRD1-lineage-step-box-plot.png')
png(
    filename = output.path,
    width = 15,
    height = 10,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

###########################################################################################
######                               met-scDRS with age                              ######
###########################################################################################
# SST THRDE
lineage_step1 = rownames(meta)[meta$fine_age_groups == '2T' & meta$newL2 == 'MGE']
lineage_step2 = rownames(meta)[meta$fine_age_groups == '3T' & meta$newL3 == 'SST-NXPH2']
lineage_step3 = rownames(meta)[meta$fine_age_groups == '1m' & meta$newL3 == 'SST-NXPH2']
lineage_step4 = rownames(meta)[meta$fine_age_groups == '4m' & meta$newL3 == 'SST-NXPH2']
lineage_step5 = rownames(meta)[meta$fine_age_groups == '7m' & meta$newL3 == 'SST-NXPH2']
lineage_step6 = rownames(meta)[meta$fine_age_groups == 'adult' & meta$newL3 == 'SST-NXPH2']

# for each step:
disease = 'PASS_Schizophrenia_Pardinas2018'
score = risk.score[[disease]]
total_lineage = list(lineage_step1, lineage_step2, lineage_step3, lineage_step4, lineage_step5, lineage_step6)
names(total_lineage) = paste0('step', seq(1, length(total_lineage)))

plot_dfs = NULL
for (step in names(total_lineage)){
    if (length(total_lineage[[step]]) > 0){
        plot_dfs[[step]] = data.frame(
            'score' = score[total_lineage[[step]], ],
            'lineage_step' = step
            )
    }
}
plot_df = Reduce('rbind', plot_dfs)

gplot <- ggplot(plot_df, aes(x = lineage_step, y = score.zscore, fill = lineage_step)) +
    geom_boxplot() +
    # scale_fill_manual(
    #     values = c(
            
    #         '1m' = '#fdc086',
    #         '4m' = '#fc8d62',
    #         '7m' = '#8da0cb',
    #         'adult' = '#66c2a5'
    #         # 'middle temporal gyrus' = '#e78ac3'
    #         )
    #     ) +
    scale_x_discrete(labels = c(
        'step1' = '2T MGE',
        'step2' = '3T SST-NXPH2',
        'step3' = '1M SST-NXPH2',
        'step4' = '4M SST-NXPH2',
        'step5' = '7M SST-NXPH2',
        'step6' = 'Adult SST-NXPH2'
        ))+
    xlab('Lineage Step') +
    ylab('met-scDRS') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size = 20)) +
    theme(legend.position = "none")

output.path <- paste0('/u/home/l/lixinzhe/project-cluo/plot/', Sys.Date(), '-SCZ-SST-NXPH2-lineage-step-box-plot.png')
png(
    filename = output.path,
    width = 15,
    height = 10,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

###########################################################################################
######                           met-scDRS with age on BP                           ######
###########################################################################################
# DRD
lineage_step1 = rownames(meta)[meta$fine_age_groups == '2T' & meta$newL3 == 'DRD1']
lineage_step2 = rownames(meta)[meta$fine_age_groups == '3T' & meta$newL3 == 'DRD1']
lineage_step3 = rownames(meta)[meta$fine_age_groups == '1m' & meta$newL3 == 'DRD1']
lineage_step4 = rownames(meta)[meta$fine_age_groups == '4m' & meta$newL3 == 'DRD1']
lineage_step5 = rownames(meta)[meta$fine_age_groups == '7m' & meta$newL3 == 'DRD1']
lineage_step6 = rownames(meta)[meta$fine_age_groups == 'adult' & meta$newL3 == 'DRD1']

# for each step:
disease = 'PASS_BIP_Mullins2021'
score = risk.score[[disease]]
total_lineage = list(lineage_step1, lineage_step2, lineage_step3, lineage_step4, lineage_step5, lineage_step6)
names(total_lineage) = paste0('step', seq(1, length(total_lineage)))

plot_dfs = NULL
for (step in names(total_lineage)){
    if (length(total_lineage[[step]]) > 0){
        plot_dfs[[step]] = data.frame(
            'score' = score[total_lineage[[step]], ],
            'lineage_step' = step
            )
    }
}
plot_df = Reduce('rbind', plot_dfs)

gplot <- ggplot(plot_df, aes(x = lineage_step, y = score.zscore, fill = lineage_step)) +
    geom_boxplot() +
    # scale_fill_manual(
    #     values = c(
            
    #         '1m' = '#fdc086',
    #         '4m' = '#fc8d62',
    #         '7m' = '#8da0cb',
    #         'adult' = '#66c2a5'
    #         # 'middle temporal gyrus' = '#e78ac3'
    #         )
    #     ) +
    scale_x_discrete(labels = c(
        'step1' = '2T DRD1',
        'step2' = '3T DRD1',
        'step3' = '1M DRD1',
        'step4' = '4M DRD1',
        'step5' = '7M DRD1',
        'step6' = 'Adult DRD1'
        ))+
    xlab('Lineage Step') +
    ylab('met-scDRS') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size = 20)) +
    theme(legend.position = "none")

output.path <- paste0('/u/home/l/lixinzhe/project-cluo/plot/', Sys.Date(), '-BP-DRD1-lineage-step-box-plot.png')
png(
    filename = output.path,
    width = 15,
    height = 10,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

###########################################################################################
######                               met-scDRS with age                              ######
###########################################################################################
# SST THRDE
lineage_step1 = rownames(meta)[meta$fine_age_groups == '2T' & meta$newL2 == 'MGE']
lineage_step2 = rownames(meta)[meta$fine_age_groups == '3T' & meta$newL3 == 'SST-NXPH2']
lineage_step3 = rownames(meta)[meta$fine_age_groups == '1m' & meta$newL3 == 'SST-NXPH2']
lineage_step4 = rownames(meta)[meta$fine_age_groups == '4m' & meta$newL3 == 'SST-NXPH2']
lineage_step5 = rownames(meta)[meta$fine_age_groups == '7m' & meta$newL3 == 'SST-NXPH2']
lineage_step6 = rownames(meta)[meta$fine_age_groups == 'adult' & meta$newL3 == 'SST-NXPH2']

# for each step:
disease = 'PASS_BIP_Mullins2021'
score = risk.score[[disease]]
total_lineage = list(lineage_step1, lineage_step2, lineage_step3, lineage_step4, lineage_step5, lineage_step6)
names(total_lineage) = paste0('step', seq(1, length(total_lineage)))

plot_dfs = NULL
for (step in names(total_lineage)){
    if (length(total_lineage[[step]]) > 0){
        plot_dfs[[step]] = data.frame(
            'score' = score[total_lineage[[step]], ],
            'lineage_step' = step
            )
    }
}
plot_df = Reduce('rbind', plot_dfs)

gplot <- ggplot(plot_df, aes(x = lineage_step, y = score.zscore, fill = lineage_step)) +
    geom_boxplot() +
    # scale_fill_manual(
    #     values = c(
            
    #         '1m' = '#fdc086',
    #         '4m' = '#fc8d62',
    #         '7m' = '#8da0cb',
    #         'adult' = '#66c2a5'
    #         # 'middle temporal gyrus' = '#e78ac3'
    #         )
    #     ) +
    scale_x_discrete(labels = c(
        'step1' = '2T MGE',
        'step2' = '3T SST-NXPH2',
        'step3' = '1M SST-NXPH2',
        'step4' = '4M SST-NXPH2',
        'step5' = '7M SST-NXPH2',
        'step6' = 'Adult SST-NXPH2'
        ))+
    xlab('Lineage Step') +
    ylab('met-scDRS') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size = 20)) +
    theme(legend.position = "none")

output.path <- paste0('/u/home/l/lixinzhe/project-cluo/plot/', Sys.Date(), '-BP-SST-NXPH2-lineage-step-box-plot.png')
png(
    filename = output.path,
    width = 15,
    height = 10,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();
