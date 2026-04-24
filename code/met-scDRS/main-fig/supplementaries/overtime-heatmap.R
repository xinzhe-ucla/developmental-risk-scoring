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
brain_traits = trait.info$Trait_Identifier[trait.info$Category == 'brain']

### look at the age specific data:
age_point = unique(meta$fine2_age_groups)
for (age in age_point){
    # cells at age point:
    cells_at_age <- rownames(meta)[meta$fine2_age_groups == age]
    
    # subset:
    meta_at_age = meta[cells_at_age, ]
    
    # get the proportion of cells that are significant in each of the disease x cell type
    significance_matrix = data.frame(matrix(NA, nrow = length(risk.score), ncol = length(unique(meta_at_age$adjusted_L3))))
    colnames(significance_matrix) = unique(meta_at_age$adjusted_L3)
    rownames(significance_matrix) = names(risk.score)
    num_matrix = significance_matrix
    
    # for each disease, get the proportion:
    for (disease in names(risk.score)){
        age_only_risk_score = risk.score[[disease]][cells_at_age, 'zscore']
        age_fdr = risk.score[[disease]][cells_at_age, 'fdr']
        
        # bind the meta to the risk score
        meta_at_age$risk_score = age_only_risk_score
        meta_at_age$fdr = age_fdr
        
        # calculate proportion:
        proprotion_sig = meta_at_age %>% group_by(adjusted_L3) %>% 
            summarize(
                n = n(),
                n_sig = sum(fdr < 0.1, na.rm = TRUE),
                proportion = n_sig / n,
                .groups = "drop"
            )
        
        # coerce into named vector:
        sig_by_celltype_in_disease = c(proprotion_sig$proportion)
        number_by_celltype_in_disease = c(proprotion_sig$n)
        names(sig_by_celltype_in_disease) = proprotion_sig$adjusted_L3
        names(number_by_celltype_in_disease) = proprotion_sig$adjusted_L3
        sig_by_celltype_in_disease = sig_by_celltype_in_disease[colnames(significance_matrix)]
        number_by_celltype_in_disease = number_by_celltype_in_disease[colnames(significance_matrix)]
        # put the result into a significant proportion matrix:
        significance_matrix[disease, ] = sig_by_celltype_in_disease
        num_matrix[disease, ] = number_by_celltype_in_disease
    }
    
    ### make the matrix ###
    significance.matrix = significance_matrix
    write.table(
        significance.matrix,
        file = paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-cell-type-brain-traits-proportion-', age, '-only.csv'),
        sep = ',',
        col.names = TRUE,
        row.names = TRUE,
        quote = FALSE
        )

    # find out the set of brain traits:
    trait.info$Trait_Identifier <- gsub('PASS_', '', trait.info$Trait_Identifier)
    trait.info$Trait_Identifier <- gsub('UKB_460K.', '', trait.info$Trait_Identifier)
    trait.info$Trait_Identifier <- gsub('cov_', '', trait.info$Trait_Identifier)
    trait.info$Trait_Identifier <- gsub('repro_', '', trait.info$Trait_Identifier)

    rownames(num_matrix) <- rownames(significance.matrix) <- gsub('PASS_', '', rownames(significance.matrix))
    rownames(num_matrix) <- rownames(significance.matrix) <- gsub('UKB_460K.', '', rownames(significance.matrix))
    rownames(num_matrix) <- rownames(significance.matrix) <- gsub('cov_', '', rownames(significance.matrix))
    rownames(num_matrix) <- rownames(significance.matrix) <- gsub('repro_', '', rownames(significance.matrix))

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
    
    # grab out the cell types:
    cell.types = unique(meta_at_age[, 'adjusted_L3'])
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
    B <- num_matrix[publication.traits, cell.type.order]
    min_cells_in_bin = 0
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
        heatmap_legend_param = heatmap.legend.param,
        cell_fun = function(j, i, x, y, w, h, fill) {
            
            if (B[i, j] <= min_cells_in_bin) {
            # grey out
            grid.rect(
                x, y, w, h,
                gp = gpar(fill = "grey85", col = NA)
            )
            } else {
            # normal heatmap cell
            grid.rect(
                x, y, w, h,
                gp = gpar(fill = fill, col = NA)
            )
            }
        }
        );
    plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 10, 70), "mm"));

    # measure the size of the heatmap:
    heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
    heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

    # use the measured width and height for drawing:
    output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-developmental-cell-type-brain-traits-proportion-', age, '-only.pdf')
    pdf(
        file = output.path,
        width = heatmap.width,
        height = heatmap.height
        );
    draw(plot, heatmap_legend_side = 'left', padding = unit(c(30, 10, 10, 70), "mm"));
    dev.off();
}



