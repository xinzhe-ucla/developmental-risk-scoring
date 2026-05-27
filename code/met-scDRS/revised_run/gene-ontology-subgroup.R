### gene-ontology.R ###############################################################################
# purpose: compute gene ontology between groups:

### PREAMBLE ######################################################################################
library(clusterProfiler);
library(dplyr)
library(tibble)
library(ggplot2)
library(grid)

# load in meta:
meta = read.table('/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz', sep = ',', header=TRUE, row.names = 1)
# label the '' as 'unknown'
meta$newL1[meta$newL1 == ''] = 'unknown'
meta$newL1[meta$newL1 == '?'] = 'unknown'
meta$newL2[meta$newL2 == ''] = 'unknown'
meta$newL2[meta$newL2 == '?'] = 'unknown'
meta$newL3[meta$newL3 == ''] = 'unknown'
meta$newL3[meta$newL3 == '?'] = 'unknown'

# load in the methylation matrix:
expr <- RcppCNPy::npyLoad("/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged_X_only.npy")
meta_cells <- read.csv("/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged_cell_meta.csv", row.names = 1)
meta_genes <- read.csv("/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged_gene_meta.csv", row.names = 1)

# format meta cells:
new_name = gsub('-0-0-0', '', rownames(meta_cells))
new_name = gsub('-1-0$', '', new_name)
new_name = gsub('-1$', '', new_name)
new_name = gsub('-1-0-0$', '', new_name)

# format the expression:
rownames(expr) = new_name
colnames(expr) = rownames(meta_genes)

# load in the score:
score = read.table("/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/PASS_Schizophrenia_Pardinas2018.score.gz", sep = '\t', header = TRUE, row.names = 1)
new_name = gsub('-0-0-0', '', rownames(score))
new_name = gsub('-1-0$', '', new_name)
new_name = gsub('-1$', '', new_name)
new_name = gsub('-1-0-0$', '', new_name)
rownames(score) = new_name

### For each group of interest, compute the correlation:
function.path <- '/u/home/l/lixinzhe/project-github/scDRS-applications/spell-book/'
source(paste0(function.path, 'score-loader.R'));
source(paste0(function.path, 'disease-score-expr-correlation.R'))
source(paste0(function.path, 'gene-ontology-caller.R'))
source(paste0(function.path, 'gs-decomposer.R'));

gs.dir <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/';
gsea.c5 <- read.gmt('/u/project/geschwind/lixinzhe/data/c5.all.v2023.1.Hs.entrez.gmt');

trait.gs <- read.table(
    file = paste0(gs.dir, 'magma_10kb_top1000_zscore.75_traits.rv1.gs'),
    sep = '\t',
    header = TRUE
    );

# split gene sets:
trait.gene.set <- lapply(trait.gs$GENESET, gs.decomposer);
names(trait.gene.set) <- trait.gs$TRAIT;

# find out the age dependent network:
meta = meta[rownames(expr), ]
meta$cell_id <- rownames(meta)
cells_of_interest <- meta %>%
  group_by(newL2) %>%
  summarise(cells = list(cell_id), .groups = "drop")

cell_type_of_interest = names(table(meta$newL2[meta$newL1 == 'Exc']))
cells_of_interest = cells_of_interest[cells_of_interest$newL2 %in% cell_type_of_interest,]

# build correlation between each sets of cells X Fine age groups and risk score
frame = cbind(score[rownames(meta), ], meta)

# produce a vector of expr:
cor_vec = rep(NA, ncol(expr))
names(cor_vec) = colnames(expr)

gene_set = risk.score = NULL
# compute the gene expression to score correlation:
for (each_combo in seq(1, nrow(cells_of_interest))){
    # get the set of cells to compute;
    cell_type_label = cells_of_interest[each_combo, 'newL2']
    cell_population = unlist(cells_of_interest[each_combo, 'cells'])
    combo_name = paste0(cell_type_label)
    
    # get universe:
    risk.score[[combo_name]] = score[cell_population, ]
    gene_set[[combo_name]] = trait.gene.set[['PASS_Schizophrenia_Pardinas2018']]
}

score.mch.cor <- expr.ds.cor(
    score = risk.score,
    expression = expr,
    gene.set = gene_set
    );

ontology.gene.num = 100
top.genes.list = {}
for (each_combo in names(gene_set)){
    top.genes.list[[each_combo]] <- names(head(sort(score.mch.cor[[each_combo]], decreasing = TRUE), ontology.gene.num));
}

bg.gene <- colnames(expr);

go_result_list = NULL
for (each_combo in names(gene_set)){
    top.genes = top.genes.list[[each_combo]]
    if (length(top.genes) > 0){
        gene.ontology.result <- gene.ontology.caller(
            x = top.genes,
            background = bg.gene,
            terms = gsea.c5,
            visualize = TRUE
            );
        readable.result <- setReadable(gene.ontology.result, 'org.Hs.eg.db', 'ENTREZID')
        
        # visualize the gene ontology for MDD:
        dot.plot <- dotplot(
            gene.ontology.result,
            showCategory=10
            ) + 
            ggtitle(paste0("gene ontology for ", each_combo)) +
            theme_classic() +
            theme(plot.title = element_text(hjust=0.5)) +
            theme(text = element_text(size = 20))

        
        # plot out the umap:
        dot.path <- '/u/home/l/lixinzhe/project-cluo/plot/dot/'
        combo_name = gsub('/', '_', each_combo)
        plot.name <- paste0(dot.path, Sys.Date(), '-', combo_name, '-gene-ontology-enrichment-dotplot.png')
        png(
            filename = plot.name,
            width = 10,
            height = 10,
            units = 'in',
            res = 400
            );
        print(dot.plot);
        dev.off();
        
        
        # save the gene ontology list:
        go_result_list[[each_combo]] = readable.result
        cor.gene = head(sort(score.mch.cor[[each_combo]], decreasing = TRUE), ontology.gene.num)
        network.plot <- cnetplot(
                readable.result,
                categorySize = "pvalue",
                color.params = list(foldChange = cor.gene, edge = TRUE),
                circular = TRUE) + 
                ggtitle(each_combo) +
                theme_classic() +
                theme(text = element_text(size = 20)) +
                theme(plot.title = element_text(hjust=0.5)) +
                theme(legend.text = element_text(size = 12))

        network.path <- '/u/home/l/lixinzhe/project-cluo/plot/network/'
        combo_name = gsub('/', '_', each_combo)
        plot.name <- paste0(network.path, Sys.Date(), '-', combo_name, '-SCZ-gene-ontology-enrichment-network-legend.png')
        png(
            filename = plot.name,
            width = 10,
            height = 10,
            units = 'in',
            res = 400
            );

        legend <- cowplot::get_legend(network.plot);
        grid.newpage()
        grid.draw(legend)
        dev.off()
        
        plot.name <- paste0(network.path, Sys.Date(), '-', combo_name, '-SCZ-gene-ontology-enrichment-network.png')
        network.plot <- network.plot + theme(legend.position = "none")
        png(
            filename = plot.name,
            width = 10,
            height = 10,
            units = 'in',
            res = 400
            );
        print(network.plot);
        dev.off();
    }

    }
