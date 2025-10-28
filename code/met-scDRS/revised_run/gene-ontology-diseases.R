### GO-enrichment.R ###############################################################################
# purpose: compute the gene ontology enrichment for GSE215353 data

### PREAMBLE ######################################################################################
# load in the libraries:
library(Seurat);
library(sceasy);
library(ggplot2);
library(data.table);
library(ggpubr);
library(clusterProfiler);
library(enrichplot);
library(cowplot);
library(grid);

# define number of top selected genes for computing the go enrichments:
ontology.gene.num <- 100;

# define specified paths:
data.dir <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/';
session.save.path <- '/u/project/pasaniuc/lixinzhe/session-info/';
system.date <- Sys.Date();
save.path <- '/u/project/pasaniuc/lixinzhe/R_saves/';
code.path <- '/u/home/l/lixinzhe/project-github/methylation-RNA-xinzhe-rotation/code/';
plotting.path <- '/u/scratch/l/lixinzhe/tmp-file/tmp-plot/';
function.path <- '/u/home/l/lixinzhe/project-github/scDRS-applications/spell-book/'

# load in the methylation matrix:
mch.fraction <- RcppCNPy::npyLoad("/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged_X_only.npy")
meta_cells <- read.csv("/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged_cell_meta.csv", row.names = 1)
meta_genes <- read.csv("/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged_gene_meta.csv", row.names = 1)

# format meta cells:
new_name = gsub('-0-0-0', '', rownames(meta_cells))
new_name = gsub('-1-0$', '', new_name)
new_name = gsub('-1$', '', new_name)
new_name = gsub('-1-0-0$', '', new_name)

# format the expression:
rownames(mch.fraction) = new_name
colnames(mch.fraction) = rownames(meta_genes)

# define the gene set path:
gs.dir <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/';

# read in the gmt file downloaded from the GSEA website for up to data enrichment estimation:
gsea.c5 <- read.gmt('/u/project/geschwind/lixinzhe/data/c5.all.v2023.1.Hs.entrez.gmt');

# load in meta:
meta = read.table('/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz', sep = ',', header=TRUE, row.names = 1)
# label the '' as 'unknown'
meta$newL1[meta$newL1 == ''] = 'unknown'
meta$newL1[meta$newL1 == '?'] = 'unknown'
meta$newL2[meta$newL2 == ''] = 'unknown'
meta$newL2[meta$newL2 == '?'] = 'unknown'
meta$newL3[meta$newL3 == ''] = 'unknown'
meta$newL3[meta$newL3 == '?'] = 'unknown'

# cehck if all the cels in the data is in the meta:
stopifnot(rownames(mch.fraction) %in% rownames(meta))

# load in the function:
source(paste0(function.path, 'gene-name-converter.R'));
source(paste0(function.path, 'gs-decomposer.R'));
source(paste0(function.path, 'score-loader.R'));
source(paste0(function.path, 'disease-score-expr-correlation.R'))
source(paste0(function.path, 'gene-ontology-caller.R'))

# load in the gene set from the 74 traits gs file:
trait.gs <- read.table(
    file = paste0(gs.dir, 'magma_10kb_top1000_zscore.75_traits.rv1.gs'),
    sep = '\t',
    header = TRUE
    );

# split gene sets:
trait.gene.set <- lapply(trait.gs$GENESET, gs.decomposer);
names(trait.gene.set) <- trait.gs$TRAIT;

# load in the scores:
score <- score.loader('/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/')

for (disease in names(score)){
    # change the index of the score:
    new_name = gsub('-0-0-0', '', rownames(score[[disease]]))
    new_name = gsub('-1-0$', '', new_name)
    new_name = gsub('-1$', '', new_name)
    new_name = gsub('-1-0-0$', '', new_name)
    rownames(score[[disease]]) = new_name
}

### ANALYSIS: GENE ONTOLOGY #######################################################################
#initiate place holder:
traits.interest <- c('PASS_ADHD_Demontis2018', "PASS_BIP_Mullins2021", "PASS_MDD_Howard2019", "PASS_Schizophrenia_Pardinas2018")
top.genes <- vector('list', length = length(traits.interest));
names(top.genes) <- traits.interest;

# get progress bar:
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(traits.interest),
    clear = FALSE
    );

# for each disease call ontology caller:  
go_result_list = NULL      
for (disease in traits.interest) {
    score.mch.cor <- expr.ds.cor(
        score = score[disease],
        expression = mch.fraction,
        gene.set = trait.gene.set[disease]
        );
    
    # extract top genes that are used for gene ontology:
    top.genes[[disease]] <- names(head(sort(score.mch.cor[[disease]], decreasing = TRUE), ontology.gene.num));
    genes = top.genes[[disease]];

    # get the background genes:
    # bg.gene <- intersect(names(trait.gene.set[[disease]]), colnames(mch.fraction));
    bg.gene <- colnames(mch.fraction)
    
    # call GO enrichment analysis:
    gene.ontology.result <- gene.ontology.caller(
        x = genes,
        background = bg.gene,
        terms = gsea.c5,
        visualize = TRUE
        );
    readable.result <- setReadable(gene.ontology.result, 'org.Hs.eg.db', 'ENTREZID')
    
    # save the results:
    saveRDS(score.mch.cor, paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', disease, '_score_mch_correlation.rds'))
    saveRDS(readable.result, paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', disease, '_readable_go_results.rds'))
    
    dot.plot <- dotplot(
        gene.ontology.result,
        showCategory=10
        ) + 
        ggtitle(paste0("gene ontology for ", disease)) +
        theme_classic() +
        theme(plot.title = element_text(hjust=0.5)) +
        theme(text = element_text(size = 20))

        
    # plot out the umap:
    dot.path <- '/u/home/l/lixinzhe/project-cluo/plot/go_disease_dot/'
    plot.name <- paste0(dot.path, Sys.Date(), '-', disease, '-gene-ontology-enrichment-dotplot.png')
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
    go_result_list[[disease]] = readable.result
    cor.gene = head(sort(score.mch.cor[[disease]], decreasing = TRUE), ontology.gene.num)
    network.plot <- cnetplot(
        readable.result,
        categorySize = "pvalue",
        color.params = list(foldChange = cor.gene, edge = TRUE),
        circular = TRUE) + 
        ggtitle(disease) +
        theme_classic() +
        theme(text = element_text(size = 20)) +
        theme(plot.title = element_text(hjust=0.5)) +
        theme(legend.text = element_text(size = 12))

    network.path <- '/u/home/l/lixinzhe/project-cluo/plot/go_disease_network/'
    plot.name <- paste0(network.path, Sys.Date(), '-', disease, '-gene-ontology-enrichment-network-legend.png')
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
        
    plot.name <- paste0(network.path, Sys.Date(), '-', disease, '-gene-ontology-enrichment-network.png')
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
    pb$tick()
    }