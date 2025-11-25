### processing script for GO is located at:
'/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/gene-ontology-age.R'


# load in the work for generating network plot:
library(Seurat);
library(sceasy);
library(ggplot2);
library(data.table);
library(ggpubr);
library(clusterProfiler);
library(enrichplot);
library(cowplot);
library(grid);

cor_collection = NULL
go_collection = NULL

traits.interest <- c('PASS_ADHD_Demontis2018', "PASS_Schizophrenia_Pardinas2018")
combos = c('1m', '2T', '3T', '4-7m', 'adult')

for (disease in c('SCZ', 'ADHD')){
    for (combo in combos){
        name = paste0(combo, '_', disease)
        cor_collection[[name]] = readRDS(paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/age/', name, '_score_mcg_correlation.rds'))
        go_collection[[name]] = readRDS(paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/age/', name, '_readable_go_results.rds'))
    }
}

###########################################################################################
######                                      PLOTTING                                 ######
###########################################################################################
# save the gene ontology list:
ontology.gene.num = 100

for (combo in names(cor_collection)){
    readable.result = go_collection[[combo]]
    cor.gene = head(sort(cor_collection[[combo]], decreasing = TRUE), ontology.gene.num)
    network.plot <- cnetplot(
        readable.result,
        showCategory = 5,
        categorySize = "pvalue",
        color.params = list(foldChange = cor.gene, edge = TRUE),
        circular = TRUE) + 
        ggtitle(combo) +
        theme_classic() +
        theme(text = element_text(size = 25)) +
        theme(plot.title = element_text(hjust=0.5)) +
        theme(legend.text = element_text(size = 12))

    network.path <- '/u/home/l/lixinzhe/project-geschwind/plot/'
    plot.name <- paste0(network.path, Sys.Date(), '-', combo, '-gene-ontology-enrichment-network-legend.png')
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
        
    plot.name <- paste0(network.path, Sys.Date(), '-', combo, '-gene-ontology-enrichment-network.png')
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
