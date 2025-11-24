### processing script for GO is located at:
'/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/gene-ontology-diseases.R'


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

for (disease in traits.interest){
    cor_collection[[disease]] = readRDS(paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', disease, '_score_mch_correlation.rds'))
    go_collection[[disease]] = readRDS(paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', disease, '_readable_go_results.rds'))
}

###########################################################################################
######                                      PLOTTING                                 ######
###########################################################################################
# save the gene ontology list:
ontology.gene.num = 100

for (disease in traits.interest){
    readable.result = go_collection[[disease]]
    cor.gene = head(sort(cor_collection[[disease]][[disease]], decreasing = TRUE), ontology.gene.num)
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

    network.path <- '/u/home/l/lixinzhe/project-geschwind/plot/'
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

}
