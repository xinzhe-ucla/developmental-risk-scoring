### for visualizing proportion plot
current_date=$(date +"%Y-%m-%d")

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/proportion-heatmap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv" \
    --field 'L1' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-initial-BICAN-mcg-l1-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/proportion-heatmap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv" \
    --field 'L2' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-initial-BICAN-mcg-l2-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/proportion-heatmap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv" \
    --field 'L3' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-initial-BICAN-mcg-L3-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/proportion-heatmap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv" \
    --field 'fine_age' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-initial-BICAN-mcg-fine-age-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/proportion-heatmap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv" \
    --field 'Region' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-initial-BICAN-mcg-region-proportion.png" \
    --plot_type "proportion"

### for umap:
current_date=$(date +"%Y-%m-%d")

Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE215353/met-scDRS-v2.0-run/visualization/significant-cell-visualization-umap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data '/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv' \
    --xaxis "UMAP_1" \
    --yaxis "UMAP_2" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/intial_run_developmental_methylation_umap/"

### add a UMAP of cell type:
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/cell_type_plot.R