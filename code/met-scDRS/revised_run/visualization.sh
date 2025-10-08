current_date=$(date +"%Y-%m-%d")

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'newL1' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-l1-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'newL2' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-l2-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'newL3' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-l3-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'fine_age_groups' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-fine-age-proportion.png" \
    --plot_type "proportion"

## UMAP:
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/' \
    --umap_coord '/u/project/cluo/heffel/BICAN3/3C/bican_2025_3Cumap.csv.gz' \
    --UMAP1 'X_umap' \
    --UMAP2 'Y_umap' \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/revised_dev_methylome/${current_date}-revised-BICAN-mcg-fine-age-proportion.png"
    
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/revised-intial_run_developmental_methylation_umap-281146/"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/exc/" \
    --meta_data '/u/project/cluo/heffel/BICAN3/3C/bican_2025_3Cumap.csv.gz' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/revised-intial_run_developmental_methylation_umap/exc/"

###########################################################################################
######                     Visualize with covariates regressed                       ######
###########################################################################################
current_date=$(date +"%Y-%m-%d")

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'newL1' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-l1-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'newL2' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-l2-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'newL3' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-l3-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'fine_age_groups' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-fine-age-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/cov-intial_run_developmental_methylation_umap-281146/"
