cp /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/UKB_IDP* /u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/visualize
cp /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/PASS_MDD_Howard2019.score.gz /u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/visualize

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/visualize" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/met-scdrs-brain-volume/"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/april-review/figure-updates/umap-IDP.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/visualize" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/met-scdrs-brain-volume/"
