# mean var logit:
cp /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/PASS_MDD_Howard2019.score.gz /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-logit-inv_std/
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-logit-inv_std/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/mean_var_logit"

# mean var arcsine:
cp /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/PASS_MDD_Howard2019.score.gz /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-arcsine-inv_std/
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-arcsine-inv_std/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/mean_var_arcsine"

# mean var library:
cp /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/PASS_MDD_Howard2019.score.gz /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-library-inv_std/
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-library-inv_std/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/mean_var_library"

### Look at mean var length logit:
cp /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/PASS_MDD_Howard2019.score.gz /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-logit-inv_std/
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-logit-inv_std/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/mean_var_length_logit"

# mean var length arcsine
cp /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/PASS_MDD_Howard2019.score.gz /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-arcsine-inv_std
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-arcsine-inv_std" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/mean_var_length_arcsine"

# mean var length library
cp /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/PASS_MDD_Howard2019.score.gz /u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-library-inv_std/
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/umap-plot.R \
    --dir "/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-library-inv_std/" \
    --meta_data '/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv' \
    --xaxis "X_umap" \
    --yaxis "Y_umap" \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/mean_var_length_library"
