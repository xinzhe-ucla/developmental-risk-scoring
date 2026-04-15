### submission.sh #################################################################################
# call scDRS:
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/april-review/intermediate-submission.sh"

# split data:
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"

head -n1 /u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs > /u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/scz_score.gs
grep 'PASS_Schizophrenia_Pardinas2018' /u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs >> /u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/scz_score.gs
gs_file="/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/scz_score.gs"

# compute scDRS:
qsub ${submission_script} \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad' \
    "${gs_file}" \
    "mean_var_length" \
    "arcsine" \
    "inv_std" \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/diagnostics/cov/' \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-arcsine-inv_std/' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_intermediate.pkl'
    
# compute scDRS:
qsub ${submission_script} \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad' \
    "${gs_file}" \
    "mean_var" \
    "arcsine" \
    "inv_std" \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/diagnostics/cov/' \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-arcsine-inv_std/' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov'  \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_intermediate.pkl'

# compute scDRS:
qsub ${submission_script} \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad' \
    "${gs_file}" \
    "mean_var" \
    "logit" \
    "inv_std" \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/diagnostics/cov/' \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-logit-inv_std/' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov'  \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_intermediate.pkl'

# compute scDRS:
qsub ${submission_script} \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad' \
    "${gs_file}" \
    "mean_var" \
    "library" \
    "inv_std" \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/diagnostics/cov/' \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var-library-inv_std/' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_intermediate.pkl'

# also investigate mean var length:
qsub ${submission_script} \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad' \
    "${gs_file}" \
    "mean_var_length" \
    "library" \
    "inv_std" \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/diagnostics/cov/' \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-library-inv_std/' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_intermediate.pkl'

# compute scDRS:
qsub ${submission_script} \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad' \
    "${gs_file}" \
    "mean_var_length" \
    "logit" \
    "inv_std" \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/diagnostics/cov/' \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/revision-stability/mean_var_length-logit-inv_std/' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_intermediate.pkl'
