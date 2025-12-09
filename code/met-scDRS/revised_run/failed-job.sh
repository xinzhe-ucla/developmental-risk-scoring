# one job failed due to not enough time, recompute it

# identify the job that failed
gs_file="/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs"
result_dir='/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/'

while read trait _; do
    if [[ ! -f "${result_dir}/${trait}.full_score.gz" ]]; then
        echo "$trait"
    fi
done < $gs_file

# turns out UKB_460K.repro_NumberChildrenEverBorn_Pooled is never computed
cd /u/home/l/lixinzhe/project-cluo_scratch/lixinzhe
head -n1 $gs_file > tmp_gs.gs
grep 'UKB_460K.repro_NumberChildrenEverBorn_Pooled' $gs_file >> tmp_gs.gs

submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/submitter.sh"
qsub ${submission_script} \
    "/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad" \
    "/u/home/l/lixinzhe/project-cluo_scratch/lixinzhe/tmp_gs.gs" \
    "mean_var_length" \
    "arcsine" \
    "inv_std" \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/diagnostics/cov/' \
    '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/' \
    '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov'
