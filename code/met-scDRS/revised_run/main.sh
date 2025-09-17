### main.sh #######################################################################################
# call scDRS:
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/submitter.sh"

# split data:
input_gs_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/"
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/parallel-splitter.R \
    --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs" \
    --output_gs "${input_gs_dir}KC_75_traits_split.gs"

for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"
    
    # compute scDRS:
    qsub ${submission_script} \
        "/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad" \
        "${gs_file}" \
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/diagnostics/' \
        '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/' \
        
    # treat the cluster nicely:
    sleep 1

done
