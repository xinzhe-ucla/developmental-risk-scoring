### main.sh #######################################################################################
# purpose: initial compute met-scDRS for merged data

### single cell level baseline ###
# without any sort of batch correction or QC or aggregation:
# call scDRS:
input_gs_dir="/u/scratch/l/lixinzhe/tmp-file/job-array/"
out_dir="/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/"
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/submitter.sh"
h5ad_file='/u/home/l/lixinzhe/project-cluo/data/2025-04-30-combined-mcg-processed.h5ad'

# split data:
Rscript /u/home/l/lixinzhe/project-github/scDRS-applications/code/met-scDRS-method/version-2.0/parallel-splitter.R \
    --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs" \
    --output_gs "${input_gs_dir}KC_75_traits_split.gs"

for gs_file in ${input_gs_dir}KC_75_traits_split.gs*; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"

    # compute scDRS:
    qsub ${submission_script} \
        "${gs_file}" \
        "${h5ad_file}" \
        "${out_dir}"

    # treat the cluster nicely:
    sleep 1

done
