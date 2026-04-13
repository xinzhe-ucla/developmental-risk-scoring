###########################################################################################
######                          Download summary statistics                          ######
###########################################################################################
# download the brain related traits:
cd /u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0124.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0125.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0126.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0127.txt.gz

curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0013.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0014.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0015.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0016.txt.gz

###########################################################################################
######                                    Create step1                               ######
###########################################################################################
# specify common file paths:
magma_dir="/u/home/l/lixinzhe/project-geschwind/software/magma"

# specify input file:
${magma_dir}/magma \
    --annotate window=10,10 \
    --snp-loc ${magma_dir}/aux/g1000_eur.bim \
    --gene-loc ${magma_dir}/aux/NCBI37.3.gene.loc \
    --out ${magma_dir}/aux/step1

out_magma_dir="/u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas/magma_output/"

###########################################################################################
######           preprocess the summary statistics and creates magma file        ######
###########################################################################################
for item in /u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas/*.txt.gz; do
    qsub /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/april-review/brain-imaging-traits/submitter.sh "$item"
done

# piece the resulting scores together:
i=0
for item in /u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas/magma_output/*-munge-output.txt; do
    echo "$i $item"
    i=$((i + 1))
    if [ "$i" -eq 1 ]; then
        head -n2 $item > "$out_magma_dir/ukb-brain-volume-traits.gs"
    else
        tail -n1 $item >> "$out_magma_dir/ukb-brain-volume-traits.gs"
    fi
done

###########################################################################################
######                                  Compute 8 traits                             ######
###########################################################################################
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/submitter.sh"

for gs_file in /u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas/magma_output/*-munge-output.txt; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"
    
    # submit:
    qsub ${submission_script} \
        "/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad" \
        "${gs_file}" \
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/diagnostics/cov/' \
        '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/' \
        '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov'
done

###########################################################################################
######                                    output an intermediate                     ######
###########################################################################################
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/april-review/intermediate-submission.sh"

for gs_file in /u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas/magma_output/0013*-munge-output.txt; do
    # for each of the gs file submit a job:
    echo "read gs file:"
    echo "$gs_file"
    
    # submit:
    qsub ${submission_script} \
        "/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad" \
        "${gs_file}" \
        "mean_var_length" \
        "arcsine" \
        "inv_std" \
        '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/diagnostics/cov/' \
        '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/cov/' \
        '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov' \
        '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_intermediate.pkl'
done