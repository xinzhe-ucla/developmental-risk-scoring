#!/bin/bash
#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=02:00:00,h_data=3G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 4
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
PATH=$PATH:$HOME/bin:/u/project/pasaniuc/lixinzhe/software/nextflow
export PATH
module load gcc
module load intel
module load java/jdk-11.0.14
module load anaconda3
conda activate default_r_base

item=$1
out_magma_dir="/u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas/magma_output/"
magma_dir="/u/home/l/lixinzhe/project-geschwind/software/magma"

echo "processing $item"
start=$(date +%s)
base=$(basename "$item" .txt.gz)

if [ -f "${out_magma_dir}${base%.txt.gz}-grch37-munge-input.txt" ]
then
    echo "${out_magma_dir}${base%.txt.gz}-grch37-munge-input.txt output already exists, skipping"
else
    # add the number of sample into the data:
    zcat "$item" \
    | awk 'BEGIN{OFS="\t"}
            NR==1 {
                $(NF+1)="N"
                $(NF+1)="pval"
                print
                next
            }
            {
                logp=$8
                p=10^(-logp)
                $(NF+1)="33219"
                $(NF+1)=p
                print
            }' > "/u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas/withN/${base%.txt.gz}_withN.txt"

    # call magma:
    /u/home/l/lixinzhe/project-geschwind/software/magma/magma \
        --bfile ${magma_dir}/aux/g1000_eur \
        --pval "/u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas/withN/${base%.txt.gz}_withN.txt" use='rsid,pval' ncol='N' \
        --gene-annot ${magma_dir}/aux/step1.genes.annot \
        --out ${out_magma_dir}${base%.txt.gz}-grch37
        
    # get time:
    end=$(date +%s)
    echo "Elapsed: $((end - start)) seconds"

    # prepare the z-score file:
    Rscript "/u/home/l/lixinzhe/project-github/scDRS-applications/code/scDRS-pipeline/get-magma-stats.R" \
        "${out_magma_dir}${base%.txt.gz}-grch37.genes.out" \
        "${magma_dir}/aux/NCBI37.3.gene.loc" \
        "${out_magma_dir}${base%.txt.gz}-grch37-munge-input.txt"

    # add new header:
    new_header="gene\tUKB_IDP${base%.txt.gz}"
    sed -i "1 s/.*/$new_header/" "${out_magma_dir}${base%.txt.gz}-grch37-munge-input.txt"
fi

# deploy gs munge:
scdrs munge-gs \
    --out-file "${out_magma_dir}${base%.txt.gz}-grch37-munge-output.txt" \
    --zscore-file "${out_magma_dir}${base%.txt.gz}-grch37-munge-input.txt" \
    --weight zscore \
    --n-max 1000