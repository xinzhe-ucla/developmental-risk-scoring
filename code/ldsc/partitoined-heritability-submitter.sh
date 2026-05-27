#!/bin/bash
#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=00:30:00,h_data=10G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 1
# Email address to notify
#$ -M lxzjason@gmail.com
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
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc

gwas=$1
annot_ldscore=$2
baseline_ldscore=$3
output=$4
python ldsc.py \
	--h2 "${gwas}" \
	--w-ld-chr /u/home/l/lixinzhe/project-geschwind/software/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--ref-ld-chr ${annot_ldscore},${baseline_ldscore} \
	--overlap-annot \
	--frqfile-chr /u/home/l/lixinzhe/project-geschwind/software/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
	--out ${output} \
	--print-coefficients

# python ldsc.py \
# 	--h2 "${gwas_dir}PGC_UKB_depression_genome_wide.sumstats.gz" \
# 	--w-ld-chr /u/home/l/lixinzhe/project-geschwind/software/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
# 	--ref-ld-chr ${annotation_dir}adult_NN-NN-VLMC.hypo_dmr_overlap.hg19.dmr.,${baseline_dir}baselineLD. \
# 	--overlap-annot \
# 	--frqfile-chr /u/home/l/lixinzhe/project-geschwind/software/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
# 	--out /u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/adult_NN-NN-VLMC.hypo_dmr_overlap \
# 	--print-coefficients