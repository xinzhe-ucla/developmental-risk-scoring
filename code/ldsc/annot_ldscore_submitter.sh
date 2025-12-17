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

bfile=$1
anno=$2
ldscore_out=$3

for annot in ${anno_out_dir}*.annot.gz; do
	python ldsc.py \
		--l2 \
		--bfile ${bfile} \
		--ld-wind-cm 1 \
		--annot ${anno} \
		--thin-annot \
		--out "${ldscore_out}" \
		--print-snps /u/home/l/lixinzhe/project-geschwind/software/ldsc/hm3_no_MHC.list.txt
done