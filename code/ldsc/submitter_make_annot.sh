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

bed_file=$1
bim_file=$2
annot_file=$3

python make_annot.py \
    --bed-file ${bed_file} \
    --bimfile "${bim_file}" \
    --annot-file "${annot_file}"