#!/bin/bash
#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/l/lixinzhe/job-log/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=25G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 5
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
conda activate default_r_base

h5ad_file="/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged.h5ad"
control_scheme="mean_var_length"
tmp_gs="/u/home/l/lixinzhe/project-geschwind/port/scratch/parallel_gs/KC_75_traits_split.gs1"
output_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/"
transform_scheme="arcsine"
weight_opt="inv_std"
diagnostic_dir="/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/"
cov_file="/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov"

met_scdrs compute_score \
    --h5ad_file ${h5ad_file} \
    --preprocess True \
    --preprocess_method inverse \
    --variance_clip 5 \
    --transformation ${transform_scheme} \
    --h5ad_species human \
    --cov_file ${cov_file} \
    --gs-file ${tmp_gs} \
    --gs_species human \
    --out_folder ${output_dir} \
    --ctrl_match_opt ${control_scheme} \
    --weight_opt ${weight_opt} \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True \
    --diagnostic True \
    --diagnostic_dir ${diagnostic_dir} \
    --verbose True \
    --intermediate "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/diagnostics/regressed_unnormalized_genebody_blacklist_allgenes_merged.h5ad"
