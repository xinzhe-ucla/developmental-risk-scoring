###########################################################################################
######                                 Running liftover                              ######
###########################################################################################
liftover_chain="/u/home/l/lixinzhe/project-geschwind/software/liftOver/aux_file/hg38ToHg19.over.chain"
mkdir -p /u/scratch/l/lixinzhe/tmp-file/DMR0526/
conda activate default_r_base
n=191
i=0
for DMR_bed in /u/home/h/hex002/project-cluo/BICAN/loop_DMR_0526/*.hypo_dmr_overlap.bed; do
    i=$((i+1))
    percent=$(( i * 100 / n ))
    filled=$(( percent / 5 ))
    empty=$(( 20 - filled ))
    printf "\r[%.*s%*s] %3d%% (%d/%d)" \
        "$filled" "####################" \
        "$empty" "" \
        "$percent" "$i" "$n"
    
    base=$(basename "$DMR_bed")
    base=${base%.*}
    hg19_DMR_bed="/u/scratch/l/lixinzhe/tmp-file/DMR0526/${base}.hg19.dmr.bed"
    CrossMap bed $liftover_chain $DMR_bed $hg19_DMR_bed 
done

###########################################################################################
######                                Running make annot                             ######
###########################################################################################
# making annotation:
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc

submission_script='/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/submitter_make_annot.sh'
anno_out_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/dmr0526-annot/'
mkdir -p $anno_out_dir
for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/DMR0526/*.hg19.dmr.bed; do
    # make annotation from chromosome 1:22
    for CHR in {1..22}; do
        base=$(basename "$hg19_DMR_bed")
        base=${base%.*}
        out="${anno_out_dir}${base}.${CHR}.annot.gz"
        if [[ -f "$out" ]]; then
            echo "[$base chr$CHR] exists -> skip"
            continue
        fi
        qsub ${submission_script} \
            "${hg19_DMR_bed}" \
            "/u/home/l/lixinzhe/project-geschwind/software/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR}.bim" \
            "${anno_out_dir}${base}.${CHR}.annot.gz"
    done
done

###########################################################################################
######                                    Running LDSC                               ######
###########################################################################################
# making annotation:
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc

submission_script='/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/annot_ldscore_submitter.sh'
anno_out_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/dmr0526-annot/'
for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/DMR0526/*.hg19.dmr.bed; do
    # make annotation from chromosome 1:22
    for CHR in {1..22}; do
        base=$(basename "$hg19_DMR_bed")
        base=${base%.*}
        annotation_file="${anno_out_dir}${base}.${CHR}.annot.gz"
        bfile="/u/home/l/lixinzhe/project-geschwind/software/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR}"
        annotation_ldscore="${anno_out_dir}${base}.${CHR}"
        if [[ -f "${annotation_ldscore}.l2.ldscore.gz" ]]; then
            echo "[$base chr$CHR] exists -> skip"
            continue
        fi
        qsub ${submission_script} \
            "${bfile}" \
            "${annotation_file}" \
            "${annotation_ldscore}"
    done
done

###########################################################################################
######                       LDSC Partitioned Heritability MDD                       ######
###########################################################################################
# do ldsc:
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc
annotation_dir="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/dmr0526-annot/"
baseline_dir="/u/home/l/lixinzhe/project-geschwind/software/ldsc/"

# for MDD:
gwas_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/gwas/'
gwas="${gwas_dir}PGC_UKB_depression_genome_wide.sumstats.gz"
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/partitoined-heritability-submitter.sh"
baselines="/u/home/l/lixinzhe/project-geschwind/software/ldsc/baselineLD."

for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/DMR0526/*.hg19.dmr.bed; do
    base=$(basename "$hg19_DMR_bed")
    annot_ldscore="${annotation_dir}${base%.*}."
    output="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/dmr0526/${base%.*}"
    
    if [[ -f "${output}.results" ]]; then
        echo "output exists -> skip"
        continue
    fi
    qsub ${submission_script} \
        "${gwas}" \
        "${annot_ldscore}" \
        "${baselines}" \
        "${output}"
done

###########################################################################################
######                  LDSC Partitioned Heritability SCZ                            ######
###########################################################################################
gwas_dir="/u/home/l/lixinzhe/project-geschwind/data/GWAS/"
gwas="${gwas_dir}Schizophrenia_pardinas2018_munged.sumstats.gz"
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/partitoined-heritability-submitter.sh"
baselines="/u/home/l/lixinzhe/project-geschwind/software/ldsc/baselineLD."

for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/DMR0526/*.hg19.dmr.bed; do
    base=$(basename "$hg19_DMR_bed")
    annot_ldscore="${annotation_dir}${base%.*}."
    output="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/dmr0526/Schizophrenia_pardinas2018/${base%.*}"
    
    if [[ -f "${output}.results" ]]; then
        echo "output exists -> skip"
        continue
    fi
    qsub ${submission_script} \
        "${gwas}" \
        "${annot_ldscore}" \
        "${baselines}" \
        "${output}"
done

###########################################################################################
######                  LDSC Partitioned Heritability ADHD                           ######
###########################################################################################
gwas_dir="/u/home/l/lixinzhe/project-geschwind/data/GWAS/"
gwas="${gwas_dir}ADHD_Demontis2018_munged.sumstats.gz"
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/partitoined-heritability-submitter.sh"
baselines="/u/home/l/lixinzhe/project-geschwind/software/ldsc/baselineLD."

for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/DMR0526/*.hg19.dmr.bed; do
    base=$(basename "$hg19_DMR_bed")
    annot_ldscore="${annotation_dir}${base%.*}."
    output="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/dmr0526/ADHD_Demontis2018/${base%.*}"
    
    if [[ -f "${output}.results" ]]; then
        echo "output exists -> skip"
        continue
    fi
    qsub ${submission_script} \
        "${gwas}" \
        "${annot_ldscore}" \
        "${baselines}" \
        "${output}"
done

###########################################################################################
######                                    Visualization                              ######
###########################################################################################

# visualize
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/visualization.R \
    --ldsc_dir '/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/dmr0526/' \
    --output "/u/home/l/lixinzhe/project-geschwind/plot/$(date +%F)-dmr0526-mdd-heritability-developmental-cell-type-by-time-point.pdf"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/visualization.R \
    --ldsc_dir '/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/dmr0526/Schizophrenia_pardinas2018/' \
    --output "/u/home/l/lixinzhe/project-geschwind/plot/$(date +%F)-dmr0526-scz-pardinas2018-heritability-developmental-cell-type-by-time-point.pdf"

# visualize:
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/visualization.R \
    --ldsc_dir '/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/dmr0526/ADHD_Demontis2018/' \
    --output "/u/home/l/lixinzhe/project-geschwind/plot/$(date +%F)-dmr0526-ADHD_Demontis2018-heritability-developmental-cell-type-by-time-point.pdf"
