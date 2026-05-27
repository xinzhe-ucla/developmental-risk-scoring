###########################################################################################
######                                 Running liftover                              ######
###########################################################################################
liftover_chain="/u/home/l/lixinzhe/project-geschwind/software/liftOver/aux_file/hg38ToHg19.over.chain"
conda activate default_r_base
n=8
i=0
for DMR_bed in /u/project/cluo/chongyua/brain_dev_snm3C/2025/pairwise_dmr/*.bed; do
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
    hg19_DMR_bed="/u/scratch/l/lixinzhe/tmp-file/DMR/pairwise/${base}.hg19.dmr.bed"
    CrossMap bed $liftover_chain $DMR_bed $hg19_DMR_bed 
done

###########################################################################################
######                                Running make annot                             ######
###########################################################################################
# making annotation:
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc

submission_script='/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/submitter_make_annot.sh'
anno_out_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/dmr_pairwise_annot/'
for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/DMR/pairwise/*.hg19.dmr.bed; do
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
anno_out_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/dmr_pairwise_annot/'
for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/DMR/pairwise/*.hg19.dmr.bed; do
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
annotation_dir="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/dmr_pairwise_annot/"
baseline_dir="/u/home/l/lixinzhe/project-geschwind/software/ldsc/"

# for MDD:
gwas_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/gwas/'
gwas="${gwas_dir}PGC_UKB_depression_genome_wide.sumstats.gz"
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/pair-wise-heritability.sh"
baselines="/u/home/l/lixinzhe/project-geschwind/software/ldsc/baselineLD."

for age_point in 1m 2T 3T 4-7m adult; do
    base="methylpy_pair_${age_point}_dms2_hypo_${age_point}_Exc_UL_L4-RORB.hg19.dmr."
    annot_ldscore="${annotation_dir}${base%.*}."
    
    background="methylpy_pair_${age_point}_dms2_hypo_${age_point}_Inh_MSN_DRD1-eccentric-CASZ1.hg19.dmr."
    background_ldscore="${annotation_dir}${background%.*}."
    
    additional_baseline="/u/home/l/lixinzhe/project-geschwind/result/ldsc/fetal_brain_baseline/annot/fetal_brain_DNase_hotspot_fdr0.01_union.bed."
    
    output="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/pairwise_with_fetal_brain/${base%.*}"
    if [[ -f "${output}.results" ]]; then
        echo "output exists -> skip"
        continue
    fi
    qsub ${submission_script} \
        "${gwas}" \
        "${annot_ldscore}" \
        "${baselines}" \
        "${background_ldscore}" \
        "${output}" \
        "${additional_baseline}"
done
