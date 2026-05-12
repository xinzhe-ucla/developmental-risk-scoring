# download E081 and E082
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc/fetal_brain_dnase_baseline/
wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E081-DNase.hotspot.fdr0.01.peaks.v2.bed.gz 
wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E082-DNase.hotspot.fdr0.01.peaks.v2.bed.gz

# corresponds to fetal male and female annotation (https://egg2.wustl.edu/roadmap/web_portal/meta.html)

#-----------------------------------
# merge male and female
#-----------------------------------
module load bedtools/2.23.0
zcat E081-DNase.hotspot.fdr0.01.peaks.v2.bed.gz \
     E082-DNase.hotspot.fdr0.01.peaks.v2.bed.gz \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  | gzip > fetal_brain_DNase_hotspot_fdr0.01_union.bed.gz

#-----------------------------------
# LDSC
#-----------------------------------
# making annotation:
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc

submission_script='/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/submitter_make_annot.sh'
anno_out_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/fetal_brain_baseline/annot/'
for hg19_DMR_bed in /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc/fetal_brain_dnase_baseline/fetal_brain_DNase_hotspot_fdr0.01_union.bed.gz; do
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
anno_out_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/fetal_brain_baseline/annot/'
for hg19_DMR_bed in /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc/fetal_brain_dnase_baseline/fetal_brain_DNase_hotspot_fdr0.01_union.bed.gz; do
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
