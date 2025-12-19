###########################################################################################
######                                       PREAMBLE                                ######
###########################################################################################
# first lets download the files needed to run partitioned heritability:
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc
wget https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.tar.gz # downloaded as file-archive zip
unzip files-archive

# inflate the tar file for 1000G genotype;
tar -xvzf 1000G_Phase3_plinkfiles.tgz

# download liftover:
cd /u/home/l/lixinzhe/project-geschwind/software/liftOver
rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/ ./
chmod +x liftOver
mkdir aux_file
cd aux_file

wget --timestamping 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' -O hg38ToHg19.over.chain.gz
gunzip hg38ToHg19.over.chain.gz

###########################################################################################
######                                 Running liftover                              ######
###########################################################################################
liftover_chain="/u/home/l/lixinzhe/project-geschwind/software/liftOver/aux_file/hg38ToHg19.over.chain"
conda activate default_r_base
n=174
i=0
for DMR_bed in /u/home/h/hex002/project-cluo/BICAN/loop_DMR/*.hypo_dmr_overlap.bed; do
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
    hg19_DMR_bed="/u/scratch/l/lixinzhe/tmp-file/${base}.hg19.dmr.bed"
    CrossMap bed $liftover_chain $DMR_bed $hg19_DMR_bed 
done
    

###########################################################################################
######                                Running make annot                             ######
###########################################################################################
# making annotation:
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc

submission_script='/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/submitter_make_annot.sh'
anno_out_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/annot/'
for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/*.hg19.dmr.bed; do
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
anno_out_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/annot/'
for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/*.hg19.dmr.bed; do
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
######                                    Munge SUMSTAT                              ######
###########################################################################################
# speicify genome wide summary statistics:
gwas='/u/home/l/lixinzhe/project-geschwind/data/MDD-GWAS/PGC_UKB_depression_genome-wide.txt'
gwas_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/gwas/'

## munge MDD statistics:
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc 
python munge_sumstats.py \
    --snp "MarkerName" \
    --N-cas 246363 \
    --N-con 561190 \
    --sumstats ${gwas} \
    --a1 A1 \
    --a2 A2 \
    --frq Freq \
    --p P \
    --signed-sumstats LogOR,0 \
    --out "${gwas_dir}PGC_UKB_depression_genome_wide"

### SCZ pardinaz 2018
cd /u/home/l/lixinzhe/project-geschwind/data/GWAS
gwas_dir="/u/home/l/lixinzhe/project-geschwind/data/GWAS/"
gwas="/u/home/l/lixinzhe/project-geschwind/data/GWAS/Schizophrenia_pardinas2018"

conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc 
python munge_sumstats.py \
    --snp "SNP" \
    --N-cas 11260 \
    --N-con 24542 \
    --sumstats ${gwas} \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --out "${gwas_dir}Schizophrenia_pardinas2018_munged"

### ADHD 
gwas_dir="/u/home/l/lixinzhe/project-geschwind/data/GWAS/"
gwas="/u/home/l/lixinzhe/project-geschwind/data/GWAS/ADHD_Demontis2018"

cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc 
python munge_sumstats.py \
    --snp "SNP" \
    --N-cas-col Nca \
    --N-con-col Nco \
    --sumstats ${gwas} \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --out "${gwas_dir}ADHD_Demontis2018_munged"

###########################################################################################
######                       LDSC Partitioned Heritability MDD                       ######
###########################################################################################
# do ldsc:
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc
annotation_dir="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/annot/"
baseline_dir="/u/home/l/lixinzhe/project-geschwind/software/ldsc/"

# for MDD:
gwas_dir='/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/gwas/'
gwas="${gwas_dir}PGC_UKB_depression_genome_wide.sumstats.gz"
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/partitoined-heritability-submitter.sh"
baselines="/u/home/l/lixinzhe/project-geschwind/software/ldsc/baselineLD."

for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/*.hg19.dmr.bed; do
    base=$(basename "$hg19_DMR_bed")
    annot_ldscore="${annotation_dir}${base%.*}."
    output="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/${base%.*}"
    
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

# visualize
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/visualization.R \
    --ldsc_dir '/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/' \
    --output "/u/home/l/lixinzhe/project-geschwind/plot/$(date +%F)-mdd-heritability-developmental-cell-type-by-time-point.pdf"

###########################################################################################
######                  LDSC Partitioned Heritability SCZ                            ######
###########################################################################################
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc
annotation_dir="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/annot/"
baseline_dir="/u/home/l/lixinzhe/project-geschwind/software/ldsc/"

gwas_dir="/u/home/l/lixinzhe/project-geschwind/data/GWAS/"
gwas="${gwas_dir}Schizophrenia_pardinas2018_munged.sumstats.gz"
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/partitoined-heritability-submitter.sh"
baselines="/u/home/l/lixinzhe/project-geschwind/software/ldsc/baselineLD."

for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/*.hg19.dmr.bed; do
    base=$(basename "$hg19_DMR_bed")
    annot_ldscore="${annotation_dir}${base%.*}."
    output="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/Schizophrenia_pardinas2018/${base%.*}"
    
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

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/visualization.R \
    --ldsc_dir '/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/Schizophrenia_pardinas2018/' \
    --output "/u/home/l/lixinzhe/project-geschwind/plot/$(date +%F)-scz-pardinas2018-heritability-developmental-cell-type-by-time-point.pdf"


###########################################################################################
######                  LDSC Partitioned Heritability ADHD                           ######
###########################################################################################
conda activate ldsc
cd /u/home/l/lixinzhe/project-geschwind/software/ldsc/ldsc
annotation_dir="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/annot/"
baseline_dir="/u/home/l/lixinzhe/project-geschwind/software/ldsc/"

gwas_dir="/u/home/l/lixinzhe/project-geschwind/data/GWAS/"
gwas="${gwas_dir}ADHD_Demontis2018_munged.sumstats.gz"
submission_script="/u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/partitoined-heritability-submitter.sh"
baselines="/u/home/l/lixinzhe/project-geschwind/software/ldsc/baselineLD."

for hg19_DMR_bed in /u/scratch/l/lixinzhe/tmp-file/*.hg19.dmr.bed; do
    base=$(basename "$hg19_DMR_bed")
    annot_ldscore="${annotation_dir}${base%.*}."
    output="/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/ADHD_Demontis2018/${base%.*}"
    
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

# visualize:
Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/ldsc/visualization.R \
    --ldsc_dir '/u/home/l/lixinzhe/project-geschwind/result/ldsc/brain-dev/partitioned_heritability/ADHD_Demontis2018/' \
    --output "/u/home/l/lixinzhe/project-geschwind/plot/$(date +%F)-ADHD_Demontis2018-heritability-developmental-cell-type-by-time-point.pdf"
