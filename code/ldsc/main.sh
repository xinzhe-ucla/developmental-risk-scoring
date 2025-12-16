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
