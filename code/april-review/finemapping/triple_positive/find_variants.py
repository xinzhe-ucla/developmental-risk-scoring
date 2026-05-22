###########################################################################################
######                                    PREAMBLE                                   ######
###########################################################################################
import pandas as pd
from datetime import date
today = date.today()
import re
import os
import subprocess
from tqdm import tqdm

# load in the gwas:
gwas = pd.read_csv(
    '/u/home/l/lixinzhe/project-cluo/data/scz_gwas/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv',
    sep="\t",
    comment = '#'
    )

# load in the finemapping results:
finemap = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/SCZ-hg19-finemapped-snps.csv', sep = ',')

# load in the gtex:
gtex = pd.read_csv(
    '/u/home/l/lixinzhe/project-geschwind/data/GTEX_eQTL/V7/Brain_Frontal_Cortex_BA9.allpairs.txt.gz',
    sep = '\t'
    )

gtex = gtex[gtex.pval_nominal < 0.05].copy()
gtex['chromosome'] = [re.sub('_.*', '',f) for f in gtex.variant_id]
gtex['stripped_id'] = gtex["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)

# split the variant id:
cols = gtex["variant_id"].str.split("_", expand=True)
gtex["chromosome"] = cols[0]
gtex["position"] = cols[1].astype(int)
gtex["ref"] = cols[2]
gtex["alt"] = cols[3]

# load in the DMR:
hypo_dmr_overlap_files = os.listdir('/u/scratch/l/lixinzhe/tmp-file/DMR/')
drd2_hypo_dmr = [f for f in hypo_dmr_overlap_files if f.endswith("hypo_dmr_overlap.hg19.dmr.bed")]
drd2_hypo_dmr = [f for f in drd2_hypo_dmr if 'DRD2' in f]
# drd2_hypo_dmr = [f for f in drd2_hypo_dmr if 'DRD2-BACH2' in f]
# drd2_hypo_dmr = ['2T_Inh-MSN-eMSN.hypo_dmr_overlap.hg19.dmr.bed'] + drd2_hypo_dmr

dmr_col = {}
for file in drd2_hypo_dmr:
    file_name = re.sub('.hypo_dmr_overlap.hg19.dmr.bed', '', file)
    dmr_col[file_name] = pd.read_table(f'/u/scratch/l/lixinzhe/tmp-file/DMR/{file}', sep = '\t', header = None)
    print(file_name)


# define a function that outputs the position:
def output_bed(df, chr, position, id, output_path):
    """
    output summary statistics of gwas or eQTL into a bed file
    
    Parameters
    ----------
    df: pd.DataFrame
        dataframe of the summary statistics, should contain postion and snp id columns
    chr: str
        column name for the chromosome of the SNP
    position: str
        column name for the position of the SNPs
    id: str
        column name for the ID of the SNP
    """
    # subset to the chromosome, position and SNP id:
    work_df = df[[chr, position, id]].copy()
    
    # clean chromosome:
    work_df['chrom'] = work_df[chr].astype(str).str.replace('^chr', '', regex = True)
    
    # get the start and end:
    work_df['start'] = work_df[position] - 1
    work_df['end'] = work_df[position]
    
    # get the id:
    work_df['ID'] = work_df[id]
    
    # get output:
    output_df = work_df[['chrom', 'start', 'end', 'ID']]
    output_df = output_df.dropna().drop_duplicates()
    output_df.to_csv(output_path, sep="\t", header=False, index=False)

def bedtool_intersect(bed1_path, bed2_path, out_path):
    cmd = [
        'bedtools',
        "intersect",
        "-a", bed1_path,
        "-b", bed2_path,
        "-u"
    ]
    with open(out_path, "w") as fout:
        subprocess.run(cmd, stdout=fout, check=True)
    if os.path.getsize(out_path) == 0:
        overlap_idx = []
    else:
        overlap = pd.read_csv(
            out_path,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "id"]
        )
        return overlap
    

###########################################################################################
######                         Check for double positives at DRD2 locus              ######
###########################################################################################
output_bed(
    finemap,
    chr = 'chromosome',
    position = 'position',
    id = 'rsid',
    output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/PGC3_SCZ_wave3_finemapped_snp.bed'
    )

output_bed(
    gtex,
    chr = 'chromosome',
    position = 'position',
    id = 'variant_id',
    output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_frontal_cortex_marginal.bed'
)

drd2_chr = '11'
drd2_start = 113_280_327 - 250000
drd2_end   = 113_346_120 + 250000

# subset to the DRD2 locus for finemap results:
# finemap_drd2 = finemap.loc[
#     (finemap["chromosome"] == drd2_chr)
#     & (finemap["POS"] >= drd2_start)
#     & (finemap["POS"] <= drd2_end)
# ].copy()

finemap_drd2 = finemap.loc[finemap.gene_symbol == 'DRD2', :]

# also subset the gtex:
# gtex_drd2 = gtex.loc[
#     (gtex["chromosome"] == drd2_chr)
#     & (gtex["position"] >= drd2_start)
#     & (gtex["position"] <= drd2_end)
# ].copy()

gtex_drd2 = gtex.loc[gtex.stripped_id == 'ENSG00000149295', :]

# output the two locus and see if there are any intersection:
output_bed(
    finemap_drd2,
    chr = 'chromosome',
    position = 'position',
    id = 'rsid',
    output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/PGC3_SCZ_wave3_finemapped_snp_drd2.bed'
    )

output_bed(
    gtex_drd2,
    chr = 'chromosome',
    position = 'position',
    id = 'variant_id',
    output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_frontal_cortex_marginal.bed'
)

result = bedtool_intersect(
    bed1_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/PGC3_SCZ_wave3_finemapped_snp_drd2.bed",
    bed2_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_frontal_cortex_marginal.bed",
    out_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/drd2_gtex_finemap_variant_intersection.bed"
    )

###########################################################################################
######                   Check for triple positives at DRD2 locus                    ######
###########################################################################################
# next check for DMR positive or not:
for file in drd2_hypo_dmr:
    # load in the data and output it by removing the chr:
    full_path = f'/u/scratch/l/lixinzhe/tmp-file/DMR/{file}'
    dmr = pd.read_table(full_path, sep = '\t', header = None)
    dmr.columns = ['chr', 'start', 'end']
    
    # remove the chr from the chromosome notation:
    dmr['chr'] = [re.sub('chr', '', f) for f in dmr.chr]
    
    # output the data:
    output_path = '/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/tmp_file.bed'
    dmr.to_csv(output_path, sep="\t", header=False, index=False)
    
    # get the result:
    result = bedtool_intersect(
        bed1_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/drd2_gtex_finemap_variant_intersection.bed",
        bed2_path = output_path,
        out_path = f"/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/drd2_intersection/{file}.intersect.bed"
        )
    print(result)

# lets see where this rs77560271 is for each data:
#   chrom      start        end          id
#0     11  113304205  113304206  rs77560271
#4-7m_Inh-MSN-DRD2-BACH2.hypo_dmr_overlap.hg19.dmr.bed

# lets search for the hit in the DMR:
snp_chrom = "11"
snp_start = 113304205
snp_end = 113304206

hit_detail = {}
for file_name, dmr in tqdm(dmr_col.items(), desc="DMR files"):
    dmr_tmp = dmr.copy()
    dmr_tmp.columns = ["chromosome", "start", "end"]

    dmr_tmp["chromosome"] = dmr_tmp["chromosome"].astype(str).str.replace("^chr", "", regex=True)
    dmr_tmp["start"] = dmr_tmp["start"].astype(int)
    dmr_tmp["end"] = dmr_tmp["end"].astype(int)

    hit = dmr_tmp[
        (dmr_tmp["chromosome"] == snp_chrom) &
        (dmr_tmp["start"] < snp_end) &
        (dmr_tmp["end"] > snp_start)
    ]

    if not hit.empty:
        hit_detail[file_name] = hit
        print(file_name)
        print(hit)
    
# ok, lets see where it is in the finemapping and gtex data:
finemap_df = finemap_drd2[finemap_drd2.rsid == 'rs77560271']
gtex_df = gtex_drd2[
    (gtex_drd2.chromosome == snp_chrom) &
    (gtex_drd2.position == snp_end)]

print('gwas hit basd on finemapping:')
print(finemap_df)

print('gtex hit:')
print(gtex_df)

print('dmr hit:')
print(hit_detail)

###########################################################################################
######                       Check for triple positive using GWS SNPS                ######
###########################################################################################
# subset to the DRD2 locus for GWAS results:
gwas_drd2 = gwas.loc[
    (gwas["CHROM"].astype(str) == drd2_chr)
    & (gwas["POS"] >= drd2_start)
    & (gwas["POS"] <= drd2_end)
    & (gwas['PVAL'] <= 5e-6)
].copy()
gwas_drd2.POS = gwas_drd2.POS.astype(int)

output_bed(
    gwas_drd2,
    chr = 'CHROM',
    position = 'POS',
    id = 'ID',
    output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/PGC3_SCZ_wave3_GWS_snp_drd2.bed'
    )

result = bedtool_intersect(
    bed1_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/PGC3_SCZ_wave3_GWS_snp_drd2.bed",
    bed2_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_frontal_cortex_marginal.bed",
    out_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/drd2_gtex_GWS_variant_intersection.bed"
    )

# next check for DMR positive or not:
for file in drd2_hypo_dmr:
    # load in the data and output it by removing the chr:
    full_path = f'/u/scratch/l/lixinzhe/tmp-file/DMR/{file}'
    dmr = pd.read_table(full_path, sep = '\t', header = None)
    dmr.columns = ['chr', 'start', 'end']
    
    # remove the chr from the chromosome notation:
    dmr['chr'] = [re.sub('chr', '', f) for f in dmr.chr]
    
    # output the data:
    output_path = '/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/tmp_file.bed'
    dmr.to_csv(output_path, sep="\t", header=False, index=False)
    
    # get the result:
    result = bedtool_intersect(
        bed1_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/drd2_gtex_GWS_variant_intersection.bed",
        bed2_path = output_path,
        out_path = f"/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/drd2_intersection/tmp.bed"
        )
    print(result)

#    chrom      start        end         id
# 0     11  113404327  113404328  rs2514226

# lets get the annotation:
snp_chrom = "11"
snp_start = 113404327
snp_end = 113404328

hit_detail = {}
for file_name, dmr in tqdm(dmr_col.items(), desc="DMR files"):
    dmr_tmp = dmr.copy()
    dmr_tmp.columns = ["chromosome", "start", "end"]

    dmr_tmp["chromosome"] = dmr_tmp["chromosome"].astype(str).str.replace("^chr", "", regex=True)
    dmr_tmp["start"] = dmr_tmp["start"].astype(int)
    dmr_tmp["end"] = dmr_tmp["end"].astype(int)

    hit = dmr_tmp[
        (dmr_tmp["chromosome"] == snp_chrom) &
        (dmr_tmp["start"] < snp_end) &
        (dmr_tmp["end"] > snp_start)
    ]

    if not hit.empty:
        hit_detail[file_name] = hit
        print(file_name)
        print(hit)
    
# ok, lets see where it is in the finemapping and gtex data:
gwas_df = gwas[gwas.ID == 'rs2514226']
gtex_df = gtex_drd2[
    (gtex_drd2.chromosome == snp_chrom) &
    (gtex_drd2.position == snp_end)]

print('gwas hit basd on pval:')
print(gwas_df)

print('gtex hit:')
print(gtex_df)

print('dmr hit:')
print(hit_detail)