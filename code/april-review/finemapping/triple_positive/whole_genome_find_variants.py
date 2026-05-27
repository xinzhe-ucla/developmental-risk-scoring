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

# define a function that outputs the position:
def output_bed(df, chr, position, id, output_path, p_value):
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
    p_value: str
        column name for te pvalue of the SNP
    """
    # subset to the chromosome, position and SNP id:
    work_df = df[[chr, position, id, p_value]].copy()
    
    # clean chromosome:
    work_df['chrom'] = work_df[chr].astype(str).str.replace('^chr', '', regex = True)
    
    # get the start and end:
    work_df['start'] = work_df[position] - 1
    work_df['end'] = work_df[position]
    
    # get the id:
    work_df['ID'] = work_df[id]
    work_df['Pval'] = work_df[p_value]
    
    # get output:
    output_df = work_df[['chrom', 'start', 'end', 'ID', 'Pval']]
    output_df = output_df.dropna().drop_duplicates()
    output_df.to_csv(output_path, sep="\t", header=False, index=False)

def bedtool_intersect(bed1_path, bed2_path, out_path):
    cmd = [
        "bedtools",
        "intersect",
        "-a", bed1_path,
        "-b", bed2_path,
        "-wa",
        "-wb"
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
    
def get_info(result_df, annotation1, annotation2, anno1_type = 'finemap', anno1_id_col = 'rsid'):
    """
    get the information using the annotations
    
    parameters
    ----------
    result_df: pd.DataFrame
        Contains chrom, start, end, id
    annotation1: pd.DataFrame
        contains the finemap annotation or the gwas annotation
    annotation2: pd.DataFrame
        contains the gtex annotation
        assume there is chromosome and position column
    anno1_type: str
        what is the type of the annotation 1, can be finemap or gwas
    anno1_id_col: str
        the column index for the annotation 1 for ID column
    """
    # get the start end chromosome location
    anno1_result_col = {}
    anno2_result_col = {}
    
    # formating:
    annotation2 = annotation2.copy()
    annotation2['chromosome'] = annotation2['chromosome'].astype(str)
    annotation2['position'] = annotation2['position'].astype(int)
    for rsid in result_df.id.unique():
        loci_result = result_df[result_df.id == rsid].copy()
        loci_chr = str(loci_result.chrom.unique()[0])
        loci_start = int(loci_result["start"].iloc[0])
        loci_end = int(loci_result["end"].iloc[0])
        
        # get the finemap result:
        if anno1_type == 'finemap':
            # subset to that snp
            finemap_df = annotation1[annotation1[anno1_id_col] == rsid]
            anno1_result_col[rsid] = finemap_df
        
        elif anno1_type == 'gwas':
            # subset to that snp;
            gwas_df = annotation1[annotation1[anno1_id_col] == rsid]
            anno1_result_col[rsid] = gwas_df
        
        else:
            raise ValueError('anno1_type must be finemap or gwas')
        
        # next get the gtex result:
        gtex_df = annotation2[
            (annotation2.chromosome == loci_chr) &
            (annotation2.position == loci_end)]
        anno2_result_col[rsid] = gtex_df
    
    if len(anno1_result_col) > 0:
        anno1_result = pd.concat(anno1_result_col)
    else:
        anno1_result = None
    if len(anno2_result_col) > 0:
        anno2_result = pd.concat(anno2_result_col)
    else:
        anno2_result = None
    
    return anno1_result, anno2_result
    
###########################################################################################
######                                 Load in the data                              ######
###########################################################################################
# load in the gwas:
gwas = pd.read_csv(
    '/u/home/l/lixinzhe/project-cluo/data/scz_gwas/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv',
    sep="\t",
    comment = '#'
    )

# load in the finemapping results:
finemap = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/SCZ-hg19-finemapped-snps.csv', sep = ',')

# load in the gtex:
gtex_files = [
    "/u/home/l/lixinzhe/project-geschwind/data/GTEX_eQTL/V7_sig_only/Brain_Caudate_basal_ganglia.marginal_significant.txt",
    "/u/home/l/lixinzhe/project-geschwind/data/GTEX_eQTL/V7_sig_only/Brain_Nucleus_accumbens_basal_ganglia.marginal_significant.txt",
    "/u/home/l/lixinzhe/project-geschwind/data/GTEX_eQTL/V7_sig_only/Brain_Putamen_basal_ganglia.marginal_significant.txt",
    '/u/home/l/lixinzhe/project-geschwind/data/GTEX_eQTL/V7_sig_only/Brain_Frontal_Cortex_BA9.marginal_significant.txt'
]

# output the gwas bed:
gwas_subset = gwas[gwas.PVAL < 5e-6]
gwas_subset['POS'] = gwas_subset['POS'].astype(int)

output_bed(
    gwas_subset,
    chr = 'CHROM',
    position = 'POS',
    id = 'ID',
    p_value = 'PVAL',
    output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/PGC3_SCZ_wave3_significant_snp.bed'
    )

output_bed(
    finemap,
    chr = 'chromosome',
    position = 'position',
    id = 'rsid',
    p_value = 'pval',
    output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/PGC3_SCZ_wave3_finemapped_snp.bed'
    )

# intersection between gwas and gtex:
gtex_col = {}
os.makedirs('/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/double_positive/', exist_ok = True)

for gtex_file in tqdm(gtex_files):
    tissue = re.sub('.*/', '', gtex_file)
    tissue = re.sub('.marginal_significant.txt', '', tissue)
    
    gtex = pd.read_csv(
        gtex_file,
        sep = '\t'
        )
    gtex_col[tissue] = gtex
    
    # output:
    output_bed(
        gtex,
        chr = 'chromosome',
        position = 'position',
        id = 'variant_id',
        p_value = 'pval_nominal',
        output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/gtex_{tissue}_marginal_significant.bed'
    )
    
    # intersection:
    finemap_double_positive_result = bedtool_intersect(
        bed1_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/PGC3_SCZ_wave3_finemapped_snp.bed",
        bed2_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/gtex_{tissue}_marginal_significant.bed',
        out_path = f"/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/double_positive/gtex_{tissue}_SCZ_finemap_intersection.bed"
        )
    
    # intersection;
    gwas_double_positive_result = bedtool_intersect(
        bed1_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/PGC3_SCZ_wave3_significant_snp.bed",
        bed2_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/gtex_{tissue}_marginal_significant.bed',
        out_path = f"/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/double_positive/gtex_{tissue}_SCZ_gwas_intersection.bed"
        )
    