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

gtex_col = {}
for gtex_file in tqdm(gtex_files):
    tissue = re.sub('.*/', '', gtex_file)
    tissue = re.sub('.marginal_significant.txt', '', tissue)
    
    gtex = pd.read_csv(
        gtex_file,
        sep = '\t'
        )
    gtex_col[tissue] = gtex

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


###########################################################################################
######            Check for triple positives at DRD2 locus with Finemap              ######
###########################################################################################
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

# output the DRD locus of finemapped result:
output_bed(
    finemap_drd2,
    chr = 'chromosome',
    position = 'position',
    id = 'rsid',
    output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/PGC3_SCZ_wave3_finemapped_snp_drd2.bed'
    )

results_col = {}
gtex_anno_col = {}
finemap_anno_col = {}
for tissue, gtex in tqdm(gtex_col.items(), desc = 'intersecting GTEX with SCZ finemapped variants'):
    gtex_drd2 = gtex.loc[gtex.stripped_id == 'ENSG00000149295', :].copy()
    
    # output:
    output_bed(
        gtex_drd2,
        chr = 'chromosome',
        position = 'position',
        id = 'variant_id',
        output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_{tissue}_marginal_significant.bed'
    )
    
    # intersect:
    result = bedtool_intersect(
        bed1_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/PGC3_SCZ_wave3_finemapped_snp_drd2.bed",
        bed2_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_{tissue}_marginal_significant.bed',
        out_path = f"/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_{tissue}_SCZ_finemap_intersection.bed"
        )
    
    # short circuit so we don't need to waste time:
    if result is None:
        continue
    
    # for each of the DMR file, also do the same thing:
    across_dmr_age_cell_type = {}
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
        os.makedirs(f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/drd2_intersection/{tissue}/', exist_ok = True)
        result = bedtool_intersect(
            bed1_path = f"/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_{tissue}_SCZ_finemap_intersection.bed",
            bed2_path = output_path,
            out_path = f"/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/drd2_intersection/{tissue}/{file}.intersect.bed"
            )
        
        # write out the cell type and age of the samples:
        if result is not None:
            time_point = re.sub('_.*', '', file)
            cell_type = re.sub('.*MSN-', '', file)
            cell_type = re.sub('.hypo_dmr_overlap.hg19.dmr.bed', '', cell_type)
            result['cell_type'] = cell_type
            result['time_point'] = time_point
            across_dmr_age_cell_type[file] = result
    
    if len(across_dmr_age_cell_type) > 0:
        result_for_tissue = pd.concat(across_dmr_age_cell_type)
        result_for_tissue['tissue'] = tissue
        # collect the results:
        results_col[tissue] = result_for_tissue
        
        # we can also obtain the annotations:
        finemap_result, gtex_tissue_result = get_info(result_for_tissue, annotation1 = finemap_drd2, annotation2 = gtex_drd2, anno1_type = 'finemap')
        gtex_tissue_result['tissue'] = tissue
        gtex_anno_col[tissue] = gtex_tissue_result
        finemap_anno_col[tissue] = finemap_result

# concatenate results:
finemap_anno = pd.concat(finemap_anno_col, ignore_index = True).drop_duplicates()
gtex_anno = pd.concat(gtex_anno_col, ignore_index = True).drop_duplicates()
intersection_result = pd.concat(results_col, ignore_index = True)

# write out the result:
os.makedirs('/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/result/finemap/', exist_ok = True)
finemap_anno.to_csv(f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/result/finemap/{today}_finemap_result_for_triple_positive_loci.csv', index = False)
gtex_anno.to_csv(f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/result/finemap/{today}_gtex_result_for_triple_positive_loci.csv', index = False)
intersection_result.to_csv(f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/result/finemap/{today}_intersection_result_for_triple_positive_loci.csv', index = False)

###########################################################################################
######               Check for triple positives at DRD2 locus with GWAS              ######
###########################################################################################
# subset tto the DRD2 locus for GWAS result:
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


results_col = {}
gtex_anno_col = {}
gwas_anno_col = {}
for tissue, gtex in tqdm(gtex_col.items(), desc = 'intersecting GTEX with SCZ finemapped variants'):
    gtex_drd2 = gtex.loc[gtex.stripped_id == 'ENSG00000149295', :].copy()
    
    # output:
    output_bed(
        gtex_drd2,
        chr = 'chromosome',
        position = 'position',
        id = 'variant_id',
        output_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_{tissue}_marginal_significant.bed'
    )
    
    # intersect:
    result = bedtool_intersect(
        bed1_path = "/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/PGC3_SCZ_wave3_GWS_snp_drd2.bed",
        bed2_path = f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_{tissue}_marginal_significant.bed',
        out_path = f"/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_{tissue}_SCZ_GWAS_intersection.bed"
        )
    
    # short circuit so we don't need to waste time:
    if result is None:
        continue
    
    # for each of the DMR file, also do the same thing:
    across_dmr_age_cell_type = {}
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
        os.makedirs(f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/drd2_intersection/{tissue}/', exist_ok = True)
        result = bedtool_intersect(
            bed1_path = f"/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/gtex_drd2_{tissue}_SCZ_GWAS_intersection.bed",
            bed2_path = output_path,
            out_path = f"/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/drd2_intersection/{tissue}/{file}_GWAS.intersect.bed"
            )
        
        # write out the cell type and age of the samples:
        if result is not None:
            time_point = re.sub('_.*', '', file)
            cell_type = re.sub('.*MSN-', '', file)
            cell_type = re.sub('.hypo_dmr_overlap.hg19.dmr.bed', '', cell_type)
            result['cell_type'] = cell_type
            result['time_point'] = time_point
            across_dmr_age_cell_type[file] = result
    
    if len(across_dmr_age_cell_type) > 0:
        result_for_tissue = pd.concat(across_dmr_age_cell_type)
        result_for_tissue['tissue'] = tissue
        # collect the results:
        results_col[tissue] = result_for_tissue
        
        # we can also obtain the annotations:
        gwas_result, gtex_tissue_result = get_info(result_for_tissue, annotation1 = gwas_drd2, annotation2 = gtex_drd2, anno1_type = 'gwas', anno1_id_col = 'ID')
        gtex_tissue_result['tissue'] = tissue
        gtex_anno_col[tissue] = gtex_tissue_result
        gwas_anno_col[tissue] = gwas_result

# concatenate results:
gwas_anno = pd.concat(gwas_anno_col, ignore_index = True).drop_duplicates()
gtex_anno = pd.concat(gtex_anno_col, ignore_index = True).drop_duplicates()
intersection_result = pd.concat(results_col, ignore_index = True)

# write out the result:
os.makedirs('/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/result/gwas/', exist_ok = True)
gwas_anno.to_csv(f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/result/gwas/{today}_gwas_result_for_triple_positive_loci.csv', index = False)
gtex_anno.to_csv(f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/result/gwas/{today}_gtex_result_for_triple_positive_loci.csv', index = False)
intersection_result.to_csv(f'/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/triple_positive/result/gwas/{today}_intersection_result_for_triple_positive_loci.csv', index = False)
