import pandas as pd
import polars as pl
# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
today = date.today()
import os
import tempfile
import subprocess
from scipy.stats import fisher_exact
import re
from tqdm import tqdm
import seaborn as sns
from pathlib import Path
from statsmodels.stats.multitest import multipletests



###########################################################################################
######                             Plot out the DRD2 region                          ######
###########################################################################################
# first read in a SCZ GWAS:
# cd 
# wget https://figshare.com/ndownloader/files/34517861

# load in the schizophrenia gwas:
# scz_gwas = pd.read_csv(
#     '/u/home/l/lixinzhe/project-cluo/data/scz_gwas/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv',
#     sep="\t",
#     comment = '#'
#     )

scz_gwas = pd.read_csv(
    '/u/home/l/lixinzhe/project-geschwind/data/GWAS/Schizophrenia_pardinas2018',
    sep = ' ',
    comment = '#'
    )
scz_gwas.columns = ['ID', 'CHROM', "POS", 'A1', 'A2', 'OR', 'SE', 'PVAL', 'DIRECTION']

# columns expected: CHR, BP, P
df = scz_gwas.dropna(subset=["CHROM", "POS", "PVAL", "ID"]).copy()[["CHROM", "POS", "PVAL", "ID"]]
df['is_gws'] = df["PVAL"] < 5e-8

# read in the overlap:
hypo_dmr_overlap_files = os.listdir('/u/home/h/hex002/project-cluo/BICAN/loop_DMR2/')
drd2_hypo_dmr = [f for f in hypo_dmr_overlap_files if f.endswith("hypo_dmr_overlap.bed")]
drd2_hypo_dmr = [f for f in drd2_hypo_dmr if 'DRD2-BACH2' in f]

dmr_col = {}
for file in drd2_hypo_dmr:
    file_name = re.sub('.hypo_dmr_overlap.bed', '', file)
    dmr_col[file_name] = pd.read_table(f'/u/home/h/hex002/project-cluo/BICAN/loop_DMR2/{file}', sep = '\t', header = None)

# prepare the gwas for liftover:
bed = df[["CHROM", "POS", "ID"]].copy()
bed["Chromosome"] = "chr" + bed["CHROM"].astype(str).str.replace("^chr", "", regex=True)
bed["Start"] = bed["POS"].astype(int) - 1   # BED 0-based start
bed["End"] = bed["POS"].astype(int)         # BED half-open end

chain_file = "/u/home/l/lixinzhe/project-geschwind/software/liftOver/aux_file/hg19ToHg38.over.chain.gz"
bed_out = bed[["Chromosome", "Start", "End", "ID"]]

with tempfile.TemporaryDirectory() as tmpdir:
    tmpdir = Path(tmpdir)

    input_bed = tmpdir / "gwas.hg19.bed"
    output_bed = tmpdir / "gwas.hg38.bed"

    # write temporary BED file
    bed_out.to_csv(
        input_bed,
        sep="\t",
        header=False,
        index=False
    )

    # run CrossMap
    cmd = [
        "CrossMap",
        "bed",
        chain_file,
        str(input_bed),
        str(output_bed)
    ]

    result = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True
    )

    # read lifted BED
    lifted = pd.read_csv(
        output_bed,
        sep="\t",
        header=None
    )

# merge the lifted bed with the original gwas:
lifted.columns = ['Chromosome', 'Start', 'End', 'ID']
lifted["CHROM"] = lifted["Chromosome"].str.replace("chr", "", regex=False)
lifted["POS_hg38"] = lifted["End"].astype(int)  # because original SNP was Start=POS-1, End=POS
df['CHROM'] = df['CHROM'].astype(str)

df_lifted = df.merge(
    lifted[["ID", "CHROM", "POS_hg38"]],
    on=["ID", "CHROM"],
    how="left"
)

###########################################################################################
######                             Interesect between DMR and GWAS                   ######
###########################################################################################
# first clean the lifted file:
df_lifted_clean = df_lifted.drop_duplicates(subset=["ID"], keep="first").copy()
df_lifted_clean = df_lifted_clean.drop_duplicates(subset=["POS_hg38"], keep="first").copy()
df_lifted_clean = df_lifted_clean.dropna().copy()
df_lifted_clean['POS_hg19'] = df_lifted_clean.POS
df_lifted_clean['POS'] = df_lifted_clean.POS_hg38

# define some function to help with the intersection:
def clean_chrom(x):
    return (
        x.astype(str)
         .str.replace("^chr", "", regex=True)
    )

def clean_chrom(x):
    return (
        x.astype(str)
         .str.replace("^chr", "", regex=True)
    )


def run_bedtools_intersect_count(
    snp_df,
    dmr,
    file_name=None,
    bedtools="bedtools"
):
    """
    Use bedtools intersect -u to mark SNPs that overlap at least one DMR.

    Assumptions:
    - snp_df has columns: chrom, POS, is_gws
    - POS is 1-based GWAS position
    - dmr has columns: chrom, start, end
    - dmr start/end are BED-style 0-based half-open coordinates
      If your DMR start/end are 1-based inclusive, see note below.
    """

    snp = snp_df.copy()
    snp["chrom"] = clean_chrom(snp["CHROM"])
    snp["POS"] = snp["POS"].astype(int)

    dmr = dmr.copy()
    dmr.columns = ["chrom", "start", "end"]
    dmr["chrom"] = clean_chrom(dmr["chrom"])
    dmr["start"] = dmr["start"].astype(int)
    dmr["end"] = dmr["end"].astype(int)

    # SNP BED: chrom, start0, end0, original_index
    snp_bed = pd.DataFrame({
        "chrom": snp["chrom"],
        "start": snp["POS"] - 1,
        "end": snp["POS"],
        "idx": snp.index.astype(str)
    })

    dmr_bed = dmr[["chrom", "start", "end"]].dropna()

    # Remove invalid intervals
    dmr_bed = dmr_bed.loc[dmr_bed["end"] > dmr_bed["start"]].copy()

    if snp_bed.empty or dmr_bed.empty:
        in_dmr = pd.Series(False, index=snp.index)
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            snp_path = os.path.join(tmpdir, "snps.bed")
            dmr_path = os.path.join(tmpdir, "dmr.bed")
            out_path = os.path.join(tmpdir, "overlap.bed")

            snp_bed.to_csv(snp_path, sep="\t", header=False, index=False)
            dmr_bed.to_csv(dmr_path, sep="\t", header=False, index=False)

            cmd = [
                bedtools,
                "intersect",
                "-a", snp_path,
                "-b", dmr_path,
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
                    names=["chrom", "start", "end", "idx"]
                )
                overlap_idx = overlap["idx"].astype(snp.index.dtype)

            in_dmr = pd.Series(False, index=snp.index)
            in_dmr.loc[overlap_idx] = True

    is_gws = snp["is_gws"].astype(bool)

    a = int((is_gws & in_dmr).sum())
    b = int((is_gws & ~in_dmr).sum())
    c = int((~is_gws & in_dmr).sum())
    d = int((~is_gws & ~in_dmr).sum())

    table = np.array([[a, b], [c, d]])

    odds_ratio, p_value = fisher_exact(table, alternative="greater")

    n_gws = int(is_gws.sum())
    n_gws_in_dmr = a
    pct_gws_in_dmr = 100 * n_gws_in_dmr / n_gws if n_gws > 0 else np.nan

    return {
        "file_name": file_name,
        "a_gws_in_dmr": a,
        "b_gws_not_in_dmr": b,
        "c_non_gws_in_dmr": c,
        "d_non_gws_not_in_dmr": d,
        "n_gws": n_gws,
        "n_gws_in_dmr": n_gws_in_dmr,
        "pct_gws_in_dmr": pct_gws_in_dmr,
        "odds_ratio": odds_ratio,
        "p_value": p_value
        }

# call the function:
results = []

for file_name, dmr in tqdm(dmr_col.items(), desc="DMR files"):
    res = run_bedtools_intersect_count(
        snp_df=df_lifted_clean,
        dmr=dmr,
        file_name=file_name
    )
    results.append(res)

# get the multiple tested pvalue:
result_df = pd.DataFrame(results)
result_df['p_value_bonf'] = multipletests(result_df.p_value, method = 'bonferroni')[1]
pd.set_option('display.max_columns', None)
print(result_df)

###########################################################################################
######                                    DRD2 region                                ######
###########################################################################################
# for DRD2 region:
drd2_chr = '11'
drd2_start = 113_409_595 - 500000
drd2_end   = 113_475_279 + 500000

# -----------------------------
# 1. Subset SNPs to DRD2 locus
# -----------------------------
snp_drd2 = df_lifted_clean.copy()
snp_drd2["CHROM"] = clean_chrom(snp_drd2["CHROM"])
snp_drd2["POS"] = snp_drd2["POS"].astype(int)
snp_drd2["is_gws"] = snp_drd2["is_gws"].astype(bool)

snp_drd2 = snp_drd2.loc[
    (snp_drd2["CHROM"] == drd2_chr)
    & (snp_drd2["POS"] >= drd2_start)
    & (snp_drd2["POS"] <= drd2_end)
].copy()

print("Number of SNPs in DRD2 region:", len(snp_drd2))
print("Number of GWS SNPs in DRD2 region:", int(snp_drd2["is_gws"].sum()))

results = []

for file_name, dmr in tqdm(dmr_col.items(), desc="DMR files"):
    dmr_drd2 = dmr.copy()
    dmr_drd2.columns = ['CHROM', 'START', 'END']
    dmr_drd2 = dmr_drd2.loc[
        (dmr_drd2['CHROM'] == 'chr11')
        & (dmr_drd2['END'] >= drd2_start)
        & (dmr_drd2['START'] <= drd2_end)
    ]
    res = run_bedtools_intersect_count(
        snp_df=snp_drd2,
        dmr=dmr_drd2,
        file_name=file_name
    )
    results.append(res)

drd2_locus_result_df = pd.DataFrame(results)
print(drd2_locus_result_df)