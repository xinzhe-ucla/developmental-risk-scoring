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
from statsmodels.stats.multitest import multipletests
import re
from tqdm import tqdm

###########################################################################################
######                             Plot out the DRD2 region                          ######
###########################################################################################
# first read in a SCZ GWAS:
# cd 
# wget https://figshare.com/ndownloader/files/34517861

# load in the schizophrenia gwas:
scz_gwas = pd.read_csv(
    '/u/home/l/lixinzhe/project-cluo/data/scz_gwas/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv',
    sep="\t",
    comment = '#'
    )
scz_finemap = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/SCZ-hg19-finemapped-snps.csv', sep = ',')

# columns expected: CHR, BP, P
df = scz_gwas.dropna(subset=["CHROM", "POS", "PVAL", "ID"]).copy()[["CHROM", "POS", "PVAL", "ID"]]
df['is_gws'] = df["ID"].isin(scz_finemap.rsid)

# read in the overlap:
hypo_dmr_overlap_files = os.listdir('/u/scratch/l/lixinzhe/tmp-file/DMR/')
drd2_hypo_dmr = [f for f in hypo_dmr_overlap_files if f.endswith("hypo_dmr_overlap.hg19.dmr.bed")]
drd2_hypo_dmr = [f for f in drd2_hypo_dmr if 'DRD2' in f]
# drd2_hypo_dmr = [f for f in drd2_hypo_dmr if 'DRD2-BACH2' in f]
drd2_hypo_dmr = ['2T_Inh-MSN-eMSN.hypo_dmr_overlap.hg19.dmr.bed'] + drd2_hypo_dmr

dmr_col = {}
for file in drd2_hypo_dmr:
    file_name = re.sub('.hypo_dmr_overlap.hg19.dmr.bed', '', file)
    dmr_col[file_name] = pd.read_table(f'/u/scratch/l/lixinzhe/tmp-file/DMR/{file}', sep = '\t', header = None)
    print(file_name)


###########################################################################################
######                                Define function for intersection               ######
###########################################################################################

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
        "p_value": p_value,
        "p_value_bonf": multipletests(p_value, method="bonferroni")[1],
    }

results = []

for file_name, dmr in tqdm(dmr_col.items(), desc="DMR files"):
    res = run_bedtools_intersect_count(
        snp_df=df,
        dmr=dmr,
        file_name=file_name
    )
    results.append(res)

pd.set_option('display.max_columns', None)
result_df = pd.DataFrame(results)
print(result_df)

###########################################################################################
######                                    DRD2 region                                ######
###########################################################################################
# for DRD2 region:
drd2_chr = '11'
drd2_start = 113_280_337 - 250000
drd2_end   = 113_346_413 + 250000

# -----------------------------
# 1. Subset SNPs to DRD2 locus
# -----------------------------
snp_drd2 = df.copy()
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
    res = run_bedtools_intersect_count(
        snp_df=snp_drd2,
        dmr=dmr,
        file_name=file_name
    )
    results.append(res)

result_df = pd.DataFrame(results)
print(result_df)