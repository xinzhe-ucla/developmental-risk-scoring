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
import seaborn as sns

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

# scz_gwas = pd.read_csv(
#     '/u/home/l/lixinzhe/project-geschwind/data/GWAS/Schizophrenia_pardinas2018',
#     sep = ' ',
#     comment = '#'
#     )
# scz_gwas.columns = ['ID', 'CHROM', "POS", 'A1', 'A2', 'OR', 'SE', 'PVAL', 'DIRECTION']

# load in the fine map results:
scz_finemap = pd.read_csv('/u/home/l/lixinzhe/project-geschwind/port/finemapped-scz/SCZ-hg19-finemapped-snps.csv', sep = ',')

# columns expected: CHR, BP, P
df = scz_gwas.dropna(subset=["CHROM", "POS", "PVAL", "ID"]).copy()[["CHROM", "POS", "PVAL", "ID"]]
df['is_gws'] = df["ID"].isin(scz_finemap.rsid)

# read in the overlap:
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
        "p_value": p_value
        }

results = []

for file_name, dmr in tqdm(dmr_col.items(), desc="DMR files"):
    res = run_bedtools_intersect_count(
        snp_df=df,
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
######                                    Visualization                              ######
###########################################################################################
# visualize the %overlap for each cell type for each locus
plot_df = result_df.copy()
plot_df['time_point'] = [re.sub('_.*', '', f) for f in plot_df.file_name]
plot_df['cell_type'] = [re.sub('.*_Inh-', '', f) for f in plot_df.file_name]
plot_df['percent_dmr_in_gws'] = plot_df['a_gws_in_dmr'] / (plot_df['a_gws_in_dmr'] + plot_df['c_non_gws_in_dmr']) * 100

# specify order:
time_order = ['2T', '3T', '1m', '4-7m', 'adult']
hue_order = sorted(plot_df["cell_type"].dropna().unique())

plt.figure(figsize=(10, 5))
plot_df_plot = plot_df.copy()

plot_df_plot["time_point"] = pd.Categorical(
    plot_df_plot["time_point"],
    categories=time_order,
    ordered=True
)

def format_p(p):
    if pd.isna(p):
        return "p=NA"
    elif p < 0.001:
        return f"p={p:.1e}"
    else:
        return f"p={p:.3f}"

plot_df_plot["p_label"] = plot_df_plot["p_value"].apply(format_p)

plt.figure(figsize=(10, 5))

ax = sns.barplot(
    data=plot_df_plot,
    x="time_point",
    y="odds_ratio",
    hue="cell_type",
    order=time_order,
    hue_order=hue_order,
    errorbar=None
)

# Add p-value labels by matching each container to one hue level
for container, cell_type in zip(ax.containers, hue_order):

    labels = []

    for time_point in time_order:
        tmp = plot_df_plot[
            (plot_df_plot["time_point"] == time_point) &
            (plot_df_plot["cell_type"] == cell_type)
        ]

        if len(tmp) == 1:
            labels.append(tmp["p_label"].iloc[0])
        else:
            labels.append("")

    ax.bar_label(
        container,
        labels=labels,
        padding=3,
        fontsize=8,
        rotation=90
    )

plt.axhline(1, linestyle="--", linewidth=1)
plt.xlabel("Time point")
plt.ylabel("Fisher odds ratio")
plt.title("SCZ finemapped credible SNP enrichment in DMRs")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Cell type", bbox_to_anchor=(1.05, 1), loc="upper left")

ymax = plot_df_plot["odds_ratio"].replace([np.inf, -np.inf], np.nan).max()
plt.ylim(0, ymax * 1.3)

plt.tight_layout()
plt.savefig(f"/u/home/l/lixinzhe/project-geschwind/plot/{today}_OR_finemapped_in_dmr_by_time_celltype.pdf", bbox_inches="tight")
plt.show()

###########################################################################################
######                                    DRD2 region                                ######
###########################################################################################
# for DRD2 region:
drd2_chr = '11'
drd2_start = 113_280_327 - 250000
drd2_end   = 113_346_120 + 250000

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
    dmr_drd2 = dmr.copy()
    dmr_drd2.columns = ['CHROM', 'START', 'END']
    dmr_drd2 = dmr_drd2.loc[
        (dmr_drd2['CHROM'] == 'chr11')
        & (dmr_drd2['START'] >= drd2_start)
        & (dmr_drd2['END'] <= drd2_end)
    ]
    res = run_bedtools_intersect_count(
        snp_df=snp_drd2,
        dmr=dmr_drd2,
        file_name=file_name
    )
    results.append(res)

drd2_locus_result_df = pd.DataFrame(results)

###########################################################################################
######                                    Visualization                              ######
###########################################################################################
# visualize the %overlap for each cell type for each locus
plot_df = drd2_locus_result_df.copy()
plot_df['time_point'] = [re.sub('_.*', '', f) for f in plot_df.file_name]
plot_df['cell_type'] = [re.sub('.*_Inh-', '', f) for f in plot_df.file_name]
plot_df['percent_dmr_in_gws'] = plot_df['a_gws_in_dmr'] / (plot_df['a_gws_in_dmr'] + plot_df['c_non_gws_in_dmr']) * 100
print(plot_df)

# specify order:
time_order = ['2T', '3T', '1m', '4-7m', 'adult']
hue_order = sorted(plot_df["cell_type"].dropna().unique())

plt.figure(figsize=(10, 5))
plot_df_plot = plot_df.copy()

plot_df_plot["time_point"] = pd.Categorical(
    plot_df_plot["time_point"],
    categories=time_order,
    ordered=True
)

def format_p(p):
    if pd.isna(p):
        return "p=NA"
    elif p < 0.001:
        return f"p={p:.1e}"
    else:
        return f"p={p:.3f}"

plot_df_plot["p_label"] = plot_df_plot["p_value"].apply(format_p)

plt.figure(figsize=(10, 5))

ax = sns.barplot(
    data=plot_df_plot,
    x="time_point",
    y="odds_ratio",
    hue="cell_type",
    order=time_order,
    hue_order=hue_order,
    errorbar=None
)

# Add p-value labels by matching each container to one hue level
for container, cell_type in zip(ax.containers, hue_order):

    labels = []

    for time_point in time_order:
        tmp = plot_df_plot[
            (plot_df_plot["time_point"] == time_point) &
            (plot_df_plot["cell_type"] == cell_type)
        ]

        if len(tmp) == 1:
            labels.append(tmp["p_label"].iloc[0])
        else:
            labels.append("")

    ax.bar_label(
        container,
        labels=labels,
        padding=3,
        fontsize=8,
        rotation=90
    )

plt.axhline(1, linestyle="--", linewidth=1)
plt.xlabel("Time point")
plt.ylabel("Fisher odds ratio")
plt.title("SCZ finemapped credible SNP enrichment in DMRs at DRD2 locus")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Cell type", bbox_to_anchor=(1.05, 1), loc="upper left")

ymax = plot_df_plot["odds_ratio"].replace([np.inf, -np.inf], np.nan).max()
plt.ylim(0, ymax * 1.3)

plt.tight_layout()
plt.savefig(f"/u/home/l/lixinzhe/project-geschwind/plot/{today}_OR_finemapped_in_DRD2_locus_dmr_by_time_celltype.pdf", bbox_inches="tight")
plt.show()
