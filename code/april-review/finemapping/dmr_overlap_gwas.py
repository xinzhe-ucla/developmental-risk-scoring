import pandas as pd
import polars as pl
# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
today = date.today()
import os
import re

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

###########################################################################################
######                              look at near the DRD2 region                     ######
###########################################################################################


# columns expected: CHR, BP, P
df = scz_gwas.dropna(subset=["CHROM", "POS", "PVAL", "ID"]).copy()[["CHROM", "POS", "PVAL", "ID"]]

def plot_locus_manhattan(
    df,
    chrom,
    output_path,
    start=None,
    end=None,
    center=None,
    window=None,
    chrom_col="CHROM",
    pos_col="POS",
    p_col="PVAL",
    id_col="ID",
    annotate_top=False,
    top_n=5,
    figsize=(10, 5),
):
    """
    Plot a local Manhattan plot for one genomic region.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain chromosome, position, and p-value columns.
    chrom : str or int
        Chromosome to plot.
    start, end : int, optional
        Region boundaries. Use these directly, or use center + window.
    center : int, optional
        Center position of region.
    window : int, optional
        Half-window around center. Example: center=113300000, window=500000
        plots [112800000, 113800000].
    chrom_col, pos_col, p_col, id_col : str
        Column names in df.
    annotate_top : bool
        Whether to label the most significant variants.
    top_n : int
        Number of top variants to annotate if annotate_top=True.
    figsize : tuple
        Figure size.

    Returns
    -------
    pandas.DataFrame
        The subset of data used for plotting.
    """

    # Make a working copy
    plot_df = df[[chrom_col, pos_col, p_col] + ([id_col] if id_col in df.columns else [])].copy()

    # Clean types
    plot_df[chrom_col] = plot_df[chrom_col].astype(str)
    chrom = str(chrom)

    plot_df[pos_col] = pd.to_numeric(plot_df[pos_col], errors="coerce")
    plot_df[p_col] = pd.to_numeric(plot_df[p_col], errors="coerce")

    # Drop bad rows
    plot_df = plot_df.dropna(subset=[chrom_col, pos_col, p_col])
    plot_df = plot_df[(plot_df[p_col] > 0) & (plot_df[p_col] <= 1)]

    # Define region
    if center is not None and window is not None:
        start = center - window
        end = center + window

    if start is None or end is None:
        raise ValueError("Provide either (start and end) or (center and window).")

    # Subset region
    plot_df = plot_df[
        (plot_df[chrom_col] == chrom) &
        (plot_df[pos_col] >= start) &
        (plot_df[pos_col] <= end)
    ].copy()

    if plot_df.empty:
        raise ValueError("No variants found in the requested region.")

    # Compute -log10(p)
    plot_df["minus_log10_p"] = -np.log10(plot_df[p_col])

    # Sort by position
    plot_df = plot_df.sort_values(pos_col)

    # Plot
    plt.figure(figsize=figsize)
    plt.scatter(plot_df[pos_col], plot_df["minus_log10_p"], s=8)
    plt.axhline(-np.log10(5e-8), linestyle="--", linewidth=1)
    plt.xlabel(f"Position on chr{chrom}")
    plt.ylabel("-log10(P)")
    plt.title(f"Local Manhattan plot: chr{chrom}:{start:,}-{end:,}")
    plt.tight_layout()

    # Optional annotation of top hits
    if annotate_top and id_col in plot_df.columns:
        top_hits = plot_df.nsmallest(top_n, p_col)
        for _, row in top_hits.iterrows():
            plt.annotate(
                row[id_col],
                (row[pos_col], row["minus_log10_p"]),
                fontsize=8,
            )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    return plot_df
    
plot_locus_manhattan(df, chrom=11, start=113475398 - 500000, end=113475398+500000, output_path=f"/u/home/l/lixinzhe/project-geschwind/plot/{today}-drd2_locus.png")

# read in the overlap:
hypo_dmr_overlap_files = os.listdir('/u/home/h/hex002/project-cluo/BICAN/loop_DMR2/')
drd2_hypo_dmr = [f for f in hypo_dmr_overlap_files if f.endswith("hypo_dmr_overlap.bed")]
drd2_hypo_dmr = [f for f in drd2_hypo_dmr if 'DRD2-BACH2' in f]

dmr_col = {}
for file in drd2_hypo_dmr:
    file_name = re.sub('.hypo_dmr_overlap.bed', '', file)
    dmr_col[file_name] = pd.read_table(f'/u/home/h/hex002/project-cluo/BICAN/loop_DMR2/{file}', sep = '\t', header = None)

###########################################################################################
######                                    overlap                                    ######
###########################################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from adjustText import adjust_text



def plot_locus_manhattan_with_dmr(
    gwas_df,
    dmr_df,
    chrom,
    output_path,
    start=None,
    end=None,
    center=None,
    window=None,
    chrom_col="CHROM",
    pos_col="POS",
    p_col="PVAL",
    id_col="ID",
    annotate_top=False,
    top_n=5,
    dpi=300,
    figsize=(10, 5),
    dmr_is_bed=True,
    shade_dmr=True,
):
    """
    Plot a regional Manhattan plot, overlay DMR intervals, and color SNPs red if they overlap a DMR.

    Parameters
    ----------
    gwas_df : pandas.DataFrame
        GWAS dataframe with columns like CHROM, POS, PVAL, ID.
    dmr_df : pandas.DataFrame
        DMR dataframe with no header:
        column 0 = chromosome, column 1 = start, column 2 = end.
    chrom : str or int
        Chromosome to plot.
    output_path : str
        File path to save the plot.
    start, end : int, optional
        Region boundaries.
    center, window : int, optional
        Alternative way to define the region.
    dmr_is_bed : bool
        If True, interpret DMR coordinates as BED (0-based, half-open) and convert to 1-based.
    shade_dmr : bool
        If True, shade DMR intervals in the background.
    """

    chrom = str(chrom)

    # Clean GWAS dataframe
    plot_df = gwas_df.copy()
    plot_df[chrom_col] = plot_df[chrom_col].astype(str).str.replace("^chr", "", regex=True)
    plot_df[pos_col] = pd.to_numeric(plot_df[pos_col], errors="coerce")
    plot_df[p_col] = pd.to_numeric(plot_df[p_col], errors="coerce")

    plot_df = plot_df.dropna(subset=[chrom_col, pos_col, p_col])
    plot_df = plot_df[(plot_df[p_col] > 0) & (plot_df[p_col] <= 1)]

    # Define region
    if center is not None and window is not None:
        start = center - window
        end = center + window

    if start is None or end is None:
        raise ValueError("Provide either (start and end) or (center and window).")

    # Subset GWAS region
    plot_df = plot_df[
        (plot_df[chrom_col] == chrom) &
        (plot_df[pos_col] >= start) &
        (plot_df[pos_col] <= end)
    ].copy()

    if plot_df.empty:
        raise ValueError("No GWAS variants found in this region.")

    plot_df["minus_log10_p"] = -np.log10(plot_df[p_col])
    plot_df = plot_df.sort_values(pos_col)

    # Clean DMR dataframe
    dmr_sub = dmr_df.iloc[:, :3].copy()
    dmr_sub.columns = ["chrom", "start", "end"]
    dmr_sub["chrom"] = dmr_sub["chrom"].astype(str).str.replace("^chr", "", regex=True)
    dmr_sub["start"] = pd.to_numeric(dmr_sub["start"], errors="coerce")
    dmr_sub["end"] = pd.to_numeric(dmr_sub["end"], errors="coerce")
    dmr_sub = dmr_sub.dropna(subset=["chrom", "start", "end"])

    # BED (0-based, half-open) -> 1-based inclusive
    if dmr_is_bed:
        dmr_sub["start"] = dmr_sub["start"] + 1
        dmr_sub['end'] = dmr_sub['end'] + 1

    # Keep only DMRs overlapping the plotting window
    dmr_sub = dmr_sub[
        (dmr_sub["chrom"] == chrom) &
        (dmr_sub["end"] >= start) &
        (dmr_sub["start"] <= end)
    ].copy()

    # Mark SNPs that overlap any DMR
    plot_df["overlap_dmr"] = False

    if not dmr_sub.empty:
        for _, row in dmr_sub.iterrows():
            mask = (plot_df[pos_col] >= row["start"]) & (plot_df[pos_col] <= row["end"])
            plot_df.loc[mask, "overlap_dmr"] = True

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    # Optional DMR shading
    if shade_dmr and not dmr_sub.empty:
        for _, row in dmr_sub.iterrows():
            ax.axvspan(row["start"], row["end"], alpha=0.15)

    # Non-overlapping SNPs
    non_overlap = plot_df[~plot_df["overlap_dmr"]]
    ax.scatter(non_overlap[pos_col], non_overlap["minus_log10_p"], s=8)

    # Overlapping SNPs in red
    overlap = plot_df[plot_df["overlap_dmr"]]
    ax.scatter(overlap[pos_col], overlap["minus_log10_p"], s=12, color="red")

    # Genome-wide significance line
    ax.axhline(-np.log10(5e-8), linestyle="--", linewidth=1)
    ax.xaxis.set_major_locator(MultipleLocator(10000))
    ax.tick_params(axis="x", rotation=90)

    # Optional annotation of top hits
    if annotate_top and id_col in plot_df.columns:
        top_hits = plot_df.loc[
            (plot_df['overlap_dmr']) & (plot_df[p_col] < 5e-8),]
        texts = []
        for _, row in top_hits.iterrows():
            texts.append(
                ax.text(
                    row[pos_col],
                    row["minus_log10_p"],
                    row[id_col],
                    fontsize=8
                )
            )
        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", lw=0.5)
        )
    
    # add all the labels:
    ax.set_xlabel(f"Position on chr{chrom}")
    ax.set_ylabel("-log10(P)")
    ax.set_title(f"chr{chrom}:{start:,}-{end:,}")

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

    return plot_df, dmr_sub


for file_name in dmr_col.keys():
    plot_locus_manhattan_with_dmr(
        df,
        dmr_df = dmr_col[file_name],
        chrom=11,
        center=113475398,
        window=250000,
        dmr_is_bed = True,
        annotate_top = True,
        output_path=f'/u/home/l/lixinzhe/project-geschwind/plot/{today}-{file_name}-drd2_locus-zero.png'
    )

###########################################################################################
######                        Count the number of significant hits in DRD            ######
###########################################################################################
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

# -----------------------------
# 1. Define DRD2 locus
# -----------------------------
drd2_chr = "11"
drd2_start = 113475398 - 250000
drd2_end   = 113475398 + 250000

# -----------------------------
# 2. Subset SNPs to DRD2 locus
# -----------------------------
snp_df = df.loc[:, ["CHROM", "POS", "PVAL", "ID"]].copy()
snp_df["CHROM"] = snp_df["CHROM"].astype(str).str.replace("^chr", "", regex=True)
snp_df["POS"] = snp_df["POS"].astype(int)

snp_df = snp_df.loc[
    (snp_df["CHROM"] == drd2_chr) &
    (snp_df["POS"] >= drd2_start) &
    (snp_df["POS"] <= drd2_end)
].copy()

snp_df["is_gws"] = snp_df["PVAL"] < 5e-8

# -----------------------------
# 3. Subset DMRs to DRD2 locus
# -----------------------------
for file_name in dmr_col.keys():
    dmr = dmr_col[file_name].copy()
    dmr.columns = ["chrom", "start", "end"]
    dmr["chrom"] = dmr["chrom"].astype(str).str.replace("^chr", "", regex=True)
    dmr["start"] = dmr["start"].astype(int)
    dmr["end"] = dmr["end"].astype(int)

    dmr = dmr.loc[
        (dmr["chrom"] == drd2_chr) &
        (dmr["end"] >= drd2_start) &
        (dmr["start"] <= drd2_end)
    ].copy()

    # -----------------------------
    # 4. Mark SNPs that fall in any DMR
    # -----------------------------
    snp_df["in_dmr"] = False

    if not dmr.empty and not snp_df.empty:
        hit = pd.Series(False, index=snp_df.index)
        for start, end in zip(dmr["start"], dmr["end"]):
            hit |= snp_df["POS"].between(start, end)
        snp_df["in_dmr"] = hit

    # -----------------------------
    # 5. Percentage of GWS SNPs in DMR within DRD2
    # -----------------------------
    gws = snp_df.loc[snp_df["is_gws"]].copy()

    n_gws = len(gws)
    n_gws_in_dmr = int(gws["in_dmr"].sum())
    pct_gws_in_dmr = 100 * n_gws_in_dmr / n_gws if n_gws > 0 else np.nan

    # print(f"DRD2 locus SNPs: {len(snp_df)}")
    # print(f"DRD2 locus DMRs: {len(dmr)}")
    # print(f"GWS SNPs in DRD2 locus: {n_gws}")
    # print(f"GWS SNPs in DMRs within DRD2 locus: {n_gws_in_dmr}")
    # print(f"Percentage: {pct_gws_in_dmr:.2f}%")

    # -----------------------------
    # 6. Fisher exact test within DRD2 locus
    # -----------------------------
    a = ((snp_df["is_gws"]) & (snp_df["in_dmr"])).sum()
    b = ((snp_df["is_gws"]) & (~snp_df["in_dmr"])).sum()
    c = ((~snp_df["is_gws"]) & (snp_df["in_dmr"])).sum()
    d = ((~snp_df["is_gws"]) & (~snp_df["in_dmr"])).sum()

    table = np.array([[a, b],
                    [c, d]])

    # print(pd.DataFrame(table,
    #                 index=["GWS", "Non-GWS"],
    #                 columns=["In_DMR", "Not_in_DMR"]))

    if table.min() >= 0:
        odds_ratio, p_value = fisher_exact(table, alternative="greater")
        print(f"{file_name}")
        print("percent:", f"{pct_gws_in_dmr:.2f}%")
        print("One-sided Fisher p-value:", p_value * 5)
    