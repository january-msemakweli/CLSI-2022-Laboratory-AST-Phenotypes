"""Manuscript tables and two figures from the paired AST parquet."""
from __future__ import annotations

import json
import os
import shutil
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

ENTEROBACTERALES = {
    "ESCHERICHIA COLI",
    "KLEBSIELLA PNEUMONIAE",
    "KLEBSIELLA OXYTOCA",
    "PROTEUS MIRABILIS",
    "ENTEROBACTER CLOACAE COMPLEX",
    "SERRATIA MARCESCENS",
    "CITROBACTER KOSERI",
    "CITROBACTER FREUNDII COMPLEX",
    "MORGANELLA MORGANII",
    "PROVIDENCIA STUARTII",
    "ENTEROBACTER AEROGENES",
    "KLEBSIELLA AEROGENES",
}

PAIR_FIGURE = [
    ("ESCHERICHIA COLI", "piperacillin_tazobactam", r"E. coli, piperacillin-tazobactam"),
    ("ESCHERICHIA COLI", "cefepime", r"E. coli, cefepime"),
    ("ESCHERICHIA COLI", "ciprofloxacin", r"E. coli, ciprofloxacin"),
    ("ESCHERICHIA COLI", "levofloxacin", r"E. coli, levofloxacin"),
    ("ESCHERICHIA COLI", "ceftriaxone", r"E. coli, ceftriaxone"),
    ("ESCHERICHIA COLI", "meropenem", r"E. coli, meropenem"),
    ("KLEBSIELLA PNEUMONIAE", "piperacillin_tazobactam", r"K. pneumoniae, piperacillin-tazobactam"),
    ("KLEBSIELLA PNEUMONIAE", "cefepime", r"K. pneumoniae, cefepime"),
    ("KLEBSIELLA PNEUMONIAE", "ciprofloxacin", r"K. pneumoniae, ciprofloxacin"),
    ("KLEBSIELLA PNEUMONIAE", "levofloxacin", r"K. pneumoniae, levofloxacin"),
    ("PSEUDOMONAS AERUGINOSA", "piperacillin_tazobactam", r"P. aeruginosa, piperacillin-tazobactam"),
    ("STAPHYLOCOCCUS AUREUS", "cefoxitin", r"S. aureus, cefoxitin"),
    ("STAPHYLOCOCCUS AUREUS", "clindamycin", r"S. aureus, clindamycin"),
    ("STAPHYLOCOCCUS AUREUS", "vancomycin", r"S. aureus, vancomycin"),
]

PAIR_TABLE = [
    ("ESCHERICHIA COLI", "piperacillin_tazobactam"),
    ("ESCHERICHIA COLI", "cefepime"),
    ("ESCHERICHIA COLI", "ceftriaxone"),
    ("ESCHERICHIA COLI", "cefazolin"),
    ("ESCHERICHIA COLI", "meropenem"),
    ("ESCHERICHIA COLI", "ciprofloxacin"),
    ("ESCHERICHIA COLI", "levofloxacin"),
    ("ESCHERICHIA COLI", "trimethoprim_sulfamethoxazole"),
    ("KLEBSIELLA PNEUMONIAE", "piperacillin_tazobactam"),
    ("KLEBSIELLA PNEUMONIAE", "cefepime"),
    ("KLEBSIELLA PNEUMONIAE", "ceftriaxone"),
    ("KLEBSIELLA PNEUMONIAE", "ciprofloxacin"),
    ("KLEBSIELLA PNEUMONIAE", "levofloxacin"),
    ("PSEUDOMONAS AERUGINOSA", "piperacillin_tazobactam"),
    ("PSEUDOMONAS AERUGINOSA", "ciprofloxacin"),
    ("STAPHYLOCOCCUS AUREUS", "oxacillin"),
    ("STAPHYLOCOCCUS AUREUS", "cefoxitin"),
    ("STAPHYLOCOCCUS AUREUS", "clindamycin"),
    ("STAPHYLOCOCCUS AUREUS", "vancomycin"),
]

CATS = ["S", "SDD", "I", "R", "NS"]
DOUBLING = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0}


def _serif():
    names = {f.name for f in font_manager.fontManager.ttflist}
    return "Times New Roman" if "Times New Roman" in names else "DejaVu Serif"


def apply_journal_style():
    serif = _serif()
    mpl.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "font.family": serif,
            "font.serif": [serif, "Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.titleweight": "normal",
            "axes.labelsize": 9,
            "axes.labelcolor": "black",
            "axes.edgecolor": "black",
            "axes.linewidth": 0.8,
            "axes.facecolor": "white",
            "axes.grid": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.color": "black",
            "ytick.color": "black",
            "text.color": "black",
            "legend.frameon": False,
            "legend.fontsize": 8,
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def pct(k, n, digits=1):
    if n == 0:
        return None
    return round(100.0 * k / n, digits)


def sir_block(s: pd.Series) -> dict:
    n = int(len(s))
    counts = {c: int((s == c).sum()) for c in CATS}
    out = {"n": n}
    for c in CATS:
        out[c] = counts[c]
        out[f"{c}_pct"] = pct(counts[c], n, 1)
    return out


def pair_sir(df: pd.DataFrame) -> dict:
    n = len(df)
    lab = sir_block(df["lab"])
    clsi = sir_block(df["clsi"])
    s_to_sdd = int(((df["lab"] == "S") & (df["clsi"] == "SDD")).sum())
    s_to_r = int(((df["lab"] == "S") & (df["clsi"] == "R")).sum())
    r_to_s = int(((df["lab"] == "R") & (df["clsi"] == "S")).sum())
    return {
        "n": n,
        "lab": lab,
        "clsi": clsi,
        "s_to_sdd": s_to_sdd,
        "s_to_sdd_pct": pct(s_to_sdd, n, 1),
        "s_to_r": s_to_r,
        "s_to_r_pct": pct(s_to_r, n, 1),
        "r_to_s": r_to_s,
        "r_to_s_pct": pct(r_to_s, n, 1),
        "lab_s_pct": lab["S_pct"],
        "clsi_s_pct": clsi["S_pct"],
        "delta_s": round(lab["S_pct"] - clsi["S_pct"], 1) if lab["S_pct"] is not None else None,
    }


def cluster_ci(df: pd.DataFrame, B: int = 1000, seed: int = 1) -> dict:
    """Patient-clustered bootstrap 95% CIs for lab %S, CLSI %S, and CLSI %SDD."""
    g = (
        df.groupby("anon_id", sort=False)
        .agg(
            n=("clsi", "size"),
            sdd=("clsi", lambda s: int((s == "SDD").sum())),
            lab_s=("lab", lambda s: int((s == "S").sum())),
            clsi_s=("clsi", lambda s: int((s == "S").sum())),
        )
        .to_numpy(dtype=np.float64)
    )
    n_pat = int(g.shape[0])
    tot = float(g[:, 0].sum())
    point = np.array(
        [
            100.0 * g[:, 2].sum() / tot,
            100.0 * g[:, 3].sum() / tot,
            100.0 * g[:, 1].sum() / tot,
        ]
    )
    rng = np.random.default_rng(seed)
    boots = np.empty((B, 3))
    for i in range(B):
        take = g[rng.integers(0, n_pat, size=n_pat)]
        t = take[:, 0].sum()
        boots[i, 0] = 100.0 * take[:, 2].sum() / t
        boots[i, 1] = 100.0 * take[:, 3].sum() / t
        boots[i, 2] = 100.0 * take[:, 1].sum() / t
    lo = np.percentile(boots, 2.5, axis=0)
    hi = np.percentile(boots, 97.5, axis=0)

    def pack(j):
        return {
            "pct": round(float(point[j]), 1),
            "lo": round(float(lo[j]), 1),
            "hi": round(float(hi[j]), 1),
        }

    return {
        "n_patients": n_pat,
        "n_rows": int(tot),
        "B": B,
        "lab_S": pack(0),
        "clsi_S": pack(1),
        "clsi_SDD": pack(2),
    }


def mic_notation(df: pd.DataFrame) -> dict:
    ineq = df["AST_inequality"].astype(str).str.strip()
    val = pd.to_numeric(df["AST_val1"], errors="coerce")
    return {
        "n": int(len(df)),
        "n_le4": int(((val <= 4) & (ineq == "<=")).sum()),
        "n_8_le": int(((val == 8) & (ineq == "<=")).sum()),
        "n_8_exact": int(((val == 8) & (ineq == "")).sum()),
        "n_16_le": int(((val == 16) & (ineq == "<=")).sum()),
        "n_16_exact": int(((val == 16) & (ineq == "")).sum()),
        "n_32_exact": int(((val == 32) & (ineq == "")).sum()),
        "n_64_exact": int(((val == 64) & (ineq == "")).sum()),
        "n_gt64": int(((val == 64) & (ineq == ">")).sum()),
        "n_ge128": int(((val == 128) & (ineq == ">=")).sum()),
    }


def mic_bin(val, ineq, panel):
    """Separate an exact doubling dilution from a censored bound.

    A result stored as <=8 or <=16 is not an MIC of 8 or 16. The bound
    stays in its own bin so the histogram does not treat it as one.
    """
    if str(panel).strip() != "broth_microdilution":
        return None
    try:
        x = float(val)
    except (TypeError, ValueError):
        return None
    if x not in DOUBLING:
        return None
    ineq = str(ineq).strip()
    if x <= 4:
        return "leq4"
    if x == 8 and ineq == "<=":
        return "le8"
    if x == 8 and ineq == "":
        return "exact8"
    if x == 16 and ineq in {"<=", "<"}:
        return "le16"
    if x == 16 and ineq == "":
        return "exact16"
    if x == 32 and ineq == "":
        return "exact32"
    if x == 64 and ineq == "":
        return "exact64"
    if (x == 64 and ineq == ">") or (x >= 128 and ineq == ">="):
        return "hi"
    return None


def save_fig(fig, stem):
    png = os.path.join(C.FIGURES_DIR, f"{stem}.png")
    pdf = os.path.join(C.FIGURES_DIR, f"{stem}.pdf")
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    print("  saved", stem, flush=True)
    return pdf


def fig_susceptibility(analytic: pd.DataFrame):
    apply_journal_style()
    labels = []
    lab_s = []
    clsi_s = []
    for org, drug, lab in PAIR_FIGURE:
        g = analytic.loc[(analytic["organism"] == org) & (analytic["antibiotic"] == drug)]
        labels.append(lab)
        lab_s.append(100.0 * (g["lab"] == "S").mean())
        clsi_s.append(100.0 * (g["clsi"] == "S").mean())

    y = np.arange(len(labels))
    h = 0.38
    fig, ax = plt.subplots(figsize=(7.2, 6.6))
    ax.barh(y + h / 2, lab_s, height=h, color="#4D4D4D", label="Laboratory-reported")
    ax.barh(y - h / 2, clsi_s, height=h, color="#B0B0B0", label="CLSI M100 2022")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Percent susceptible")
    ax.set_xlim(0, 100)
    ax.invert_yaxis()
    fig.legend(
        loc="upper center",
        ncol=2,
        bbox_to_anchor=(0.58, 1.02),
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return save_fig(fig, "Figure1_Susceptibility")


def fig_ptz_mic(analytic: pd.DataFrame):
    apply_journal_style()
    order = ["leq4", "le8", "exact8", "le16", "exact16", "exact32", "exact64", "hi"]
    tick = [r"$\leq$4", r"$\leq$8", "8", r"$\leq$16", "16", "32", "64", r"$>$64"]
    panels = [
        ("ESCHERICHIA COLI", r"A.  $\mathit{Escherichia\ coli}$"),
        ("KLEBSIELLA PNEUMONIAE", r"B.  $\mathit{Klebsiella\ pneumoniae}$"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.8), sharey=False)
    # Light gray: at or below the susceptible breakpoint, including the <=8 bound.
    # White with a hatch: <=16, a bound that can sit in susceptible or SDD.
    # Mid gray: exact 16, the SDD dilution. Dark: resistant.
    colors = [
        "#BDBDBD",
        "#BDBDBD",
        "#BDBDBD",
        "#FFFFFF",
        "#6B6B6B",
        "#2F2F2F",
        "#2F2F2F",
        "#2F2F2F",
    ]
    hatches = ["", "", "", "///", "", "", "", ""]
    for ax, (org, title) in zip(axes, panels):
        g = analytic.loc[
            (analytic["organism"] == org)
            & (analytic["antibiotic"] == "piperacillin_tazobactam")
        ].copy()
        g["bin"] = [
            mic_bin(v, i, p)
            for v, i, p in zip(g["AST_val1"], g["AST_inequality"], g["AST_panel"])
        ]
        g = g.loc[g["bin"].notna()]
        counts = g["bin"].value_counts()
        heights = [int(counts.get(b, 0)) for b in order]
        bars = ax.bar(
            np.arange(len(order)),
            heights,
            color=colors,
            width=0.82,
            edgecolor="#4D4D4D",
            linewidth=0.4,
        )
        for bar, hatch in zip(bars, hatches):
            bar.set_hatch(hatch)
        ax.set_xticks(np.arange(len(order)))
        ax.set_xticklabels(tick, fontsize=7)
        ax.set_title(title, loc="left", fontsize=9)
        ax.set_xlabel("Piperacillin-tazobactam MIC (µg/mL)")
        n = int(len(g))
        ax.text(0.98, 0.96, f"n = {n:,}", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    axes[0].set_ylabel("Number of results")
    from matplotlib.patches import Patch

    fig.legend(
        handles=[
            Patch(facecolor="#BDBDBD", edgecolor="#4D4D4D", label=r"Susceptible ($\leq$8/4)"),
            Patch(facecolor="#FFFFFF", edgecolor="#4D4D4D", hatch="///", label=r"Censored $\leq$16"),
            Patch(facecolor="#6B6B6B", edgecolor="#4D4D4D", label="SDD (exact 16/4)"),
            Patch(facecolor="#2F2F2F", edgecolor="#4D4D4D", label=r"Resistant ($\geq$32/4)"),
        ],
        loc="upper center",
        ncol=4,
        bbox_to_anchor=(0.5, 1.04),
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return save_fig(fig, "Figure2_PTZ_MIC")


def main():
    path = os.path.join(C.DERIVED_DIR, "ast_paired.parquet")
    cols = [
        "analytic",
        "organism",
        "antibiotic",
        "lab",
        "clsi",
        "site",
        "AST_val1",
        "AST_inequality",
        "AST_panel",
        "anon_id",
        "order_proc_id_coded",
        "order_time_jittered_utc_shifted",
    ]
    print("reading parquet", flush=True)
    ast = pd.read_parquet(path, columns=cols)
    analytic = ast.loc[ast["analytic"]].copy()
    del ast
    n_patients = int(analytic["anon_id"].nunique())
    n_cultures = int(analytic["order_proc_id_coded"].nunique())
    analytic["_t"] = pd.to_datetime(
        analytic["order_time_jittered_utc_shifted"], errors="coerce"
    )
    first_keys = (
        analytic.sort_values(["anon_id", "organism", "_t", "order_proc_id_coded"])
        .drop_duplicates(["anon_id", "organism"], keep="first")[
            ["anon_id", "organism", "order_proc_id_coded"]
        ]
    )
    first = analytic.merge(
        first_keys, on=["anon_id", "organism", "order_proc_id_coded"], how="inner"
    )
    print(
        f"first isolate {len(first):,} patients {first['anon_id'].nunique():,}",
        flush=True,
    )
    print(
        f"analytic {len(analytic):,} patients {n_patients:,} cultures {n_cultures:,}",
        flush=True,
    )

    summary = {}
    summary["all"] = pair_sir(analytic)
    for site in ["Urine", "Blood", "Respiratory"]:
        summary[f"site_{site.lower()}"] = pair_sir(analytic.loc[analytic["site"] == site])
    for org, lab in C.ORGANISM_LABEL.items():
        key = "org_" + lab.replace(". ", "_").replace(".", "")
        summary[key] = pair_sir(analytic.loc[analytic["organism"] == org])
        summary[key]["label"] = lab

    ptz_ent = analytic.loc[
        (analytic["antibiotic"] == "piperacillin_tazobactam")
        & (analytic["organism"].isin(ENTEROBACTERALES))
    ]
    summary["ptz_enterobacterales"] = pair_sir(ptz_ent)
    ptz_pa = analytic.loc[
        (analytic["antibiotic"] == "piperacillin_tazobactam")
        & (analytic["organism"] == "PSEUDOMONAS AERUGINOSA")
    ]
    summary["ptz_pa"] = pair_sir(ptz_pa)

    pair_rows = []
    for org, drug in PAIR_TABLE:
        g = analytic.loc[(analytic["organism"] == org) & (analytic["antibiotic"] == drug)]
        rec = pair_sir(g)
        rec["organism"] = org
        rec["organism_label"] = C.ORGANISM_LABEL[org]
        rec["antibiotic"] = drug
        rec["drug_label"] = C.DRUG_LABEL[drug]
        pair_rows.append(rec)
        summary[f"pair_{org}_{drug}"] = rec

    for org in ["ESCHERICHIA COLI", "KLEBSIELLA PNEUMONIAE"]:
        g = analytic.loc[
            (analytic["organism"] == org)
            & (analytic["antibiotic"] == "piperacillin_tazobactam")
        ]
        summary[f"ptz_{org}"] = pair_sir(g)

    other_ent = ENTEROBACTERALES - {"ESCHERICHIA COLI", "KLEBSIELLA PNEUMONIAE"}
    summary["ptz_other_enterobacterales"] = pair_sir(
        analytic.loc[
            (analytic["antibiotic"] == "piperacillin_tazobactam")
            & (analytic["organism"].isin(other_ent))
        ]
    )

    summary["first_all"] = pair_sir(first)
    ptz_ent_first = first.loc[
        (first["antibiotic"] == "piperacillin_tazobactam")
        & (first["organism"].isin(ENTEROBACTERALES))
    ]
    summary["first_ptz_enterobacterales"] = pair_sir(ptz_ent_first)
    ptz_bmd = ptz_ent.loc[
        ptz_ent["AST_panel"].astype(str).str.strip() == "broth_microdilution"
    ]
    summary["ptz_enterobacterales_bmd"] = pair_sir(ptz_bmd)
    summary["ptz_enterobacterales_bmd"]["n_parent"] = int(len(ptz_ent))
    for panel_name, key in [
        ("kirby_bauer", "ptz_enterobacterales_disk"),
        ("gradient_diffusion", "ptz_enterobacterales_gradient"),
    ]:
        summary[key] = pair_sir(
            ptz_ent.loc[ptz_ent["AST_panel"].astype(str).str.strip() == panel_name]
        )
    fox_first = first.loc[
        (first["organism"] == "STAPHYLOCOCCUS AUREUS")
        & (first["antibiotic"] == "cefoxitin")
    ]
    summary["first_sa_cefoxitin"] = pair_sir(fox_first)
    summary["ci_ptz"] = cluster_ci(ptz_ent)
    summary["ci_ptz_first"] = cluster_ci(ptz_ent_first)
    summary["ci_ptz_bmd"] = cluster_ci(ptz_bmd)
    summary["ptz_mic_notation"] = mic_notation(ptz_bmd)

    # Compact CSV for Table 1
    t1_rows = []
    order_keys = [
        ("All results", "all"),
        ("Urine", "site_urine"),
        ("Blood", "site_blood"),
        ("Respiratory", "site_respiratory"),
        ("E. coli", "org_E_coli"),
        ("K. pneumoniae", "org_K_pneumoniae"),
        ("P. aeruginosa", "org_P_aeruginosa"),
        ("S. aureus", "org_S_aureus"),
        ("PTZ, selected Enterobacterales", "ptz_enterobacterales"),
    ]
    # org keys as built
    summary["org_E_coli"] = pair_sir(analytic.loc[analytic["organism"] == "ESCHERICHIA COLI"])
    summary["org_K_pneumoniae"] = pair_sir(
        analytic.loc[analytic["organism"] == "KLEBSIELLA PNEUMONIAE"]
    )
    summary["org_P_aeruginosa"] = pair_sir(
        analytic.loc[analytic["organism"] == "PSEUDOMONAS AERUGINOSA"]
    )
    summary["org_S_aureus"] = pair_sir(
        analytic.loc[analytic["organism"] == "STAPHYLOCOCCUS AUREUS"]
    )

    for label, key in order_keys:
        rec = summary[key]
        t1_rows.append(
            {
                "stratum": label,
                "n": rec["n"],
                "lab_S": rec["lab"]["S_pct"],
                "lab_SDD": rec["lab"]["SDD_pct"],
                "lab_I": rec["lab"]["I_pct"],
                "lab_R": rec["lab"]["R_pct"],
                "clsi_S": rec["clsi"]["S_pct"],
                "clsi_SDD": rec["clsi"]["SDD_pct"],
                "clsi_I": rec["clsi"]["I_pct"],
                "clsi_R": rec["clsi"]["R_pct"],
            }
        )
    t1 = pd.DataFrame(t1_rows)
    t1_path = os.path.join(C.TABLES_DIR, "table1_categories.csv")
    t1.to_csv(t1_path, index=False)
    print("wrote", t1_path, flush=True)

    t2 = pd.DataFrame(
        [
            {
                "organism": r["organism_label"],
                "drug": r["drug_label"],
                "n": r["n"],
                "lab_S": r["lab"]["S_pct"],
                "clsi_S": r["clsi"]["S_pct"],
                "clsi_SDD": r["clsi"]["SDD_pct"],
                "lab_R": r["lab"]["R_pct"],
                "clsi_R": r["clsi"]["R_pct"],
                "s_to_sdd": r["s_to_sdd"],
                "s_to_r": r["s_to_r"],
                "r_to_s": r["r_to_s"],
            }
            for r in pair_rows
        ]
    )
    t2_path = os.path.join(C.TABLES_DIR, "table2_pairs.csv")
    t2.to_csv(t2_path, index=False)
    print("wrote", t2_path, flush=True)

    t3_spec = [
        ("E. coli", "ptz_ESCHERICHIA COLI"),
        ("K. pneumoniae", "ptz_KLEBSIELLA PNEUMONIAE"),
        ("Other selected Enterobacterales", "ptz_other_enterobacterales"),
        ("All selected Enterobacterales", "ptz_enterobacterales"),
        ("P. aeruginosa", "ptz_pa"),
    ]
    t3 = pd.DataFrame(
        [
            {
                "group": lab,
                "n": summary[k]["n"],
                "lab_S": summary[k]["lab"]["S_pct"],
                "clsi_S": summary[k]["clsi"]["S_pct"],
                "clsi_SDD": summary[k]["clsi"]["SDD_pct"],
                "lab_I": summary[k]["lab"]["I_pct"],
                "clsi_I": summary[k]["clsi"]["I_pct"],
                "lab_R": summary[k]["lab"]["R_pct"],
                "clsi_R": summary[k]["clsi"]["R_pct"],
                "s_to_sdd": summary[k]["s_to_sdd"],
                "s_to_sdd_pct": summary[k]["s_to_sdd_pct"],
            }
            for lab, k in t3_spec
        ]
    )
    t3_path = os.path.join(C.TABLES_DIR, "table3_ptz.csv")
    t3.to_csv(t3_path, index=False)
    print("wrote", t3_path, flush=True)

    slim = {
        "_cohort": {
            "n": int(len(analytic)),
            "patients": n_patients,
            "cultures": n_cultures,
        }
    }
    for k, v in summary.items():
        if isinstance(v, dict) and "lab" in v and "S_pct" in v["lab"]:
            slim[k] = {
                "n": v["n"],
                "lab_S": v["lab"]["S_pct"],
                "lab_SDD": v["lab"]["SDD_pct"],
                "lab_I": v["lab"]["I_pct"],
                "lab_R": v["lab"]["R_pct"],
                "lab_NS": v["lab"]["NS_pct"],
                "clsi_S": v["clsi"]["S_pct"],
                "clsi_SDD": v["clsi"]["SDD_pct"],
                "clsi_I": v["clsi"]["I_pct"],
                "clsi_R": v["clsi"]["R_pct"],
                "clsi_NS": v["clsi"]["NS_pct"],
                "s_to_sdd": v.get("s_to_sdd"),
                "s_to_sdd_pct": v.get("s_to_sdd_pct"),
                "s_to_r": v.get("s_to_r"),
                "s_to_r_pct": v.get("s_to_r_pct"),
                "r_to_s": v.get("r_to_s"),
                "r_to_s_pct": v.get("r_to_s_pct"),
                "delta_s": v.get("delta_s"),
            }
    slim["_checks"] = {
        "first_n": summary["first_all"]["n"],
        "first_lab_S": summary["first_all"]["lab"]["S_pct"],
        "first_clsi_S": summary["first_all"]["clsi"]["S_pct"],
        "first_clsi_SDD": summary["first_all"]["clsi"]["SDD_pct"],
        "first_ptz_n": summary["first_ptz_enterobacterales"]["n"],
        "first_ptz_lab_S": summary["first_ptz_enterobacterales"]["lab"]["S_pct"],
        "first_ptz_clsi_S": summary["first_ptz_enterobacterales"]["clsi"]["S_pct"],
        "first_ptz_clsi_SDD": summary["first_ptz_enterobacterales"]["clsi"]["SDD_pct"],
        "ptz_bmd_n": summary["ptz_enterobacterales_bmd"]["n"],
        "ptz_bmd_parent_n": summary["ptz_enterobacterales_bmd"].get("n_parent"),
        "ptz_bmd_clsi_SDD": summary["ptz_enterobacterales_bmd"]["clsi"]["SDD_pct"],
        "ptz_bmd_clsi_S": summary["ptz_enterobacterales_bmd"]["clsi"]["S_pct"],
        "ptz_bmd_lab_S": summary["ptz_enterobacterales_bmd"]["lab"]["S_pct"],
        "first_sa_fox_n": summary["first_sa_cefoxitin"]["n"],
        "first_sa_fox_lab_S": summary["first_sa_cefoxitin"]["lab"]["S_pct"],
        "first_sa_fox_clsi_S": summary["first_sa_cefoxitin"]["clsi"]["S_pct"],
        "ptz_disk_n": summary["ptz_enterobacterales_disk"]["n"],
        "ptz_disk_sdd": summary["ptz_enterobacterales_disk"]["clsi"]["SDD_pct"],
        "ptz_grad_n": summary["ptz_enterobacterales_gradient"]["n"],
        "ptz_grad_sdd": summary["ptz_enterobacterales_gradient"]["clsi"]["SDD_pct"],
        "ci_ptz": summary["ci_ptz"],
        "ci_ptz_first": summary["ci_ptz_first"],
        "ci_ptz_bmd": summary["ci_ptz_bmd"],
        "ptz_mic_notation": summary["ptz_mic_notation"],
    }
    num_path = os.path.join(C.TABLES_DIR, "numbers.json")
    with open(num_path, "w", encoding="utf-8") as f:
        json.dump(slim, f, indent=2)
    print("wrote", num_path, flush=True)
    print("CHECKS", json.dumps(slim["_checks"], indent=2), flush=True)

    print("\nTABLE 1", flush=True)
    print(t1.to_string(index=False), flush=True)
    print("\nTABLE 2", flush=True)
    print(t2.to_string(index=False), flush=True)
    print("\nTABLE 3", flush=True)
    print(t3.to_string(index=False), flush=True)

    pdf1 = fig_susceptibility(analytic)
    pdf2 = fig_ptz_mic(analytic)
    ms_dir = os.path.join(C.PROJECT_DIR, "manuscript", "dmid")
    if os.path.isdir(ms_dir):
        shutil.copy2(pdf1, os.path.join(ms_dir, "Figure_1.pdf"))
        shutil.copy2(pdf2, os.path.join(ms_dir, "Figure_2.pdf"))
        for stale in ("Figure_3.pdf", "Figure_4.pdf"):
            p = os.path.join(ms_dir, stale)
            if os.path.exists(p):
                os.remove(p)
                print("removed", stale, flush=True)
        print("copied figures into manuscript folder", flush=True)


if __name__ == "__main__":
    main()
