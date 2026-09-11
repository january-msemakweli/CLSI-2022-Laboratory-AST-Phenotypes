"""Build the paired laboratory vs CLSI 2022 AST parquet from ARMD-MGB."""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C


def map_pheno(val):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    t = str(val).strip().lower().replace("_", " ").replace("-", " ")
    t = " ".join(t.split())
    if t == "susceptible":
        return "S"
    if t in ("susceptible dose dependent", "sdd"):
        return "SDD"
    if t == "intermediate":
        return "I"
    if t == "resistant":
        return "R"
    if t in ("non susceptible", "nonsusceptible"):
        return "NS"
    return None


def load_micro():
    path = C.armd(C.MICRO_FILE)
    print(f"reading {path}", flush=True)
    df = pd.read_csv(
        path,
        usecols=C.MICRO_COLS,
        dtype=str,
        keep_default_na=False,
        na_values=[],
    )
    print(f"  micro rows={len(df):,}", flush=True)
    return df


def build_cohort(micro: pd.DataFrame) -> pd.DataFrame:
    micro = micro.copy()
    micro["lab"] = micro["AST_pheno"].map(map_pheno)
    micro["clsi"] = micro["CLSI_2022_pheno"].map(map_pheno)
    ast = micro.loc[micro["lab"].notna() & micro["clsi"].notna()].copy()
    ast["prelim"] = ast["prelim_AST"].str.strip() == "X"
    ast["analytic"] = ~ast["prelim"]
    ast["site"] = ast["culture_description"].map(lambda x: C.SITE_LABEL.get(x, x))
    ast["org_label"] = ast["organism"].map(
        lambda x: C.ORGANISM_LABEL.get(x, x.title())
    )
    ast["drug_label"] = ast["antibiotic"].map(
        lambda x: C.DRUG_LABEL.get(x, x.replace("_", "-"))
    )
    n = int(ast["analytic"].sum())
    print(f"  paired={len(ast):,} analytic={n:,}", flush=True)
    return ast


def main():
    micro = load_micro()
    ast = build_cohort(micro)
    del micro
    ast_path = os.path.join(C.DERIVED_DIR, "ast_paired.parquet")
    ast.to_parquet(ast_path, index=False)
    print(f"  wrote {ast_path}", flush=True)


if __name__ == "__main__":
    main()
