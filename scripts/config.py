"""Shared paths and labels for the paired AST analysis."""

from __future__ import annotations

import os

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE_DIR = os.path.dirname(PROJECT_DIR)
ARMD_DIR = os.path.join(
    WORKSPACE_DIR, "physionet.org", "files", "armd-mgb", "1.0.0"
)

DERIVED_DIR = os.path.join(PROJECT_DIR, "derived")
TABLES_DIR = os.path.join(PROJECT_DIR, "tables")
FIGURES_DIR = os.path.join(PROJECT_DIR, "figures")

for _d in (DERIVED_DIR, TABLES_DIR, FIGURES_DIR):
    os.makedirs(_d, exist_ok=True)


def armd(name: str) -> str:
    return os.path.join(ARMD_DIR, name)


MICRO_FILE = "microbiology_cohort_deid_tj_updated.csv"

MICRO_COLS = [
    "anon_id",
    "order_proc_id_coded",
    "order_time_jittered_utc_shifted",
    "culture_description",
    "organism",
    "prelim_AST",
    "AST_panel",
    "AST_inequality",
    "AST_val1",
    "antibiotic",
    "AST_pheno",
    "CLSI_2022_pheno",
]

ORGANISM_LABEL = {
    "ESCHERICHIA COLI": "E. coli",
    "KLEBSIELLA PNEUMONIAE": "K. pneumoniae",
    "PSEUDOMONAS AERUGINOSA": "P. aeruginosa",
    "STAPHYLOCOCCUS AUREUS": "S. aureus",
}

DRUG_LABEL = {
    "piperacillin_tazobactam": "Piperacillin-tazobactam",
    "ceftriaxone": "Ceftriaxone",
    "cefazolin": "Cefazolin",
    "cefepime": "Cefepime",
    "meropenem": "Meropenem",
    "ciprofloxacin": "Ciprofloxacin",
    "levofloxacin": "Levofloxacin",
    "trimethoprim_sulfamethoxazole": "Trimethoprim-sulfamethoxazole",
    "oxacillin": "Oxacillin",
    "cefoxitin": "Cefoxitin",
    "vancomycin": "Vancomycin",
    "clindamycin": "Clindamycin",
}

SITE_LABEL = {
    "URINE": "Urine",
    "BLOOD": "Blood",
    "RESPIRATORY_TRACT": "Respiratory",
}
