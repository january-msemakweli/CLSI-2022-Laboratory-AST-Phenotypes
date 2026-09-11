# Laboratory-reported versus CLSI M100 2022 AST phenotypes (ARMD-MGB)

Analysis code for a result-level comparison of laboratory-reported antimicrobial
susceptibility categories with CLSI M100 2022 interpretations of the same
recorded MIC or zone in the Antibiotic Resistance Microbiology Dataset from
Mass General Brigham (ARMD-MGB) v1.0.0.

This archive accompanies the manuscript *Impact of CLSI M100 2022 Breakpoint
Changes on the Susceptibility Reporting for Enterobacterales and Staphylococcus
aureus*.

ARMD-MGB source files are not included.

## Data availability

This repository contains **no patient-level data**. ARMD-MGB v1.0.0 is available
to credentialed researchers from PhysioNet
(https://physionet.org/content/armd-mgb/1.0.0/) and cannot be redistributed
under the data use agreement. The `derived/` directory (paired AST parquet) is
patient-level and is excluded by `.gitignore`. Only code and non-identifiable
aggregate outputs (tables and figures) are included.

## Repository layout

```
scripts/            Python analysis (`01_analysis.py`, `02_figures.py`, `config.py`)
tables/             Non-identifiable aggregate tables
figures/            Publication figures (PDF and PNG)
requirements.txt    Python dependencies
```

## Reproducing the analysis

1. Obtain ARMD-MGB v1.0.0 from PhysioNet and place the CSV files at
   `../physionet.org/files/armd-mgb/1.0.0/` relative to this project (see
   `scripts/config.py`, `ARMD_DIR`). Adjust that path if your layout differs.
2. Create the environment:

```bash
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

3. Run the pipeline from the project root:

```bash
python scripts/01_analysis.py
python scripts/02_figures.py
```

`01_analysis.py` writes `derived/ast_paired.parquet` (git-ignored).
`02_figures.py` writes the three manuscript tables to `tables/` and the two
figures to `figures/`.

## Data use

ARMD-MGB is de-identified and available after PhysioNet credentialing and
required training. A PhysioNet data-use agreement applies. This code does not
redistribute those files.

## Citation

Zenodo: [https://doi.org/10.5281/zenodo.22699153](https://doi.org/10.5281/zenodo.22699153)

GitHub: [https://github.com/january-msemakweli/CLSI-2022-Laboratory-AST-Phenotypes](https://github.com/january-msemakweli/CLSI-2022-Laboratory-AST-Phenotypes)

See `CITATION.cff`.

## License

MIT. See `LICENSE`. ARMD-MGB data are governed by the PhysioNet Credentialed
Health Data License.
