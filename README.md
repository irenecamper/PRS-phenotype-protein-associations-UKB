# 🌿 carpe-diem-ukb

> Polygenic risk score (PRS) analysis pipeline for UK Biobank data — association with phenotypes and plasma proteomics.

---

## Overview

`carpe-diem-ukb` is a reproducible R pipeline for:

1. Loading and harmonising UK Biobank (UKB) phenotypic, genetic, and proteomic data
2. Running PRS × phenotype association analyses (linear & logistic regression)
3. Running PRS × plasma proteome analyses (Olink NPX data)
4. Producing forest plots, volcano plots, and summary Excel workbooks

The analysis entry-point is `notebook.Rmd`. All helper functions live in `R/utils.R`. Package management is handled by [`renv`](https://rstudio.github.io/renv/).

---
## Set up

```
git clone https://github.com/irenecamper/PRS-phenotype-protein-associations-UKB
cd PRS-phenotype-protein-associations-UKB
```

## Repository structure

```
carpe-diem-ukb/
├── notebook.Rmd
│   └── Main analysis notebook (knit to HTML)
├── R/
│   ├── utils.R
│   │   └── Helper functions for models, plots, and IO
│   └── packages.R
│       └── renv bootstrap and package loading
├── inputs/
│   ├── ukb_metadata.csv
│   │   └── Trait definitions, ICD-10 codes, and transformations
│   ├── dummy_data_participant.csv
│   │   └── Example/schema participant input
│   ├── dummy_PRS_pred_allchr.sscore
│   │   └── Example/schema PRS input
│   ├── dummy_data_diag_hesin.csv
│   │   └── Example/schema diagnosis input
│   ├── generate_example_inputs.R
│   │   └── Script to generate dummy/example .rds files
│   ├── dummy_protein_panel_map.rds
│   │   └── Example/schema protein panel map
│   └── dummy_pro_batch_fastingtime.rds
│       └── Example/schema protein input
├── figures/
│   └── Output figures [git-ignored]
├── results/
│   └── Output Excel tables [git-ignored]
├── renv.lock
│   └── Exact R package versions
└── .gitignore
```

> **Private inputs** (raw UKB data) go in `inputs/` and are listed in `.gitignore`.

---

## Quickstart

### 1. Clone and set up the R environment

```r
# install renv if needed
install.packages("renv")

# restore exact package versions from renv.lock
renv::restore()
```

### 2. Provide your input data

Place the following files in `inputs/` (see `inputs/examples/` for expected column layouts):

| File | Description |
|------|-------------|
| `data_participant.csv` | UKB participant phenotype table (see column list below) |
| `PRS_pred_allchr.sscore` | PLINK2 `.sscore` output with `IID` and `SCORE1_SUM` |
| `data_diag_hesin.csv` | Long-format ICD-10 diagnoses (`eid`, `diag_icd10`) |
| `pro_batch_fastingtime.rds` | Olink NPX proteomic data with batch and fasting time |
| `protein_panel_map.rds` | Olink assay → panel mapping |

To generate example `.rds` schema files locally:

```r
source("inputs/generate_dummy_inputs.R")
```

### 3. Run the analysis

Open `notebook.Rmd` in RStudio and knit, or run from the command line:

```bash
Rscript -e "rmarkdown::render('notebook.Rmd')"
```

### Adding a new trait

To add a new phenotype, append a row to `inputs/ukb_metadata.csv` following the schema above — no R code changes required. For a new continuous UKB field:

```
direct;p12345_i0;my_trait;My trait label;Trait group;qn;;;continuous
```

For a new ICD-10 binary flag:

```
icd10_flag;diag_icd10;has_X;Disease X;Disease group;binary;X01,X02;0.0;binary
```

---

## Data access

Raw data files are **not** included in this repository and must never be committed (enforced via `.gitignore`).

---

## Reproducibility

- R package versions are locked in `renv.lock`; restore with `renv::restore()`
- All outputs (figures, tables) are stamped with the run date in their filenames

## R version

This project was developed with R >= 4.5.0.
Minor version differences should not affect results.