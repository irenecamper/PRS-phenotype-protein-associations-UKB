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

---

## Input file schemas

### `data_participant.csv`

One row per participant. Key columns (full list in `inputs/examples/data_participant.csv`):

| Column | Description |
|--------|-------------|
| `eid` | UK Biobank participant identifier |
| `p21001_i0` | BMI |
| `p31` | Sex (0 = female, 1 = male) |
| `p21003_i0` | Age at assessment |
| `p22009_a1`–`p22009_a20` | Genetic principal components 1–20 |
| `p4080_i0_a0` / `p4079_i0_a0` | Systolic / diastolic blood pressure |
| `p30690_i0`, `p30760_i0`, etc. | Blood biochemistry (cholesterol, HDL, etc.) |
| `p6149_i0` | Dental prosthetics (includes dentures) |
| `p26414` / `p26411` | Education / income Townsend scores |

### `PRS_pred_allchr.sscore`

Tab-separated PLINK2 score file:

```
#FID    IID    SCORE1_SUM
0       1000001    0.00412
```

### `data_diag_hesin.csv`

Long-format hospital episode statistics (one row per diagnosis):

```
eid,diag_icd10
1000001,E11
1000001,I25
```

### `protein_panel_map.rds`

R data frame with columns: `coding`, `Assay`, `Panel`, `Panel_clean`.

### `pro_batch_fastingtime.rds`

R data frame with columns: `eid`, `ins_index`, `PlateID`, `Batch`, `Fasting_time`, then one column per protein (lower-case Olink assay name, e.g. `il6`, `gdf15`, …).

---

## Metadata file (`inputs/ukb_metadata.csv`)

**This file is the single source of truth for all traits in the pipeline.** It is semicolon-delimited and committed to the repo — editing it is the only thing you need to do to add, rename, or remove a phenotype.

### Column definitions

| Column | Values | Description |
|--------|--------|-------------|
| `source_type` | `direct` \| `icd10_flag` | How the variable is derived. `direct` = taken straight from a UKB field column; `icd10_flag` = binary flag built by searching `data_diag_hesin.csv` for ICD-10 codes |
| `ukb_field` | UKB field ID or `diag_icd10` | For `direct` rows: the exact column name in `data_participant.csv` that gets renamed to `variable`. For `icd10_flag` rows: the column in the diagnosis file to search (always `diag_icd10`) |
| `variable` | string | Internal R variable name used throughout the pipeline (after renaming) |
| `label` | string | Human-readable label shown on plot axes and table headers |
| `trait_group` | string | Phenotype category used for colour-coding and faceting in forest / volcano plots |
| `transformation` | `none` \| `qn` \| `log` \| `boxcox` \| `binary` | Transformation applied to the outcome before regression. `qn` = quantile-normal; `binary` = coerce to 0/1 numeric |
| `icd10_values` | comma-separated ICD-10 codes | Only used for `icd10_flag` rows. Codes are matched at 3-character level (e.g. `E11` matches `E110`, `E119`, …) |
| `fill_missing` | `0` \| _(empty)_ | For `icd10_flag` rows: participants not found in the diagnosis file are assumed disease-free (`0`). Leave empty if missingness should be preserved as `NA` |
| `trait_type` | `continuous` \| `binary` \| `none` | Determines which regression model is used (`lm` vs `glm(family=binomial)`). Rows with `none` are covariates only (e.g. genetic PCs) |

### How the pipeline uses this file

```
load_maps(METADATA_CSV)          # builds all lookup maps
  ├── rename_map       # UKB field ID  →  internal variable name
  ├── var_name_map     # variable      →  plot label
  ├── var_trans_map    # variable      →  transformation method
  ├── var_category_map # variable      →  trait_group
  └── var_code_map     # variable      →  UKB field or ICD-10 codes

safe_icd_join(...)               # creates binary ICD-10 flag columns
  └── uses source_type == "icd10_flag" rows + icd10_values + fill_missing
```

### Full trait inventory

#### Continuous outcomes (`trait_type = continuous`)

| `variable` | `label` | `trait_group` | UKB field | Transformation |
|------------|---------|---------------|-----------|----------------|
| `bmi` | BMI | Adiposity | p21001_i0 | qn |
| `waist` | Waist | Adiposity | p48_i0 | qn |
| `hip` | Hip | Adiposity | p49_i0 | qn |
| `body_fat_perc` | Body fat % | Adiposity | p23099_i0 | qn |
| `glucose` | Glucose | Glycemic | p30740_i0 | qn |
| `hba1c` | HbA1c | Glycemic | p30750_i0 | qn |
| `chol` | Total cholesterol | Lipid | p30690_i0 | qn |
| `ldl` | LDL cholesterol | Lipid | p23405_i0 | qn |
| `hdl` | HDL cholesterol | Lipid | p30760_i0 | qn |
| `tgs` | Triglycerides | Lipid | p30870_i0 | qn |
| `vldl` | VLDL cholesterol | Lipid | p23403_i0 | qn |
| `alat` | ALAT | Liver | p30620_i0 | qn |
| `ggt` | GGT | Liver | p30730_i0 | qn |
| `alp` | ALP | Liver | p30610_i0 | qn |
| `pdff` | PDFF | Liver | p40061_i2 | qn |
| `ct1` | Liver iron corrected T1 | Liver | p40062_i2 | qn |
| `sbp` | SBP | Cardiovascular | p4080_i0_a0 | qn |
| `dbp` | DBP | Cardiovascular | p4079_i0_a0 | qn |
| `crp` | C-reactive protein | Cardiovascular | p30710_i0 | qn |
| `education` | Education score | Sociodemographic factor | p26414 | qn |
| `income` | Income score | Sociodemographic factor | p26411 | qn |

#### Binary outcomes (`trait_type = binary`)

| `variable` | `label` | `trait_group` | Source | ICD-10 codes |
|------------|---------|---------------|--------|--------------|
| `has_cirrhosis` | Cirrhosis | Metabolic disease | ICD-10 flag | K74 |
| `has_ald` | Alcohol related Liver disease | Metabolic disease | ICD-10 flag | K70 |
| `has_fatty_Liver` | Fatty (change of) Liver | Metabolic disease | ICD-10 flag | K76.0 |
| `has_t2d` | Type 2 diabetes | Metabolic disease | ICD-10 flag | E11 |
| `has_mi` | Myocardial infarction | Cardiovascular disease | ICD-10 flag | I21, I22 |
| `has_stroke` | Stroke | Cardiovascular disease | ICD-10 flag | I63, I64 |
| `has_cad` | Coronary artery disease | Cardiovascular disease | ICD-10 flag | I25 |
| `has_dentures` | Dentures | Dental caries | Direct (p6149_i0) | — |

#### Covariates / identifiers (`trait_type = none`)

| `variable` | `label` | UKB field |
|------------|---------|-----------|
| `eid` | PATID | eid |
| `sex` | Sex | p31 |
| `age` | Age | p21003_i0 |
| `pc1`–`pc20` | Genetic PC 1–20 | p22009_a1–p22009_a20 |

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