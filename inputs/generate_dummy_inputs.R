#!/usr/bin/env Rscript
# generate_example_inputs.R
#
# Run this script ONCE to create the two example .rds files that live in
# inputs.  The real files (inputs/) are private UKB data and are
# listed in .gitignore.
#
# Usage:
#   Rscript inputs/generate_dummy_inputs.R
# or from inside R:
#   source("inputs/generate_dummy_inputs.R")

set.seed(42)

out_dir <- file.path("inputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── 1. protein_panel_map.rds ──────────────────────────────────────────────────
# Columns: coding | Assay | Panel | Panel_clean
# 'Assay' is the lower-case protein name used as column header in the batch
# data. 'coding' is the gene-level identifier (same or HGNC alias).

protein_panel_map <- data.frame(
  coding     = c("il6", "tnf", "crp", "igf1", "gdf15",
                 "nt_probnp", "troponin_i", "vegfa",
                 "nfl", "gfap", "snap25", "uchl1",
                 "psa", "ca125", "cea", "afp"),
  Assay      = c("il6", "tnf", "crp", "igf1", "gdf15",
                 "nt_probnp", "troponin_i", "vegfa",
                 "nfl", "gfap", "snap25", "uchl1",
                 "psa", "ca125", "cea", "afp"),
  Panel      = c(rep("Inflammation", 4),
                 rep("Cardiometabolic", 4),
                 rep("Neurology", 4),
                 rep("Oncology", 4)),
  Panel_clean = c(rep("Inflammation", 4),
                  rep("Cardiometabolic", 4),
                  rep("Neurology", 4),
                  rep("Oncology", 4)),
  stringsAsFactors = FALSE
)

saveRDS(protein_panel_map,
        file.path(out_dir, "dummy_protein_panel_map.rds"))
message("Saved: ", file.path(out_dir, "dummy_protein_panel_map.rds"))

# ── 2. pro_batch_fastingtime.rds ──────────────────────────────────────────────
# Columns: eid | ins_index | PlateID | Batch | Fasting_time | <protein cols>
# One row per participant; protein values are approximately normal on the
# NPX (log2) scale used by Olink.

n <- 3  # keep tiny — this is just a schema example

proteins <- protein_panel_map$Assay

pro_batch <- data.frame(
  eid          = c(1000001L, 1000002L, 1000003L),
  ins_index    = c(0L, 0L, 0L),
  PlateID      = c("P001", "P001", "P002"),
  Batch        = c(1L, 1L, 2L),
  Fasting_time = c(8.5, 12.0, 6.0)
)

# Add fake NPX protein columns (Gaussian noise around 5 NPX units)
protein_mat <- matrix(
  rnorm(n * length(proteins), mean = 5, sd = 1.2),
  nrow = n,
  dimnames = list(NULL, proteins)
)

pro_batch <- cbind(pro_batch, as.data.frame(protein_mat))

saveRDS(pro_batch,
        file.path(out_dir, "dummy_pro_batch_fastingtime.rds"))
message("Saved: ", file.path(out_dir, "dummy_pro_batch_fastingtime.rds"))

message("\nDone. Example .rds files written to '", out_dir, "/'.")
