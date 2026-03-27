############################################################
# packages.R
# purpose: restore the renv environment and load all packages
############################################################


# ── 0. Ensure renv is available ───────────────────────────────────────────────

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}


# ── 1. Activate & restore renv ───────────────────────────────────────────────

if (file.exists("renv/activate.R")) {
  source("renv/activate.R")   # activate project library
  
  if (requireNamespace("renv", quietly = TRUE)) {
    renv::restore(prompt = FALSE)
  } else {
    warning("[packages.R] renv not available after activation.")
  }
  
} else {
  message(
    "[packages.R] renv not initialised. ",
    "Run `renv::init()` once to set up reproducible package management."
  )
}


# ── 2. Required packages ─────────────────────────────────────────────────────

required_packages <- c(
  # data wrangling
  "dplyr",
  "tidyr",
  "purrr",
  # file I/O
  "readxl",
  "writexl",
  "openxlsx",
  # visualisation
  "ggplot2",
  "ggrepel",
  "patchwork",
  # statistics
  "MASS",
  "VGAM",
  # utilities
  "lubridate",
  # "fuzzyjoin",
  "styler"
)


# ── 3. Install missing packages (fallback) ───────────────────────────────────

.install_if_missing <- function(pkgs) {
  installed <- rownames(installed.packages())
  missing   <- setdiff(pkgs, installed)
  
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  }
}

.install_if_missing(required_packages)


# ── 4. Load packages ─────────────────────────────────────────────────────────

invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

message("All packages loaded.")