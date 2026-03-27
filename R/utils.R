make_theme <- function(base_size = 16) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(fill = NA, colour = "black"),
      plot.title         = ggplot2::element_text(size = base_size + 2, hjust = 0.5),
      legend.position    = "right",
      legend.title       = ggplot2::element_text(size = base_size),
      legend.text        = ggplot2::element_text(size = base_size - 2),
      axis.title         = ggplot2::element_text(size = base_size),
      axis.text          = ggplot2::element_text(size = base_size - 2),
      strip.text         = ggplot2::element_text(size = base_size, hjust = 0, vjust = 0.5, angle = 0),
      strip.background   = ggplot2::element_blank()
    )
}

load_maps <- function(metadata_csv) {
  meta <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE, sep = ";")
  list(
    metadata         = meta,
    rename_map       = setNames(
      meta$variable[meta$source_type == "direct"],
      meta$ukb_field[meta$source_type == "direct"]
    ),
    var_name_map     = setNames(meta$label,            meta$variable),
    var_trans_map    = setNames(meta$transformation,   meta$variable),
    var_category_map = setNames(meta$trait_group,      meta$variable),
    var_code_map     = setNames(
      ifelse(meta$source_type == "direct", meta$ukb_field, meta$icd10_values),
      meta$variable
    ),
    var_type_map     = setNames(meta$label, meta$trait_type)
  )
}

normalize_icd10_3 <- function(x) {
  x <- toupper(trimws(x))
  x <- sub(" .*", "", x)
  sub("^([A-Z][0-9][0-9]).*$", "\\1", x)
}

split_codes <- function(values) {
  vals <- trimws(unlist(strsplit(values, ",", fixed = TRUE)))
  vals <- toupper(vals)
  vals[vals != ""]
}

make_icd_flags_from_metadata <- function(data_icd10, metadata_csv, id_col = "eid") {
  meta <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE, sep = ";")
  meta <- meta[meta$source_type == "icd10_flag", , drop = FALSE]

  if (is.null(data_icd10) || nrow(data_icd10) == 0 || nrow(meta) == 0)
    return(data.frame(eid = numeric(0)))

  if (!id_col %in% names(data_icd10))
    stop("missing id column in data_icd10: ", id_col)

  source_cols <- unique(meta$ukb_field)
  if (length(source_cols) != 1)
    stop("all icd10_flag rows must use the same ukb_field source column.")

  code_col <- source_cols[1]
  if (!code_col %in% names(data_icd10)) {
    if (code_col == "diag_icd10" && "icd10_code" %in% names(data_icd10)) {
      code_col <- "icd10_code"
    } else {
      stop("ICD10 source column not found in data_icd10: ", source_cols[1])
    }
  }

  data_icd10[[code_col]] <- normalize_icd10_3(data_icd10[[code_col]])
  out <- data.frame(eid = unique(as.numeric(data_icd10[[id_col]])))

  for (i in seq_len(nrow(meta))) {
    vals <- split_codes(meta$icd10_values[i])
    ids  <- unique(as.numeric(data_icd10[[id_col]][data_icd10[[code_col]] %in% vals]))
    out[[meta$variable[i]]] <- as.integer(out$eid %in% ids)
  }
  out
}

fill_missing_flags_from_metadata <- function(data, metadata_csv) {
  meta <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE, sep = ";")
  meta <- meta[
    meta$source_type == "icd10_flag" & !is.na(meta$fill_missing) & meta$fill_missing == 0,
    , drop = FALSE
  ]
  for (v in meta$variable) {
    if (!v %in% names(data)) {
      data[[v]] <- 0L
    } else {
      data[[v]][is.na(data[[v]])] <- 0L
    }
  }
  data
}

safe_icd_join <- function(data, data_icd10, metadata_csv, id_col = "eid") {
  meta     <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE, sep = ";")
  icd_vars <- meta$variable[meta$source_type == "icd10_flag"]

  if (length(icd_vars) == 0) return(data)

  if (!is.null(data_icd10) && nrow(data_icd10) > 0) {
    icd_flags <- make_icd_flags_from_metadata(data_icd10, metadata_csv, id_col)
    if (nrow(icd_flags) > 0) {
      data[[id_col]] <- as.numeric(data[[id_col]])
      data <- dplyr::left_join(data, icd_flags, by = id_col)
    } else {
      for (v in icd_vars) if (!v %in% names(data)) data[[v]] <- 0L
    }
  } else {
    for (v in icd_vars) if (!v %in% names(data)) data[[v]] <- 0L
  }

  fill_missing_flags_from_metadata(data, metadata_csv)
}

make_fig_filename <- function(datatype = NULL, plot_type, ext = "png") {
  stamp <- format(Sys.time(), "%Y_%m_%d")
  parts <- Filter(function(x) !is.null(x) && nzchar(x), list(datatype, plot_type, stamp))
  paste0(paste(parts, collapse = "_"), ".", ext)
}

make_table_filename <- function(datatype, ext = "xlsx") {
  stamp <- format(Sys.time(), "%Y_%m_%d")
  paste0(datatype, "_", stamp, ".", ext)
}

transform_variable <- function(x, method = "log") {
  switch(method,
    none   = x,
    binary = as.numeric(x),
    qn     = stats::qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))),
    log = {
      if (any(x <= 0, na.rm = TRUE)) x <- x + abs(min(x, na.rm = TRUE)) + 1
      log(x)
    },
    boxcox = {
      if (any(x <= 0, na.rm = TRUE)) x <- x + abs(min(x, na.rm = TRUE)) + 1
      bc     <- MASS::boxcox(x ~ 1, plotit = FALSE)
      lambda <- bc$x[which.max(bc$y)]
      ifelse(lambda == 0, log(x), (x^lambda - 1) / lambda)
    },
    stop("Invalid transformation method: '", method, "'")
  )
}

scale_variable <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x)
  as.numeric(scale(x))
}

scale_from_map <- function(x, varname, trans_map) {
  method <- unname(trans_map[varname])
  if (!is.na(method) && length(method) > 0 && method == "binary") return(as.numeric(x))
  if (!is.numeric(x)) return(x)
  scale_variable(x)
}

.apply_transform <- function(x, varname, trans_map) {
  method <- unname(trans_map[varname])
  if (length(method) == 0 || is.na(method) || method == "none") return(x)
  if (method == "binary") {
    if (is.factor(x)) x <- as.character(x)
    return(as.numeric(x))
  }
  if (!is.numeric(x)) {
    warning("Variable '", varname, "' is not numeric but transformation '", method, "' was requested. Returning unchanged.")
    return(x)
  }
  transform_variable(x, method = method)
}

run_lm_grid <- function(df,
                        outcomes,
                        predictors,
                        covariates   = c("age", "sex"),
                        idcol,
                        conf_level   = 0.95,
                        trait_map    = setNames(character(0), character(0)),
                        trans_map    = setNames(character(0), character(0)),
                        panel_lookup = NULL,
                        model_filter = c("both", "linear", "logistic")) {

  model_filter <- match.arg(model_filter)

  empty_out <- function() {
    data.frame(
      outcome = character(0), predictor = character(0),
      trait_category = character(0),
      model = character(0), N = integer(0),
      beta = numeric(0), R2 = numeric(0),
      pvalue = numeric(0), pvalue_fdr = numeric(0),
      ci_upper = numeric(0), ci_lower = numeric(0),
      se = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  if (length(outcomes) == 0 || length(predictors) == 0) return(empty_out())

  req <- unique(c(outcomes, predictors, covariates, idcol))
  if (!all(req %in% names(df))) {
    warning("missing columns: ", paste(setdiff(req, names(df)), collapse = ", "))
    return(empty_out())
  }

  fit_one <- function(dat, out, pred) {
    rhs_terms      <- unique(c(pred, covariates))
    panel_i        <- NULL
    sample_age_var <- NULL

    if (!is.null(panel_lookup)) {
      panel_i <- panel_lookup[[out]]
      if (!is.null(panel_i)) {
        sample_age_var <- paste0("sample_age_", panel_i)
        if (sample_age_var %in% names(dat)) rhs_terms <- unique(c(rhs_terms, sample_age_var))
      }
    }

    needed <- unique(c(out, rhs_terms, idcol))
    dat2   <- dat[, needed, drop = FALSE]
    dat2   <- dat2[!duplicated(dat2[[idcol]]), , drop = FALSE]

    dat2[[out]]  <- .apply_transform(dat2[[out]],  out,  trans_map)
    dat2[[pred]] <- .apply_transform(dat2[[pred]], pred, trans_map)
    for (cv in rhs_terms[rhs_terms != pred])
      dat2[[cv]] <- .apply_transform(dat2[[cv]], cv, trans_map)

    dat2 <- dat2[stats::complete.cases(dat2[, c(out, rhs_terms), drop = FALSE]), , drop = FALSE]
    if (nrow(dat2) < 3) return(NULL)

    is_binary <- isTRUE(unname(trans_map[out]) == "binary")

    if (model_filter == "linear"   &&  is_binary) return(NULL)
    if (model_filter == "logistic" && !is_binary) return(NULL)

    if (is_binary) {
      vals <- sort(unique(dat2[[out]]))
      if (!all(vals %in% c(0, 1))) { warning("binary outcome '", out, "' contains non-0/1 values; skipping."); return(NULL) }
      if (length(vals) < 2)        { warning("binary outcome '", out, "' has only one class; skipping.");    return(NULL) }
    }

    fml <- stats::as.formula(paste(out, "~", paste(rhs_terms, collapse = " + ")))
    fit <- try(
      if (is_binary) stats::glm(fml, data = dat2, family = stats::binomial())
      else           stats::lm(fml,  data = dat2),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) return(NULL)

    sm    <- summary(fit)
    cf    <- stats::coef(sm)
    if (!(pred %in% rownames(cf))) return(NULL)

    b     <- cf[pred, "Estimate"]
    se    <- cf[pred, "Std. Error"]
    p_col <- if (is_binary) "Pr(>|z|)" else "Pr(>|t|)"
    if (!p_col %in% colnames(cf)) return(NULL)
    p     <- cf[pred, p_col]

    ci <- try(
      if (is_binary) stats::confint.default(fit, level = conf_level)
      else           stats::confint(fit,          level = conf_level),
      silent = TRUE
    )
    lo <- if (!inherits(ci, "try-error") && pred %in% rownames(ci)) ci[pred, 1] else NA_real_
    up <- if (!inherits(ci, "try-error") && pred %in% rownames(ci)) ci[pred, 2] else NA_real_

    r2_val <- if (is_binary) NA_real_ else as.numeric(sm$r.squared)

    trait_cat <- if (!is.null(panel_i) && panel_i %in% names(trait_map)) {
      unname(trait_map[panel_i])
    } else if (out %in% names(trait_map)) {
      unname(trait_map[out])
    } else {
      NA_character_
    }

    model_type <- if (is_binary) "logistic" else "linear"

    .fmt_term <- function(v) {
      tr <- unname(trans_map[v])
      if (!is.na(tr) && length(tr) > 0 && tr != "none") paste0(v, " [", tr, "]") else v
    }

    model_string <- paste0(
      "[", model_type, "] ",
      .fmt_term(out), " ~ ",
      paste(vapply(rhs_terms, .fmt_term, character(1)), collapse = " + ")
    )

    data.frame(
      outcome        = out,         predictor    = pred,
      trait_category = trait_cat,   model_type   = model_type,
      model          = model_string, N           = as.integer(stats::nobs(fit)),
      beta           = as.numeric(b), R2         = r2_val,
      pvalue         = as.numeric(p), pvalue_fdr = NA_real_,
      ci_upper       = as.numeric(up), ci_lower  = as.numeric(lo),
      se             = as.numeric(se),
      stringsAsFactors = FALSE
    )
  }

  res <- vector("list", length(outcomes) * length(predictors))
  k   <- 0L
  for (out in outcomes) for (pred in predictors) {
    k        <- k + 1L
    res[[k]] <- fit_one(df, out, pred)
  }

  results <- dplyr::bind_rows(res) %>% dplyr::arrange(dplyr::desc(beta))
  if (nrow(results) == 0) return(empty_out())

  results
}

summarise_population_characteristics_wide <- function(df,
                                                       binary_vars,
                                                       continuous_vars,
                                                       sex_col = "sex",
                                                       age_col = "age",
                                                       id_col  = "eid") {
  stopifnot(is.data.frame(df))

  vars_needed  <- unique(c(binary_vars, continuous_vars, sex_col, age_col, id_col))
  missing_vars <- setdiff(vars_needed, names(df))
  if (length(missing_vars) > 0)
    stop("Missing columns in df: ", paste(missing_vars, collapse = ", "))

  summarise_continuous <- function(v) {
    x <- df[[v]]
    data.frame(
      variable = v, type = "continuous",
      n_nonmissing = sum(!is.na(x) & !is.nan(x)),
      mean   = mean(x, na.rm = TRUE),
      median = stats::median(x, na.rm = TRUE),
      iqr    = stats::IQR(x, na.rm = TRUE),
      min    = min(x, na.rm = TRUE),
      max    = max(x, na.rm = TRUE),
      n_0 = NA, n_1 = NA, count = NA, pct = NA,
      stringsAsFactors = FALSE
    )
  }

  summarise_binary <- function(v) {
    x       <- df[[v]]
    x_valid <- x[!is.na(x) & !is.nan(x)]
    data.frame(
      variable = v, type = "binary",
      n_nonmissing = length(x_valid),
      mean = NA, median = NA, iqr = NA, min = NA, max = NA,
      n_0  = sum(x_valid == 0, na.rm = TRUE),
      n_1  = sum(x_valid == 1, na.rm = TRUE),
      count = NA, pct = NA,
      stringsAsFactors = FALSE
    )
  }

  sex_tab <- df %>%
    dplyr::distinct(.data[[id_col]], .data[[sex_col]]) %>%
    dplyr::filter(!is.na(.data[[sex_col]])) %>%
    dplyr::count(.data[[sex_col]], name = "count") %>%
    dplyr::mutate(
      pct          = 100 * count / sum(count),
      variable     = paste0("sex_", .data[[sex_col]]),
      type         = "sex_distribution",
      n_nonmissing = sum(count),
      mean = NA, median = NA, iqr = NA, min = NA, max = NA,
      n_0 = NA, n_1 = NA
    ) %>%
    dplyr::select(variable, type, n_nonmissing, mean, median, iqr, min, max, n_0, n_1, count, pct)

  age_vals <- df %>%
    dplyr::distinct(.data[[id_col]], .data[[age_col]]) %>%
    dplyr::pull(.data[[age_col]])

  age_tab <- data.frame(
    variable = age_col, type = "age_distribution",
    n_nonmissing = sum(!is.na(age_vals) & !is.nan(age_vals)),
    mean   = mean(age_vals, na.rm = TRUE),
    median = stats::median(age_vals, na.rm = TRUE),
    iqr    = stats::IQR(age_vals, na.rm = TRUE),
    min    = min(age_vals, na.rm = TRUE),
    max    = max(age_vals, na.rm = TRUE),
    n_0 = NA, n_1 = NA, count = NA, pct = NA,
    stringsAsFactors = FALSE
  )

  dplyr::bind_rows(
    lapply(continuous_vars, summarise_continuous),
    lapply(binary_vars,     summarise_binary),
    sex_tab,
    age_tab
  )
}

extract_generic_formula <- function(model_vec) {
  model_vec <- unique(stats::na.omit(model_vec))
  if (length(model_vec) == 0) return(NULL)

  x     <- sub("^\\[.*?\\]\\s*", "", model_vec[1])
  rhs   <- sub("^.*?~\\s*", "", x)
  terms <- trimws(strsplit(rhs, "\\s*\\+\\s*")[[1]])

  terms <- ifelse(grepl("^sample_age_", terms), "sample age", terms)
  terms[terms == "z_fasting_time"] <- "fasting time sds"

  pc_terms <- grep("^pc[0-9]+$", terms, value = TRUE)
  if (length(pc_terms) > 1) {
    pc_nums <- sort(as.integer(sub("^pc", "", pc_terms)))
    if (identical(pc_nums, seq(min(pc_nums), max(pc_nums)))) {
      terms <- terms[!terms %in% pc_terms]
      terms <- c(terms, paste0("pc", min(pc_nums), "-pc", max(pc_nums)))
    }
  }

  paste("model: outcome ~", paste(unique(terms), collapse = " + "))
}

save_figure <- function(p, outputpath, fname, width, height, dpi) {
  if (!dir.exists(outputpath)) dir.create(outputpath, recursive = TRUE)
  out_file <- file.path(outputpath, fname)
  ggplot2::ggsave(filename = out_file, plot = p, width = width / dpi, height = height / dpi, dpi = dpi)
  message("saved: ", out_file)
  invisible(out_file)
}

convert_highest <- function(x, qual_map) {
  if (is.na(x) || x == "") return(NA_real_)
  parts <- trimws(strsplit(x, "\\|")[[1]])
  vals  <- unname(qual_map[parts])
  vals  <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NA_real_)
  max(vals)
}

make_forest_plot <- function(model_results,
                             var_map      = setNames(character(0), character(0)),
                             trait_map    = setNames(character(0), character(0)),
                             facet_map    = setNames(character(0), character(0)),
                             trait_colors = character(0),
                             cohort       = "",
                             datatype     = "",
                             plot_type    = "forest",
                             outputpath   = "", 
                             p_nominal    = 0.05,
                             fdr_cutoff   = 0.05,
                             facet_by_predictor = FALSE,
                             caption      = NULL,
                             dpi          = 300,
                             width        = 3200,
                             height       = 3200,
                             save         = TRUE,
                             ext          = "png") {

  stopifnot(is.data.frame(model_results))

  needed <- c("outcome", "predictor", "trait_category", "beta", "ci_lower", "ci_upper", "pvalue", "pvalue_fdr", "model", "N")
  miss <- setdiff(needed, names(model_results))
  if (length(miss) > 0) stop("model_results missing: ", paste(miss, collapse = ", "))

  df <- model_results %>%
    dplyr::mutate(
      trait_category = dplyr::case_when(
        outcome %in% names(trait_map) ~ unname(trait_map[outcome]),
        TRUE ~ trait_category
      ),
      facet_group = dplyr::case_when(
        trait_category %in% names(facet_map) ~ unname(facet_map[trait_category]),
        TRUE ~ "other"
      ),
      outcome   = dplyr::if_else(outcome   %in% names(var_map), unname(var_map[outcome]),   outcome),
      predictor = dplyr::if_else(predictor %in% names(var_map), unname(var_map[predictor]), predictor)
    ) %>%
    dplyr::arrange(facet_group, beta) %>%
    dplyr::mutate(
      predictor      = factor(predictor, levels = unique(predictor)),
      outcome        = factor(outcome,   levels = unique(outcome)),
      signif_cat = dplyr::case_when(
        .data$pvalue_fdr < fdr_cutoff ~ "FDR p-value < 0.05",
        .data$pvalue     < p_nominal  ~ "p-value < 0.05",
        TRUE                          ~ "not significant"
      ),
      signif_cat      = factor(signif_cat, levels = c("FDR p-value < 0.05", "p-value < 0.05", "not significant")),
      facet_group    = factor(facet_group,    levels = unique(facet_group)),
      trait_category = factor(trait_category, levels = unique(trait_category))
    )

  present_trait_categories <- intersect(names(trait_colors), unique(as.character(stats::na.omit(df$trait_category))))
  df$trait_category <- factor(df$trait_category, levels = present_trait_categories)

  xintercept_val <- 0
  xlab_val <- "Effect size / log(OR) ± 95% CI"

  cap <- if (!is.null(caption)) caption else extract_generic_formula(model_results$model)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = beta, y = outcome)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lower, xmax = ci_upper, color = trait_category), width = 0.15, linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(color = trait_category, size = signif_cat)) +
    ggplot2::geom_vline(xintercept = xintercept_val, linetype = "dashed", color = "grey60") +
    ggplot2::scale_size_manual(name = "Significance", values = c("FDR p-value < 0.05" = 4, "p-value < 0.05" = 2, "not significant" = 1), drop = FALSE) +
    ggplot2::scale_color_manual(name = "Phenotype category", values = trait_colors, drop = FALSE) +
    ggplot2::labs(x = xlab_val, y = "Outcome") +
    make_theme() +
    ggplot2::facet_grid(facet_group ~ ., scales = "free_y", space = "free_y")

  if (facet_by_predictor)
    p <- p + ggplot2::facet_grid(facet_group ~ predictor, scales = "free_y", space = "free_y")

  if (isTRUE(save)) {
    fname <- make_fig_filename(datatype = datatype, plot_type = plot_type, ext = ext)
    save_figure(p, outputpath, fname, width, height, dpi)
  }

  p
}

plot_volcano_prs_proteins <- function(
    results_df,
    prs_name_map,
    cov_formula  = NULL,
    top_n        = 5,
    p_nominal    = 0.05,
    fdr_cutoff   = 0.05,
    x_breaks     = NULL,
    x_step       = 0.1,
    size_map     = c("FDR p-value < 0.05" = 4, "p-value < 0.05" = 2, "not significant" = 1),
    color_map    = PANEL_COLORS,
    datatype     = "",
    plot_type    = "volcano",
    outputpath   = "",
    dpi          = 300,
    width        = 4800,
    height       = 1200,
    save         = TRUE,
    ext          = "png") {

  stopifnot(is.data.frame(results_df))
  req     <- c("beta", "pvalue", "pvalue_fdr", "outcome", "predictor")
  missing <- setdiff(req, names(results_df))
  if (length(missing) > 0) stop("results_df missing: ", paste(missing, collapse = ", "))

  plot_df <- results_df %>%
    dplyr::mutate(neg_log10_p = -log10(.data$pvalue), prs_name = .data$predictor) %>%
    dplyr::group_by(.data$prs_name) %>%
    dplyr::arrange(.data$pvalue, .by_group = TRUE) %>%
    dplyr::mutate(label_top = ifelse(dplyr::row_number() <= top_n, .data$outcome, NA_character_)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      signif_cat = dplyr::case_when(
        .data$pvalue_fdr < fdr_cutoff ~ "FDR p-value < 0.05",
        .data$pvalue     < p_nominal  ~ "p-value < 0.05",
        TRUE                          ~ "not significant"
      ),
      signif_cat      = factor(signif_cat, levels = c("FDR p-value < 0.05", "p-value < 0.05", "not significant")),
      predictor_label = dplyr::recode(.data$prs_name, !!!prs_name_map)
    ) %>%
    dplyr::mutate(
    trait_category = dplyr::case_when(
      .data$signif_cat == "not significant" ~ "not significant",
      TRUE ~ .data$trait_category
    )
    )


  xvals <- plot_df$beta[is.finite(plot_df$beta)]
  if (length(xvals) == 0) stop("No finite beta values for volcano plot.")

  if (is.null(x_breaks)) {
    x_breaks <- seq(
      floor(min(xvals)   / x_step) * x_step,
      ceiling(max(xvals) / x_step) * x_step,
      by = x_step
    )
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$beta, y = .data$neg_log10_p, color = .data$trait_category, size = .data$signif_cat)) +
    ggplot2::geom_hline(yintercept = -log10(p_nominal), linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_manual(values = color_map, name = "Olink panel") +
    ggplot2::scale_size_manual(values  = size_map,  name = "Significance") +
    ggrepel::geom_label_repel(ggplot2::aes(label = .data$label_top), fill = "white", box.padding = 0.5, size = 3, na.rm = TRUE, show.legend = FALSE) +
    ggplot2::labs(x = "Effect size (PRS)", y = "-log10(p-value)") +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = function(x) sprintf("%.1f", x)) +
    make_theme() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = -90, hjust = 0, vjust = 0), legend.position = "right", plot.caption = ggplot2::element_text(size = 8), aspect.ratio = 1)

  if (isTRUE(save)) {
    fname <- make_fig_filename(datatype = datatype, plot_type = plot_type, ext = ext)
    save_figure(p, outputpath, fname, width, height, dpi)
  }

  p
}

plot_binary_by_prs_decile <- function(data,
                                      var,
                                      prs_var = "prs",
                                      n = 10,
                                      position_type = c("stack", "fill", "dodge"),
                                      output_dir = NULL,
                                      datatype = NULL,
                                      plot_type = NULL,
                                      width = 2400,
                                      height = 2400,
                                      dpi = 300) {
  
  position_type <- match.arg(position_type)
  
  required_cols <- c(prs_var, var)
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  df_plot <- data %>%
    dplyr::select(
      prs_value = dplyr::all_of(prs_var),
      binary_value = dplyr::all_of(var)
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      prs_decile = dplyr::ntile(prs_value, n),
      prs_decile = factor(
        prs_decile,
        levels = 1:n,
        labels = paste0("D", 1:n)
      ),
      binary_label = dplyr::case_when(
        binary_value == 0 ~ "No",
        binary_value == 1 ~ "Yes",
        TRUE ~ NA_character_
      ),
      binary_label = factor(binary_label, levels = c("No", "Yes"))
    ) %>%
    tidyr::drop_na(binary_label)
  
  count_df <- df_plot %>%
    dplyr::count(prs_decile, binary_label)
  
  legend_title <- gsub("_", " ", toupper(var))
  
  p <- ggplot2::ggplot(
    count_df,
    ggplot2::aes(
      x = prs_decile,
      y = n,
      fill = binary_label
    )
  ) +
    ggplot2::geom_col(
      position = position_type,
      width = 0.7
    ) +
    ggplot2::scale_fill_manual(
      values = c("No" = "#BDBDBD", "Yes" = "#E78AC3"),
      name = legend_title
    ) +
    ggplot2::labs(
      x = "PRS decile",
      y = if (position_type == "fill") "Proportion" else "Number of individuals"
    ) +
    make_theme() +
    ggplot2::theme(
      aspect.ratio = 1
    )
  
  if (position_type == "fill") {
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent)
  }
  
  if (!is.null(output_dir) && !is.null(datatype) && !is.null(plot_type)) {
  
  fname <- make_fig_filename(
    datatype = datatype,
    plot_type = plot_type,
    ext = "png"
  )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  file_out <- file.path(output_dir, fname)
  
  ggplot2::ggsave(
    filename = file_out,
    plot = p,
    width = width/dpi,
    height = height/dpi,
    dpi = dpi
  )
}
  
  return(p)
}

plot_protein_by_prs_decile <- function(data,
                                       protein,
                                       prs_var = "prs",
                                       id_var = "eid",
                                       output_dir = NULL,
                                       datatype = NULL,
                                       plot_type = NULL,
                                       width = 2400,
                                       height = 2400,
                                       dpi = 300,
                                       LOWER = -3,
                                       UPPER = 3,
                                       n = 10) {
  
  required_cols <- c(id_var, protein, prs_var)
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  df <- data %>%
    dplyr::select(
      dplyr::all_of(id_var),
      protein_value = dplyr::all_of(protein),
      prs_value = dplyr::all_of(prs_var)
    ) %>%
    tidyr::drop_na()
  
  df_plot <- df %>%
    dplyr::mutate(
      prs_decile = dplyr::ntile(prs_value, n),
      prs_decile = factor(
        prs_decile,
        levels = 1:n,
        labels = paste0("D", 1:n)
      )
    )
  
  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = prs_decile,
      y = protein_value
    )
  ) +
    ggplot2::geom_jitter(
      width = 0.08,
      height = 0,
      size = 0.2,
      alpha = 0.02,
      color = "black"
    ) +
    ggplot2::stat_summary(
      fun.data = function(x) {
        data.frame(
          y = median(x),
          ymin = quantile(x, 0.25),
          ymax = quantile(x, 0.75)
        )
      },
      geom = "pointrange",
      linewidth = 0.6,
      color = "black"
    ) +
    ggplot2::labs(
      x = "PRS decile",
      y = "Median \u00B1 IQR",
      subtitle = paste0("Protein: ", protein),
    ) +
    make_theme() +
    ggplot2::coord_cartesian(ylim = c(LOWER, UPPER)) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      aspect.ratio = 1
    )
  
  if (!is.null(output_dir) && !is.null(datatype) && !is.null(plot_type)) {
  
  fname <- make_fig_filename(
    datatype = datatype,
    plot_type = plot_type,
    ext = "png"
  )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  file_out <- file.path(output_dir, fname)
  
  ggplot2::ggsave(
    filename = file_out,
    plot = p,
    width = width/dpi,
    height = height/dpi,
    dpi = dpi
  )
}
  
  return(p)
}
