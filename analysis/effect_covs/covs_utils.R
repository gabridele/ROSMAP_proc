# Shared helpers for analyses in analysis/effect_covs/.
#
# Scientific model formulas intentionally remain in the individual analysis
# scripts. These helpers centralize repeated data preparation, FD filtering,
# significance labels, pairwise-contrast annotations, and console formatting.
#
# Expected source order in entry-point scripts:
#   1. analysis/paths.R
#   2. analysis/effect_covs/covs_utils.R

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggeffects)
library(lme4)
library(lmerTest)
library(lubridate)
library(emmeans)
library(ggpubr)
library(readxl)
library(grDevices)
library(visreg)

emm_options(pbkrtest.limit = 10000)
emm_options(lmerTest.limit = 50000)

fd_threshold <- FD_THRESHOLD

EC_MSEX_LEVELS <- c(
  "female",
  "male"
)

EC_SITE_LEVELS <- c(
  "BNK",
  "UC",
  "MG",
  "RIRC"
)

EC_DX_LEVELS <- c(
  "NCI",
  "MCI",
  "AD",
  "other"
)

EC_EYES_LEVELS <- c(
  "closed",
  "open"
)

EC_SYN_LEVELS <- c(
  "not SyN",
  "SyN"
)

target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

network_colors <- c(
  "Vis" = "#9B59B6",
  "SomMot" = "#6C8EBF",
  "Default" = "#D36B78",
  "Limbic" = "#C9D39A",
  "DorsAttn" = "#3C8D2F",
  "SalVentAttn" = "#C84CCF",
  "Cont" = "#E5B53A"
)

# Between-network columns
target_combos <- c(
  "Cont_to_Default", "Cont_to_DorsAttn",
  "Cont_to_Limbic", "Cont_to_SalVentAttn", "Cont_to_SomMot",
  "Cont_to_Vis", "Default_to_DorsAttn", "Default_to_Limbic",
  "Default_to_SalVentAttn", "Default_to_SomMot", "Default_to_Vis",
  "DorsAttn_to_Limbic", "DorsAttn_to_SalVentAttn", "DorsAttn_to_SomMot",
  "DorsAttn_to_Vis", "Limbic_to_SalVentAttn", "Limbic_to_SomMot",
  "Limbic_to_Vis", "SalVentAttn_to_SomMot", "SalVentAttn_to_Vis",
  "SomMot_to_Vis"
)

between_network_colors <- c(
  "Cont_to_Default" = "#DC9059",
  "Cont_to_DorsAttn" = "#91A135",
  "Cont_to_Limbic" = "#D7C46A",
  "Cont_to_SalVentAttn" = "#D78185",
  "Cont_to_SomMot" = "#A9A27D",
  "Cont_to_Vis" = "#C08778",
  "Default_to_DorsAttn" = "#887C54",
  "Default_to_Limbic" = "#CE9F89",
  "Default_to_SalVentAttn" = "#CE5CA3",
  "Default_to_SomMot" = "#A07D9C",
  "Default_to_Vis" = "#B76297",
  "DorsAttn_to_Limbic" = "#83B065",
  "DorsAttn_to_SalVentAttn" = "#826D7F",
  "DorsAttn_to_SomMot" = "#548E77",
  "DorsAttn_to_Vis" = "#6B7373",
  "Limbic_to_SalVentAttn" = "#C990B4",
  "Limbic_to_SomMot" = "#9BB1AD",
  "Limbic_to_Vis" = "#B296A8",
  "SalVentAttn_to_SomMot" = "#9A6DC7",
  "SalVentAttn_to_Vis" = "#B252C3",
  "SomMot_to_Vis" = "#8474BB"
)

covariates_to_run <- c(
  "msex",
  "site",
  "eyes",
  "syn_bin",
  "dcfdx"
)

ec_factor_with_levels <- function(
    x,
    levels,
    variable_name
) {

  observed <- unique(
    stats::na.omit(
      as.character(x)
    )
  )

  unexpected <- setdiff(
    observed,
    levels
  )

  if (length(unexpected) > 0) {
    stop(
      "Unexpected level(s) in ",
      variable_name,
      ": ",
      paste(
        unexpected,
        collapse = ", "
      )
    )
  }

  factor(
    x,
    levels = levels
  )
}

# Return conventional significance stars for a vector of p/q values.
ec_sig_from_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

# Add a numeric session index parsed from BIDS-style session labels.
ec_add_session_number <- function(data, session_col = "ses_id") {
  if (!session_col %in% names(data)) {
    stop("Session column not found: ", session_col)
  }

  data$ses_num <- as.numeric(
    stringr::str_extract(as.character(data[[session_col]]), "\\d+")
  )
  data
}

# Keep the earliest numeric session for each participant.
ec_select_baseline <- function(
    data,
    subject_col = "sub_id",
    session_col = "ses_id"
) {
  if (!"ses_num" %in% names(data)) {
    data <- ec_add_session_number(data, session_col = session_col)
  }

  data %>%
    dplyr::group_by(.data[[subject_col]]) %>%
    dplyr::arrange(ses_num, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

# Apply the prespecified publication motion-exclusion rule.
# FD_THRESHOLD is defined once in analysis/paths.R.
ec_apply_fd_filter <- function(
    data,
    threshold = FD_THRESHOLD,
    fd_col = "mean_FD"
) {
  if (!fd_col %in% names(data)) {
    stop("FD column not found: ", fd_col)
  }

  data %>%
    dplyr::filter(
      !is.na(.data[[fd_col]]),
      .data[[fd_col]] < threshold
    ) %>%
    droplevels()
}

# Convert model variables to stable types and apply the
# prespecified factor/reference-level ordering.
ec_prepare_model_variables <- function(
    data,
    longitudinal = FALSE
) {

  if (
    longitudinal &&
    "sub_id" %in% names(data)
  ) {
    data$sub_id <- factor(
      data$sub_id
    )
  }

  if ("mean_FD" %in% names(data)) {
    data$mean_FD <- as.numeric(
      data$mean_FD
    )
  }

  if ("age_scandate" %in% names(data)) {
    data$age_scandate <- as.numeric(
      data$age_scandate
    )
  }


  if ("msex" %in% names(data)) {

    data$msex <- ec_factor_with_levels(
      data$msex,
      EC_MSEX_LEVELS,
      "msex"
    )
  }


  if ("site" %in% names(data)) {

    data$site <- ec_factor_with_levels(
      data$site,
      EC_SITE_LEVELS,
      "site"
    )
  }


  if ("dcfdx" %in% names(data)) {

    data$dcfdx <- ec_factor_with_levels(
      data$dcfdx,
      EC_DX_LEVELS,
      "dcfdx"
    )
  }


  if ("syn_bin" %in% names(data)) {

    data$syn_bin <- ec_factor_with_levels(
      data$syn_bin,
      EC_SYN_LEVELS,
      "syn_bin"
    )
  }


  if ("eyes" %in% names(data)) {

    data$eyes <- ec_factor_with_levels(
      data$eyes,
      EC_EYES_LEVELS,
      "eyes"
    )
  }


  data
}

# Pivot connectivity measures to long format and apply one explicit
# complete-case definition. This replaces several subtly different make_long()
# copies. required_vars should contain every variable used by the model.
ec_make_long <- function(
    data,
    measure_cols,
    measure_name,
    value_name,
    required_vars,
    verbose = FALSE,
    site_levels = NULL
) {
  missing_measure_cols <- setdiff(measure_cols, names(data))
  if (length(missing_measure_cols)) {
    stop(
      "Missing connectivity columns: ",
      paste(missing_measure_cols, collapse = ", ")
    )
  }

  long_data <- tidyr::pivot_longer(
    data,
    cols = dplyr::all_of(measure_cols),
    names_to = measure_name,
    values_to = value_name
  )

  long_data[[measure_name]] <- factor(
    long_data[[measure_name]],
    levels = measure_cols
  )

  if (!is.null(site_levels) && "site" %in% names(long_data)) {
    long_data$site <- factor(long_data$site, levels = site_levels)
  }

  required_vars <- unique(c(value_name, required_vars))
  missing_required <- setdiff(required_vars, names(long_data))
  if (length(missing_required)) {
    stop(
      "Missing variables required by the model: ",
      paste(missing_required, collapse = ", ")
    )
  }

  if (verbose) {
    cat("\nRows after pivot:", nrow(long_data), "\n")
    missing_counts <- vapply(
      required_vars,
      function(column) sum(is.na(long_data[[column]])),
      integer(1)
    )
    print(data.frame(variable = names(missing_counts), missing = missing_counts))
  }

  keep <- stats::complete.cases(
    long_data[, required_vars, drop = FALSE]
  )
  long_data <- droplevels(long_data[keep, , drop = FALSE])

  if (verbose) {
    cat("\nRows after model filter:", nrow(long_data), "\n")
  }

  long_data
}

# Create alternating significance-bracket positions above/below each facet.
# This supports both within-network (`network`) and between-network
# (`network_combo`) plots.
ec_make_pairwise_brackets <- function(
    plot_data,
    pairwise_results,
    covariate,
    measure_col,
    value_col,
    group_col = covariate,
    p_col = "q_across",
    sig_col = "sig_across"
) {
  if (!all(c(measure_col, value_col, group_col) %in% names(plot_data))) {
    stop("Plot data are missing required bracket columns.")
  }

  x_raw <- plot_data[[group_col]]
  x_levels <- if (is.factor(x_raw)) {
    levels(x_raw)
  } else {
    unique(as.character(x_raw))
  }

  y_positions <- plot_data %>%
    dplyr::group_by(.data[[measure_col]]) %>%
    dplyr::summarise(
      y_min = min(.data[[value_col]], na.rm = TRUE),
      y_max = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_range = y_max - y_min,
      y_range = ifelse(!is.finite(y_range) | y_range == 0, 1, y_range)
    )

  annotations <- pairwise_results %>%
    dplyr::filter(.data[[p_col]] < 0.05) %>%
    tidyr::separate(
      contrast,
      into = c("group1", "group2"),
      sep = " - ",
      remove = FALSE
    ) %>%
    dplyr::mutate(
      group1_clean = stringr::str_remove(group1, paste0("^", covariate)),
      group2_clean = stringr::str_remove(group2, paste0("^", covariate)),
      group1 = ifelse(group1 %in% x_levels, group1, group1_clean),
      group2 = ifelse(group2 %in% x_levels, group2, group2_clean)
    ) %>%
    dplyr::select(-group1_clean, -group2_clean)

  annotations <- dplyr::left_join(
    annotations,
    y_positions,
    by = measure_col
  ) %>%
    dplyr::group_by(.data[[measure_col]]) %>%
    dplyr::arrange(.data[[p_col]], .by_group = TRUE) %>%
    dplyr::mutate(
      bracket_number = dplyr::row_number(),
      bracket_rank = ceiling(bracket_number / 2),
      y.position = ifelse(
        bracket_number %% 2 == 1,
        y_max + bracket_rank * 0.08 * y_range,
        y_min - bracket_rank * 0.08 * y_range
      ),
      label = .data[[sig_col]]
    ) %>%
    dplyr::ungroup()

  annotations
}

# Print pairwise statistics without rounding the underlying result object.
ec_print_pairwise_table <- function(
    pairwise_results,
    title
) {

  cat("\n============================================================\n")
  cat(title, "\n")
  cat("============================================================\n")

  pairwise_results %>%
    dplyr::mutate(
      estimate = round(
        estimate,
        4
      ),

      se = round(
        se,
        4
      ),

      statistic = round(
        as.numeric(statistic),
        3
      ),

      p_raw = signif(
        p_raw,
        3
      ),

      p_tukey = signif(
        p_tukey,
        3
      ),

      q_across = signif(
        q_across,
        3
      )
    ) %>%
    tibble::as_tibble() %>%
    print(
      n = Inf
    )
}

# Print FD-model statistics without modifying values subsequently written to CSV.
ec_print_fd_model_table <- function(model_results, title, measure_col = NULL) {
  cat("\n============================================================\n")
  cat(title, "\n")
  cat("============================================================\n")

  display <- model_results

  if (is.null(measure_col)) {
    measure_col <- intersect(c("network", "network_combo"), names(display))[1]
  }

  preferred <- c(
    measure_col,
    "beta_adjusted", "t_val_adjusted",
    "p_adjusted", "q_adjusted", "sig_adjusted"
  )
  preferred <- preferred[!is.na(preferred) & preferred %in% names(display)]
  if (length(preferred)) display <- display[, preferred, drop = FALSE]

  if ("beta_adjusted" %in% names(display)) {
    display$beta_adjusted <- round(display$beta_adjusted, 4)
  }
  if ("t_val_adjusted" %in% names(display)) {
    display$t_val_adjusted <- round(display$t_val_adjusted, 3)
  }
  for (column in intersect(c("p_adjusted", "q_adjusted"), names(display))) {
    display[[column]] <- signif(display[[column]], 3)
  }

  print(tibble::as_tibble(display), n = Inf)
}

# Fit the same scientific model independently for each connectivity measure.
# The formula is always supplied by the entry-point script so the model remains
# visible next to the analysis it represents.
ec_fit_models_by_measure <- function(
    data_long,
    measure_levels,
    measure_col,
    model_formula,
    mixed = FALSE,
    optimizer = "bobyqa"
) {
  models <- purrr::map(
    measure_levels,
    function(measure) {
      df_measure <- data_long[
        data_long[[measure_col]] == measure,
        ,
        drop = FALSE
      ]
      df_measure <- droplevels(df_measure)

      if (mixed) {
        if (is.null(optimizer)) {
          lmerTest::lmer(
            model_formula,
            data = df_measure
          )
        } else {
          lmerTest::lmer(
            model_formula,
            data = df_measure,
            control = lme4::lmerControl(optimizer = optimizer)
          )
        }
      } else {
        stats::lm(model_formula, data = df_measure)
      }
    }
  )

  stats::setNames(models, measure_levels)
}

# Attach fitted values from an already-fitted model list. For mixed models the
# default deliberately preserves the historical behavior of fitted(), which
# includes estimated random effects. Set fixed_only = TRUE for population-level
# predictions instead.
ec_get_fitted_data <- function(
    data_long,
    models,
    measure_levels,
    measure_col,
    prediction_col = "predicted_conn",
    fixed_only = FALSE
) {
  out <- purrr::map_dfr(
    measure_levels,
    function(measure) {
      df_measure <- data_long[
        data_long[[measure_col]] == measure,
        ,
        drop = FALSE
      ]
      df_measure <- droplevels(df_measure)
      model <- models[[measure]]

      predicted <- if (
        fixed_only && inherits(model, "merMod")
      ) {
        stats::predict(model, newdata = df_measure, re.form = NA)
      } else {
        stats::fitted(model)
      }

      df_measure[[prediction_col]] <- as.numeric(predicted)
      df_measure
    }
  )

  out[[measure_col]] <- factor(out[[measure_col]], levels = measure_levels)
  out
}

# Run Tukey-adjusted estimated-marginal-mean contrasts from a fitted model list.
# p_adj is intentionally retained as the historical column name; it represents
# the emmeans/Tukey adjustment within each fitted connectivity outcome.
ec_pairwise_from_models <- function(
    models,
    covariate,
    measure_col,
    measure_levels
) {

  results <- purrr::imap_dfr(
    models,
    function(model, measure) {

      emm_args <- list(
        object = model,
        specs = covariate
      )

      if (inherits(model, "merMod")) {
        emm_args$lmer.df <- "satterthwaite"
      }


      emm <- do.call(
        emmeans::emmeans,
        emm_args
      )


      contrasts_raw <- as.data.frame(
        emmeans::pairs(
          emm,
          adjust = "none"
        )
      )


      contrasts_tukey <- as.data.frame(
        emmeans::pairs(
          emm,
          adjust = "tukey"
        )
      )


      stat_col <- intersect(
        c("t.ratio", "z.ratio"),
        names(contrasts_raw)
      )[1]


      if (is.na(stat_col)) {
        stop(
          "Could not identify test-statistic column for ",
          measure
        )
      }


      out <- tibble::tibble(
        measure = measure,
        covariate = covariate,
        contrast = contrasts_raw$contrast,

        estimate = contrasts_raw$estimate,
        se = contrasts_raw$SE,

        df = if (
          "df" %in% names(contrasts_raw)
        ) {
          contrasts_raw$df
        } else {
          NA_real_
        },

        statistic = contrasts_raw[[stat_col]],

        p_raw = contrasts_raw$p.value,
        p_tukey = contrasts_tukey$p.value
      )


      names(out)[1] <- measure_col

      out
    }
  )


  # BH-FDR across all connectivity outcomes and pairwise
  # contrasts for this covariate/analysis dataset.
  results$q_across <- stats::p.adjust(
    results$p_raw,
    method = "BH"
  )


  results$sig_tukey <- ec_sig_from_p(
    results$p_tukey
  )

  results$sig_across <- ec_sig_from_p(
    results$q_across
  )


  results[[measure_col]] <- factor(
    results[[measure_col]],
    levels = measure_levels
  )


  results
}

# Extract visreg partial residuals from an already-fitted mixed-model list.
ec_partial_residuals_from_models <- function(
    models,
    covariate,
    measure_col,
    measure_levels,
    value_col = "value"
) {
  out <- purrr::imap_dfr(
    models,
    function(model, measure) {
      v <- visreg::visreg(
        model,
        covariate,
        plot = FALSE,
        predict = list(re.form = NA)
      )

      residual_col <- intersect(
        c("visreg_res", "visregRes"),
        names(v$res)
      )[1]
      if (is.na(residual_col)) {
        stop(
          "Could not find partial residuals for ",
          covariate,
          " in ",
          measure
        )
      }

      result <- data.frame(
        measure = measure,
        group = v$res[[covariate]],
        value = v$res[[residual_col]],
        stringsAsFactors = FALSE
      )
      names(result)[1] <- measure_col
      names(result)[3] <- value_col
      tibble::as_tibble(result)
    }
  )

  out[[measure_col]] <- factor(out[[measure_col]], levels = measure_levels)
  out
}

# Extract the mean-FD coefficient from each fitted model and apply BH correction
# across the connectivity measures in that analysis family.
ec_fd_results_from_models <- function(
    models,
    measure_col,
    measure_levels
) {
  results <- purrr::imap_dfr(
    models,
    function(model, measure) {
      coefs <- summary(model)$coefficients
      if (!"mean_FD" %in% rownames(coefs)) {
        stop("mean_FD coefficient not found for ", measure)
      }

      out <- data.frame(
        measure = measure,
        beta_adjusted = coefs["mean_FD", "Estimate"],
        t_val_adjusted = coefs["mean_FD", "t value"],
        p_adjusted = coefs["mean_FD", "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
      names(out)[1] <- measure_col
      tibble::as_tibble(out)
    }
  )

  results$q_adjusted <- stats::p.adjust(results$p_adjusted, method = "BH")
  results$sig_adjusted <- ec_sig_from_p(results$q_adjusted)
  results$label <- sprintf(
    "β = %.3f %s\nq = %.3g",
    results$beta_adjusted,
    results$sig_adjusted,
    results$q_adjusted
  )
  results[[measure_col]] <- factor(
    results[[measure_col]],
    levels = measure_levels
  )
  results
}

# Generate ggeffects mean-FD predictions from a fitted model list.
ec_fd_predictions_from_models <- function(
    models,
    measure_col,
    measure_levels,
    fixed_only = FALSE
) {
  predictions <- purrr::imap_dfr(
    models,
    function(model, measure) {
      args <- list(
        model = model,
        terms = "mean_FD [all]"
      )
      if (fixed_only && inherits(model, "merMod")) {
        args$type <- "fixed"
      }

      out <- do.call(ggeffects::predict_response, args) %>%
        as.data.frame()
      out[[measure_col]] <- measure
      out
    }
  )

  predictions[[measure_col]] <- factor(
    predictions[[measure_col]],
    levels = measure_levels
  )
  predictions
}

# Common FD-effect plot for the compact baseline/between-network scripts.
# The longitudinal within-network analysis keeps its more specialized plotting
# function because it has publication-specific axis and panel formatting.
ec_plot_fd_effects <- function(
    data_long,
    predictions,
    model_results,
    title,
    measure_col,
    value_col,
    palette,
    y_label,
    label_box = TRUE,
    label_size = 2.5,
    strip_size = NULL,
    rotate_x = FALSE
) {
  label_pos <- data_long %>%
    dplyr::group_by(.data[[measure_col]]) %>%
    dplyr::summarise(
      x = max(mean_FD, na.rm = TRUE),
      y = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(model_results, by = measure_col)

  p <- ggplot2::ggplot(
    data_long,
    ggplot2::aes(
      x = mean_FD,
      y = .data[[value_col]],
      color = .data[[measure_col]]
    )
  ) +
    ggplot2::geom_point(alpha = 0.30, size = 1) +
    ggplot2::geom_ribbon(
      data = predictions,
      ggplot2::aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = .data[[measure_col]]
      ),
      inherit.aes = FALSE,
      alpha = if (label_box) 0.22 else 0.25
    ) +
    ggplot2::geom_line(
      data = predictions,
      ggplot2::aes(
        x = x,
        y = predicted,
        color = .data[[measure_col]]
      ),
      inherit.aes = FALSE,
      linewidth = 1
    )

  if (label_box) {
    p <- p + ggplot2::geom_label(
      data = label_pos,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = label_size,
      color = "black",
      fill = "white",
      alpha = 0.85,
      linewidth = 0
    )
  } else {
    p <- p + ggplot2::geom_text(
      data = label_pos,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = label_size,
      color = "black"
    )
  }

  facet_formula <- stats::as.formula(paste("~", measure_col))
  strip_text <- if (is.null(strip_size)) {
    ggplot2::element_text(face = "bold")
  } else {
    ggplot2::element_text(face = "bold", size = strip_size)
  }
  x_text <- if (rotate_x) {
    ggplot2::element_text(angle = 45, hjust = 1)
  } else {
    ggplot2::element_text()
  }

  p +
    ggplot2::facet_wrap(facet_formula, scales = "fixed") +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      strip.text = strip_text,
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = x_text
    ) +
    ggplot2::labs(
      title = title,
      x = "mean_FD",
      y = y_label
    )
}

# Generic violin/box/jitter distribution used by categorical covariate plots.
# It supports both model-fitted distributions and visreg partial residuals.
ec_plot_factor_distribution <- function(
    plot_data,
    covariate,
    pairwise_results,
    title,
    subtitle,
    y_label,
    measure_col,
    value_col,
    group_col,
    palette,
    y_limits = NULL,
    y_breaks = NULL,
    facet_ncol = NULL,
    facet_axes_all = FALSE,
    x_text_angle = 45,
    x_text_size = NULL,
    y_text_size = NULL,
    base_size = 13
) {
  pairwise_annot <- ec_make_pairwise_brackets(
    plot_data = plot_data,
    pairwise_results = pairwise_results,
    covariate = covariate,
    measure_col = measure_col,
    value_col = value_col,
    group_col = group_col
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[group_col]],
      y = .data[[value_col]],
      color = .data[[measure_col]],
      fill = .data[[measure_col]]
    )
  ) +
    ggplot2::geom_violin(
      alpha = 0.18,
      linewidth = 0.3,
      trim = FALSE
    ) +
    ggplot2::geom_boxplot(
      width = 0.30,
      alpha = 0.75,
      outlier.shape = NA,
      linewidth = 0.35
    ) +
    ggplot2::geom_jitter(
      width = 0.12,
      alpha = 0.09,
      size = 0.4
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(
      breaks = if (is.null(y_breaks)) ggplot2::waiver() else y_breaks,
      expand = ggplot2::expansion(mult = c(0.18, 0.18))
    ) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::theme_minimal(base_size = base_size)

  facet_formula <- stats::as.formula(paste("~", measure_col))
  if (facet_axes_all) {
    p <- p + ggplot2::facet_wrap(
      facet_formula,
      ncol = facet_ncol,
      scales = "fixed",
      axes = "all",
      axis.labels = "all"
    )
  } else {
    p <- p + ggplot2::facet_wrap(
      facet_formula,
      ncol = facet_ncol,
      scales = "fixed"
    )
  }

  x_theme <- ggplot2::element_text(
    angle = x_text_angle,
    hjust = if (x_text_angle == 0) 0.5 else 1,
    size = x_text_size
  )
  y_theme <- ggplot2::element_text(size = y_text_size)

  p <- p +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = x_theme,
      axis.text.y = y_theme
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = covariate,
      y = y_label
    )

  if (nrow(pairwise_annot) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      pairwise_annot,
      label = "label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0,
      size = 4,
      bracket.size = 0.35,
      hide.ns = TRUE
    )
  }

  p
}

# Save an editable vector pair with consistent dimensions.
ec_save_pdf_svg <- function(
    plot,
    out_dir,
    filename,
    width = 13,
    height = 7,
    cairo = FALSE
) {
  pdf_args <- list(
    filename = file.path(out_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height
  )
  if (cairo) pdf_args$device <- grDevices::cairo_pdf
  do.call(ggplot2::ggsave, pdf_args)

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(filename, ".svg")),
    plot = plot,
    device = svglite::svglite,
    width = width,
    height = height
  )
}
