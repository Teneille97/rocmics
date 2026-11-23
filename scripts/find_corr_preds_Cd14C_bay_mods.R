
# =============================================================================
# MALLE ERW: DIAGNOSTIC TABLES FOR CORRELATION-FILTERED PREDICTORS
#
# Purpose:
#   1. Reproduce the original Bayesian predictor preprocessing for every C
#      and d14C response.
#   2. Identify the predictor paired with each excluded predictor at the
#      exact iteration when it was removed.
#   3. Calculate correlations between each excluded predictor and ALL
#      predictors retained after the complete correlation-filtering process.
#   4. Save CSV tables. No Bayesian models are fitted.
#
# Outputs:
#   outputs/bayesian_C_d14C_predictions/tables/
#     correlation_exclusion_pairs_at_removal.csv
#     correlation_exclusion_vs_all_final_retained_predictors.csv
#     correlation_exclusion_status_by_response.csv
#
# Run from the RStudio project root.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# -----------------------------------------------------------------------------
# 0. Paths and settings
# -----------------------------------------------------------------------------

input_file <- here::here(
  "outputs", "integrated_geochemistry_for_bayesian.csv"
)

out_dir <- here::here(
  "outputs", "bayesian_C_d14C_predictions_withftir/tables"
)

table_dir <- file.path(out_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

min_predictor_observed_fraction <- 0.70
max_abs_predictor_correlation <- 0.95

# Predictors specifically being investigated.
target_dropped_predictors <- c(
  "Tessier_Acetic_K_mg_g",
  "FTIR_delta_carbonyl",
  "FTIR_delta_poorly_crystalline",
  "H_Oct2023"
)

# -----------------------------------------------------------------------------
# 1. Read and prepare the integrated dataset
# -----------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

integrated <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  na = c("", "NA", "#VALUE!")
)

if (!all(c("sample", "Treatment") %in% names(integrated))) {
  stop("The integrated CSV must contain 'sample' and 'Treatment'.")
}

normalise_treatment <- function(x) {
  z <- stringr::str_to_lower(stringr::str_squish(as.character(x)))
  z <- stringr::str_replace_all(z, "ü", "u")
  z <- stringr::str_replace_all(z, "[^a-z0-9]+", "")
  
  dplyr::case_when(
    stringr::str_detect(z, "^control$") ~ "Control",
    stringr::str_detect(z, "^lime") ~ "Lime",
    stringr::str_detect(z, "^bolsdorfer") ~ "Bolsdorfer",
    stringr::str_detect(z, "^huhnerberg") ~ "Huhnerberg",
    stringr::str_detect(z, "^eifelgold") ~ "Eifelgold",
    TRUE ~ NA_character_
  )
}

allowed_treatments <- c(
  "Control", "Lime_2", "Bolsdorfer_50",
  "Huhnerberg_50", "Eifelgold_50"
)

dat <- integrated %>%
  filter(Treatment %in% allowed_treatments) %>%
  mutate(
    sample = suppressWarnings(as.integer(as.character(sample))),
    Treatment_model = normalise_treatment(Treatment)
  ) %>%
  filter(!is.na(sample), !is.na(Treatment_model))

if (nrow(dat) == 0L) {
  stop("No observations remain after treatment filtering.")
}

# -----------------------------------------------------------------------------
# 2. Reconstruct the original candidate predictor list
# -----------------------------------------------------------------------------

explicit_geochem_predictors <- c(
  "Agg stab", "GMD_h2o", "GMD_oxalate", "MWD_h2o", "MWD_oxalate",
  "D_h2o", "D_oxalate", "Dh2o", "Doxalate",
  "Fe_ox_over_DCB", "Al_ox_over_DCB", "Mn_ox_over_DCB", "Si_ox_over_DCB",
  "Tessier_Acetic_Ca_mg_g", "Tessier_Acetic_K_mg_g",
  "Tessier_Acetic_Na_mg_g", "Tessier_Acetic_Mg_mg_g",
  "Tessier_Hydroxylamine_Fe_mg_g", "Tessier_Hydroxylamine_Al_mg_g",
  "Tessier_Hydroxylamine_Mn_mg_g", "Tessier_Hydroxylamine_Si_mg_g",
  "soil_pH", "soil_EC",
  "pH_Oct2023", "EC_Oct2023", "H_Oct2023", "Ca_Oct2023",
  "K_Oct2023", "Mg_Oct2023", "Na_Oct2023", "Al_Oct2023",
  "Fe_Oct2023", "Mn_Oct2023", "CEC_Oct2023",
  "Base_saturation_Oct2023"
)

numeric_cols <- names(dat)[
  vapply(dat, is.numeric, logical(1))
]

pattern_predictors <- numeric_cols[
  stringr::str_starts(numeric_cols, "Tessier_") |
    stringr::str_starts(numeric_cols, "CEC") |
    stringr::str_starts(numeric_cols, "FTIR_delta") |
    stringr::str_starts(numeric_cols, "Fe_ox") |
    stringr::str_starts(numeric_cols, "Al_ox") |
    stringr::str_starts(numeric_cols, "Mn_ox") |
    stringr::str_starts(numeric_cols, "Si_ox")
]

candidate_predictors <- unique(c(
  intersect(explicit_geochem_predictors, names(dat)),
  intersect(pattern_predictors, names(dat))
))

response_like <- names(dat)[
  stringr::str_detect(
    names(dat),
    "^C_(Bulk|MAOM|oPOM|fPOM)$"
  ) |
    stringr::str_detect(names(dat), "^d14C")
]

metadata_cols <- c(
  "sample", "Treatment", "Treatment_model", "Basalt_type_exp",
  "Dose_resp_exp", "Date", "year", "depth", "fraction", "SOM_fraction"
)

candidate_predictors <- setdiff(
  candidate_predictors,
  c(response_like, metadata_cols)
)

candidate_predictors <- candidate_predictors[
  !stringr::str_detect(
    candidate_predictors,
    regex(
      "dataset|dose_resp_exp|basalt_type_exp|treatment|^sample$",
      ignore_case = TRUE
    )
  )
]

candidate_predictors <- candidate_predictors[
  vapply(
    candidate_predictors,
    function(v) is.numeric(dat[[v]]),
    logical(1)
  )
]

if (length(candidate_predictors) == 0L) {
  stop("No candidate predictors were found.")
}

missing_targets <- setdiff(
  target_dropped_predictors,
  candidate_predictors
)

if (length(missing_targets)) {
  warning(
    "These target predictors are not in the reconstructed candidate list: ",
    paste(missing_targets, collapse = ", ")
  )
}

message("Candidate predictors: ", length(candidate_predictors))

# -----------------------------------------------------------------------------
# 3. Define responses
# -----------------------------------------------------------------------------

c_responses <- c(
  "C_Bulk", "C_MAOM", "C_oPOM", "C_fPOM"
)

d14c_responses <- c(
  "d14C_2026_Bulk", "d14C_2026_MAOM",
  "d14C_2026_oPOM", "d14C_2026_fPOM"
)

responses <- intersect(
  c(c_responses, d14c_responses),
  names(dat)
)

if (length(responses) == 0L) {
  stop("None of the expected C or d14C response columns were found.")
}

# -----------------------------------------------------------------------------
# 4. Reproduce the original correlation filtering and record every step
# -----------------------------------------------------------------------------

pair_results <- list()
retained_correlation_results <- list()
status_results <- list()

for (response in responses) {
  
  message("Processing: ", response)
  
  # Use exactly the observations available for this response.
  y_raw <- suppressWarnings(as.numeric(dat[[response]]))
  keep_y <- is.finite(y_raw)
  
  work <- dat[keep_y, , drop = FALSE]
  
  if (sum(keep_y) < 8L) {
    warning("Skipping ", response, ": fewer than 8 observations.")
    next
  }
  
  # Apply the original missingness and variance filters.
  eligible <- candidate_predictors[
    vapply(candidate_predictors, function(v) {
      x <- suppressWarnings(as.numeric(work[[v]]))
      
      observed_fraction <- mean(is.finite(x))
      sx <- stats::sd(x, na.rm = TRUE)
      
      is.finite(observed_fraction) &&
        observed_fraction >= min_predictor_observed_fraction &&
        is.finite(sx) &&
        sx > 0
    }, logical(1))
  ]
  
  if (length(eligible) == 0L) {
    warning("No eligible predictors for ", response)
    next
  }
  
  # Median-impute and standardise exactly as in the Bayesian script.
  X <- as.data.frame(lapply(eligible, function(v) {
    
    x <- suppressWarnings(as.numeric(work[[v]]))
    
    med <- stats::median(
      x[is.finite(x)],
      na.rm = TRUE
    )
    
    x[!is.finite(x)] <- med
    
    as.numeric(scale(x))
  }))
  
  names(X) <- eligible
  
  # Remove any columns that became non-finite or constant.
  valid <- vapply(
    X,
    function(x) {
      all(is.finite(x)) &&
        is.finite(stats::sd(x)) &&
        stats::sd(x) > 0
    },
    logical(1)
  )
  
  X <- X[, valid, drop = FALSE]
  
  if (ncol(X) < 2L) {
    warning("Fewer than two predictors remain for ", response)
    next
  }
  
  # Preserve the complete eligible predictor matrix so that we can later
  # calculate correlations between excluded variables and all final survivors.
  X_initial <- X
  
  # ---------------------------------------------------------------------------
  # Iteratively reproduce the original correlation-filtering algorithm.
  # ---------------------------------------------------------------------------
  
  dropped_correlated <- character(0)
  step_number <- 0L
  
  repeat {
    
    if (ncol(X) < 2L) break
    
    cm <- suppressWarnings(
      stats::cor(X, use = "pairwise.complete.obs")
    )
    
    diag(cm) <- 0
    
    max_cor <- max(abs(cm), na.rm = TRUE)
    
    if (!is.finite(max_cor) ||
        max_cor <= max_abs_predictor_correlation) {
      break
    }
    
    # Same tie-breaking and removal rule as the original model script.
    idx <- which(
      abs(cm) == max_cor,
      arr.ind = TRUE
    )[1, ]
    
    paired_predictor <- colnames(cm)[idx[1]]
    dropped_predictor <- colnames(cm)[idx[2]]
    
    signed_r <- cm[idx[1], idx[2]]
    
    step_number <- step_number + 1L
    
    dropped_correlated <- c(
      dropped_correlated,
      dropped_predictor
    )
    
    pair_results[[length(pair_results) + 1L]] <- tibble(
      response = response,
      response_type = ifelse(
        stringr::str_starts(response, "C_"),
        "C",
        "d14C"
      ),
      n_observations = nrow(work),
      exclusion_step = step_number,
      dropped_predictor = dropped_predictor,
      paired_predictor = paired_predictor,
      Pearson_r = signed_r,
      absolute_Pearson_r = abs(signed_r),
      correlation_direction = ifelse(
        signed_r > 0,
        "Positive",
        "Negative"
      )
    )
    
    # Remove the same predictor as the original Bayesian workflow.
    X <- X[
      ,
      setdiff(names(X), dropped_predictor),
      drop = FALSE
    ]
  }
  
  final_retained <- names(X)
  
  message(
    "  Final retained predictors: ", length(final_retained),
    "; correlation-filtered predictors: ", length(dropped_correlated)
  )
  
  # ---------------------------------------------------------------------------
  # 5. Correlate each target excluded variable with ALL final survivors.
  # ---------------------------------------------------------------------------
  
  for (target in target_dropped_predictors) {
    
    if (!target %in% names(X_initial)) {
      
      status_results[[length(status_results) + 1L]] <- tibble(
        response = response,
        response_type = ifelse(
          stringr::str_starts(response, "C_"),
          "C",
          "d14C"
        ),
        predictor = target,
        status = "Not eligible after missingness/variance filtering",
        n_observations = nrow(work)
      )
      
      next
    }
    
    was_dropped <- target %in% dropped_correlated
    
    status_results[[length(status_results) + 1L]] <- tibble(
      response = response,
      response_type = ifelse(
        stringr::str_starts(response, "C_"),
        "C",
        "d14C"
      ),
      predictor = target,
      status = ifelse(
        was_dropped,
        "Removed by correlation filtering",
        "Not removed by correlation filtering"
      ),
      n_observations = nrow(work)
    )
    
    if (!was_dropped) next
    
    # Calculate correlations with every predictor that survived the complete
    # iterative filtering process.
    for (retained in final_retained) {
      
      x_target <- X_initial[[target]]
      x_retained <- X_initial[[retained]]
      
      pearson_r <- suppressWarnings(
        stats::cor(
          x_target,
          x_retained,
          method = "pearson",
          use = "complete.obs"
        )
      )
      
      spearman_rho <- suppressWarnings(
        stats::cor(
          x_target,
          x_retained,
          method = "spearman",
          use = "complete.obs"
        )
      )
      
      retained_correlation_results[[length(retained_correlation_results) + 1L]] <- tibble(
        response = response,
        response_type = ifelse(
          stringr::str_starts(response, "C_"),
          "C",
          "d14C"
        ),
        n_observations = nrow(work),
        excluded_predictor = target,
        retained_predictor = retained,
        Pearson_r = pearson_r,
        absolute_Pearson_r = abs(pearson_r),
        Spearman_rho = spearman_rho,
        absolute_Spearman_rho = abs(spearman_rho)
      )
    }
  }
}

# -----------------------------------------------------------------------------
# 6. Combine and save tables
# -----------------------------------------------------------------------------

if (length(pair_results) > 0L) {
  
  pairs_table <- bind_rows(pair_results) %>%
    arrange(
      dropped_predictor,
      response_type,
      response,
      exclusion_step
    )
  
  # Check whether the paired predictor was itself subsequently removed.
  pairs_table <- pairs_table %>%
    rowwise() %>%
    mutate(
      paired_predictor_later_dropped = any(
        pair_results %>%
          bind_rows() %>%
          filter(
            response == .data$response,
            dropped_predictor == .data$paired_predictor
          ) %>%
          nrow() > 0
      )
    ) %>%
    ungroup()
  
  readr::write_csv(
    pairs_table,
    file.path(
      table_dir,
      "correlation_exclusion_pairs_at_removal.csv"
    )
  )
  
} else {
  
  warning("No correlation exclusions were recorded.")
}

if (length(retained_correlation_results) > 0L) {
  
  retained_correlations <- bind_rows(
    retained_correlation_results
  ) %>%
    arrange(
      excluded_predictor,
      response_type,
      response,
      desc(absolute_Pearson_r)
    )
  
  readr::write_csv(
    retained_correlations,
    file.path(
      table_dir,
      "correlation_exclusion_vs_all_final_retained_predictors.csv"
    )
  )
  
} else {
  
  warning("No correlations with final retained predictors were calculated.")
}

if (length(status_results) > 0L) {
  
  status_table <- bind_rows(status_results) %>%
    arrange(predictor, response_type, response)
  
  readr::write_csv(
    status_table,
    file.path(
      table_dir,
      "correlation_exclusion_status_by_response.csv"
    )
  )
}

# -----------------------------------------------------------------------------
# 7. Create a compact summary: strongest final retained correlate per target
# -----------------------------------------------------------------------------

if (length(retained_correlation_results) > 0L) {
  
  strongest_correlates <- bind_rows(
    retained_correlation_results
  ) %>%
    group_by(response, response_type, excluded_predictor) %>%
    slice_max(
      order_by = absolute_Pearson_r,
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    arrange(
      excluded_predictor,
      response_type,
      response
    )
  
  readr::write_csv(
    strongest_correlates,
    file.path(
      table_dir,
      "strongest_retained_correlate_for_each_excluded_predictor.csv"
    )
  )
}

# -----------------------------------------------------------------------------
# 8. Report
# -----------------------------------------------------------------------------

message("\n============================================================")
message("CORRELATION EXCLUSION DIAGNOSTICS COMPLETED")
message("============================================================")
message("Tables saved to: ", table_dir)
message("")
message("1. correlation_exclusion_pairs_at_removal.csv")
message("   Exact predictor pair responsible for each exclusion.")
message("")
message("2. correlation_exclusion_vs_all_final_retained_predictors.csv")
message("   Pearson and Spearman correlations with every final survivor.")
message("")
message("3. correlation_exclusion_status_by_response.csv")
message("   Whether each target was eligible and whether it was excluded.")
message("")
message("4. strongest_retained_correlate_for_each_excluded_predictor.csv")
message("   Strongest final retained correlate for each excluded predictor.")
message("============================================================")
