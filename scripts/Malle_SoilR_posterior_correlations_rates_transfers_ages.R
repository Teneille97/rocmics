# =============================================================================
# Malle ERW: Posterior-draw correlations between SoilR rates, transfer
# coefficients, pool ages, system transit time, and treatment-level geochemistry
#
# Run from the RStudio project root.
#
# Uses existing MCMC_3p_*.Rdata files; no model refitting is performed.
# The four derived persistence metrics follow the definitions in
# Malle_Bayesian_SoilR_Geochemistry(3).R:
#   MAOM_age, oPOM_age, fPOM_age, system_transit_time
# and use burn-in = 20,000, 1,000 posterior parameter draws per treatment,
# and an age grid from 0 to 500 years in 1-year increments.
#
# The correlation heatmaps combine:
#   - decomposition rates: kf, ki, ks
#   - transfer coefficients: alpha_fi, alpha_is
#   - derived pool mean ages: fPOM_age, oPOM_age, MAOM_age
#   - system mean transit time: system_transit_time
#
# IMPORTANT:
# There are only FIVE independent treatment-level observations. Posterior
# draws propagate uncertainty in SoilR-derived quantities but do not increase
# the number of independent treatments. Results are exploratory, not a
# Bayesian model of correlation and not evidence of causality. No p-values.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(SoilR)
})

set.seed(20261028)

# -----------------------------------------------------------------------------
# 0. Paths and settings
# -----------------------------------------------------------------------------

integrated_file <- here::here("outputs", "integrated_geochemistry_for_bayesian.csv")
mcmc_dir <- here::here("mod_runs")
out_dir <- here::here("outputs", "soilR_posterior_correlations_age_transfer")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

burnin <- 20000L
n_posterior_draws <- 1000L       # matches the supplied age/transit-time script
n_correlation_draws <- 5000L
age_grid <- seq(0, 500, by = 1)
mean_C_inputs <- 100

mcmc_files <- c(
  Control = "MCMC_3p_Control.Rdata",
  Lime = "MCMC_3p_Lime.Rdata",
  Bolsdorfer = "MCMC_3p_Bolsdorfer.Rdata",
  Huhnerberg = "MCMC_3p_Huhnerberg.Rdata",
  Eifelgold = "MCMC_3p_Eifelgold.Rdata"
)
model_treatments <- names(mcmc_files)
soilR_parameters <- c("kf", "ki", "ks", "alpha_fi", "alpha_is")
age_outcomes <- c("fPOM_age", "oPOM_age", "MAOM_age", "system_transit_time")
all_outcomes <- c(soilR_parameters, age_outcomes)

response_columns <- c(
  "C_Bulk", "C_MAOM", "C_oPOM", "C_fPOM",
  "d14C_2026_Bulk", "d14C_2026_MAOM", "d14C_2026_oPOM", "d14C_2026_fPOM"
)
metadata_patterns <- c(
  "^sample$", "^Treatment($|_)", "^Treatment_model$", "^Basalt_type",
  "^App_rate$", "^Dose$", "^Date$", "^Year$", "^Depth$", "^Rep$",
  "^Plot$", "^plot_id$", "^sample_id$", "^date$", "^year$", "^depth$",
  "^d14C", "^Delta14C", "^Fm$"
)

# -----------------------------------------------------------------------------
# 1. Read and summarise treatment-level geochemistry
# -----------------------------------------------------------------------------

if (!file.exists(integrated_file)) {
  stop("Cannot find ", integrated_file,
       ". Expected outputs/integrated_geochemistry_for_bayesian.csv.")
}
integrated <- readr::read_csv(integrated_file, show_col_types = FALSE)
if (!all(c("sample", "Treatment") %in% names(integrated))) {
  stop("The integrated CSV must contain columns named sample and Treatment.")
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

integrated <- integrated %>%
  filter(Treatment %in% c(
    "Control", "Lime_2", "Bolsdorfer_50", "Huhnerberg_50", "Eifelgold_50"
  )) %>%
  mutate(Treatment_model = normalise_treatment(Treatment))

if (anyNA(integrated$Treatment_model)) {
  stop("Unmapped treatment labels: ",
       paste(unique(integrated$Treatment[is.na(integrated$Treatment_model)]),
             collapse = ", "))
}
treatment_check <- integrated %>%
  count(Treatment_model, name = "n_plots") %>%
  arrange(match(Treatment_model, model_treatments))
print(treatment_check, n = Inf)
if (nrow(treatment_check) != 5L || any(treatment_check$n_plots != 3L) ||
    !setequal(treatment_check$Treatment_model, model_treatments)) {
  stop("Expected exactly five modeled treatments with three plot rows each.")
}

metadata_regex <- paste(metadata_patterns, collapse = "|")
num_cols <- names(integrated)[vapply(integrated, is.numeric, logical(1))]
predictor_candidates <- setdiff(num_cols, c("sample", response_columns))
predictor_candidates <- predictor_candidates[
  !stringr::str_detect(predictor_candidates,
                       regex(metadata_regex, ignore_case = TRUE))
]
predictor_candidates <- predictor_candidates[
  !stringr::str_detect(
    predictor_candidates,
    regex("^(C_total|C_stock|SOC_stock|SIC_stock|CaCO3_stock|carbon_stock)$",
          ignore_case = TRUE)
  )
]
if (length(predictor_candidates) == 0L) {
  stop("No numeric geochemical predictor columns remain after exclusions.")
}

geochem_summary <- integrated %>%
  group_by(Treatment_model) %>%
  summarise(
    across(all_of(predictor_candidates), list(
      mean = ~ if (all(!is.finite(.x))) NA_real_ else mean(.x[is.finite(.x)]),
      sd = ~ if (sum(is.finite(.x)) < 2L) NA_real_ else sd(.x[is.finite(.x)]),
      n = ~ sum(is.finite(.x))
    ), .names = "{.col}__{.fn}"),
    .groups = "drop"
  ) %>%
  arrange(match(Treatment_model, model_treatments))
readr::write_csv(geochem_summary,
                 file.path(out_dir, "geochemistry_treatment_summary.csv"))

# -----------------------------------------------------------------------------
# 2. Helpers to reproduce SoilR ages and transit time from each parameter draw
# -----------------------------------------------------------------------------

trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

mean_age <- function(age, density) {
  d <- as.numeric(density)
  d[!is.finite(d) | d < 0] <- 0
  area <- trapz(age, d)
  if (!is.finite(area) || area <= 0) return(NA_real_)
  trapz(age, age * d) / area
}

build_A_u <- function(pars) {
  pars <- as.numeric(pars)
  if (length(pars) < 5L) {
    stop("Expected five parameters: kf, ki, ks, alpha_fi, alpha_is.")
  }
  kf <- pars[1]
  ki <- pars[2]
  ks <- pars[3]
  alpha_fi <- pars[4]
  alpha_is <- pars[5]
  A <- diag(-c(kf, ki, ks))
  A[2, 1] <- kf * alpha_fi
  A[3, 2] <- ki * alpha_is
  list(A = A, u = matrix(c(mean_C_inputs, 0, 0), ncol = 1))
}

one_persistence_draw <- function(pars) {
  au <- build_A_u(pars)
  sa <- SoilR::systemAge(A = au$A, u = au$u, a = age_grid)
  tt <- SoilR::transitTime(A = au$A, u = au$u, a = age_grid)
  if (is.null(sa$poolAgeDensity) || ncol(sa$poolAgeDensity) < 3L ||
      is.null(tt$transitTimeDensity)) {
    stop("Unexpected SoilR systemAge()/transitTime() output structure.")
  }
  tibble(
    fPOM_age = mean_age(age_grid, sa$poolAgeDensity[, 1]),
    oPOM_age = mean_age(age_grid, sa$poolAgeDensity[, 2]),
    MAOM_age = mean_age(age_grid, sa$poolAgeDensity[, 3]),
    system_transit_time = mean_age(age_grid, tt$transitTimeDensity)
  )
}

get_par_draws <- function(path, n_draws = n_posterior_draws) {
  if (!file.exists(path)) stop("MCMC file not found: ", path)
  e <- new.env(parent = emptyenv())
  obj_names <- load(path, envir = e)
  if (!"MCMC" %in% obj_names || is.null(e$MCMC$pars)) {
    stop(path, " must contain an object named MCMC with MCMC$pars.")
  }
  pars <- as.data.frame(e$MCMC$pars)
  if (all(soilR_parameters %in% names(pars))) {
    pars <- pars[, soilR_parameters, drop = FALSE]
  } else {
    if (ncol(pars) < length(soilR_parameters)) {
      stop("Fewer than five parameter columns in ", path)
    }
    warning("Parameter names not all found in ", basename(path),
            "; using first five columns in assumed order.")
    pars <- pars[, seq_along(soilR_parameters), drop = FALSE]
    names(pars) <- soilR_parameters
  }
  if (nrow(pars) <= burnin) {
    stop("burnin = ", burnin, " leaves no draws in ", basename(path),
         ". Check chain length and adjust burnin.")
  }
  pars <- pars[(burnin + 1L):nrow(pars), , drop = FALSE]
  pars <- pars[complete.cases(pars) &
                 apply(pars, 1, function(x) all(is.finite(x))), , drop = FALSE]
  if (!nrow(pars)) stop("No finite posterior draws in ", path)
  idx <- sample.int(nrow(pars), size = min(n_draws, nrow(pars)),
                    replace = n_draws > nrow(pars))
  pars[idx, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 3. Sample rate/transfer parameters and derive age/transit-time distributions
# -----------------------------------------------------------------------------

persistence_draws <- purrr::map_dfr(names(mcmc_files), function(tr) {
  path <- file.path(mcmc_dir, mcmc_files[[tr]])
  message("Sampling saved posterior parameters and calculating ages for ", tr)
  pars <- get_par_draws(path)

  ages <- purrr::map_dfr(seq_len(nrow(pars)), function(i) {
    one_persistence_draw(pars[i, ])
  })
  bind_cols(
    tibble(
      Treatment_model = tr,
      draw = seq_len(nrow(pars))
    ),
    pars %>% as_tibble(),
    ages
  )
}) %>%
  arrange(match(Treatment_model, model_treatments), draw)

readr::write_csv(
  persistence_draws,
  file.path(out_dir, "soilR_parameter_and_age_posterior_draws.csv")
)

# -----------------------------------------------------------------------------
# 4. Build treatment-specific posterior distributions for each outcome
# -----------------------------------------------------------------------------

outcome_draw_lists <- setNames(lapply(all_outcomes, function(outcome) {
  setNames(lapply(model_treatments, function(tr) {
    persistence_draws %>%
      filter(Treatment_model == tr) %>%
      pull(all_of(outcome)) %>%
      .[is.finite(.)]
  }), model_treatments)
}), all_outcomes)

x_means <- geochem_summary %>%
  select(Treatment_model, ends_with("__mean")) %>%
  pivot_longer(-Treatment_model, names_to = "predictor_mean", values_to = "value") %>%
  mutate(predictor = stringr::str_remove(predictor_mean, "__mean$")) %>%
  select(Treatment_model, predictor, value) %>%
  pivot_wider(names_from = Treatment_model, values_from = value)

valid_predictors <- x_means %>%
  pivot_longer(
    cols = all_of(model_treatments),
    names_to = "Treatment_model",
    values_to = "value"
  ) %>%
  group_by(predictor) %>%
  summarise(
    n_treatments_available = sum(is.finite(value)),
    n_unique = n_distinct(value[is.finite(value)]),
    .groups = "drop"
  ) %>%
  filter(n_treatments_available >= 4L, n_unique >= 2L) %>%
  pull(predictor)

if (!length(valid_predictors)) {
  stop("No predictors have sufficient finite variation across treatments.")
}
message("Eligible geochemical predictors: ", length(valid_predictors))

x_by_predictor <- setNames(lapply(valid_predictors, function(pred) {
  row <- x_means %>% filter(predictor == pred) %>% slice(1)
  vals <- as.numeric(unlist(row[1, model_treatments], use.names = FALSE))
  names(vals) <- model_treatments
  vals
}), valid_predictors)

# -----------------------------------------------------------------------------
# 5. Draw-propagated Spearman, Kendall and Pearson correlations
# -----------------------------------------------------------------------------

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 4L || dplyr::n_distinct(x) < 2L ||
      dplyr::n_distinct(y) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(x, y, method = method, use = "complete.obs"))
}

methods <- c("spearman", "kendall", "pearson")
correlation_draws <- list()
idx_out <- 1L

for (outcome in all_outcomes) {
  message("Correlating outcome: ", outcome)
  for (predictor in valid_predictors) {
    x <- x_by_predictor[[predictor]][model_treatments]
    if (sum(is.finite(x)) < 4L ||
        dplyr::n_distinct(x[is.finite(x)]) < 2L) next

    y_draw_matrix <- vapply(model_treatments, function(tr) {
      vals <- outcome_draw_lists[[outcome]][[tr]]
      if (!length(vals)) return(rep(NA_real_, n_correlation_draws))
      sample(vals, size = n_correlation_draws, replace = TRUE)
    }, numeric(n_correlation_draws))
    colnames(y_draw_matrix) <- model_treatments

    for (method in methods) {
      rho <- vapply(seq_len(n_correlation_draws), function(i) {
        safe_cor(x, y_draw_matrix[i, ], method = method)
      }, numeric(1))
      rho <- rho[is.finite(rho)]
      if (!length(rho)) next
      correlation_draws[[idx_out]] <- tibble(
        outcome = outcome,
        outcome_group = case_when(
          outcome %in% c("kf", "ki", "ks") ~ "Decomposition rate",
          outcome %in% c("alpha_fi", "alpha_is") ~ "Transfer coefficient",
          outcome %in% c("fPOM_age", "oPOM_age", "MAOM_age") ~ "Pool mean age",
          outcome == "system_transit_time" ~ "System transit time",
          TRUE ~ "Other"
        ),
        predictor = predictor,
        method = method,
        draw = seq_along(rho),
        correlation = rho
      )
      idx_out <- idx_out + 1L
    }
  }
}
correlation_draws <- bind_rows(correlation_draws)

if (!nrow(correlation_draws)) {
  stop("No valid correlation draws were produced; inspect predictor variation.")
}

correlation_summary <- correlation_draws %>%
  group_by(outcome, outcome_group, predictor, method) %>%
  summarise(
    n_correlation_draws = n(),
    correlation_mean = mean(correlation),
    correlation_median = median(correlation),
    correlation_sd = sd(correlation),
    correlation_q025 = unname(quantile(correlation, 0.025)),
    correlation_q975 = unname(quantile(correlation, 0.975)),
    Pr_correlation_positive = mean(correlation > 0),
    Pr_correlation_negative = mean(correlation < 0),
    Pr_abs_correlation_gt_0_5 = mean(abs(correlation) > 0.5),
    .groups = "drop"
  ) %>%
  arrange(method, match(outcome, all_outcomes), desc(abs(correlation_median)))

readr::write_csv(
  correlation_summary,
  file.path(out_dir, "soilR_rates_transfer_ages_geochemistry_correlation_summary.csv")
)
readr::write_csv(
  correlation_draws,
  file.path(out_dir, "soilR_rates_transfer_ages_geochemistry_correlation_draws.csv")
)

# -----------------------------------------------------------------------------
# 6. Heatmaps: combined outcomes and age/transit-time-only view
# -----------------------------------------------------------------------------

outcome_labels <- c(
  kf = "kf: fPOM decomposition rate",
  ki = "ki: oPOM decomposition rate",
  ks = "ks: MAOM decomposition rate",
  alpha_fi = "alpha_fi: fPOM -> oPOM transfer",
  alpha_is = "alpha_is: oPOM -> MAOM transfer",
  fPOM_age = "fPOM mean age",
  oPOM_age = "oPOM mean age",
  MAOM_age = "MAOM mean age",
  system_transit_time = "System mean transit time"
)

plot_heatmap <- function(method_name, outcome_subset, filename_stub, title_text) {
  plot_dat <- correlation_summary %>%
    filter(method == method_name, outcome %in% outcome_subset) %>%
    mutate(
      outcome = factor(outcome, levels = rev(outcome_subset),
                      labels = rev(unname(outcome_labels[outcome_subset]))),
      predictor = factor(predictor, levels = valid_predictors)
    )

  p <- ggplot(plot_dat, aes(x = predictor, y = outcome,
                            fill = correlation_median)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1), oob = scales::squish,
      name = ifelse(method_name == "spearman", "Median rho",
                    ifelse(method_name == "kendall", "Median tau", "Median r"))
    ) +
    labs(
      title = title_text,
      subtitle = paste0(
        "Five treatment-level observations; median ", method_name,
        " correlation across sampled SoilR posterior draws"
      ),
      x = "Geochemical predictor (treatment mean)",
      y = NULL,
      caption = paste(
        "Exploratory only: posterior draws propagate SoilR uncertainty,",
        "but do not increase the number of independent treatments."
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      plot.title = element_text(face = "bold", colour = "black"),
      plot.subtitle = element_text(colour = "black"),
      plot.caption = element_text(colour = "black"),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black"),
      axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1)
    )

  width <- max(12, min(24, 0.28 * length(valid_predictors) + 5))
  height <- max(5.5, 0.48 * length(outcome_subset) + 2.5)
  ggsave(file.path(out_dir, paste0(filename_stub, ".png")), p,
         width = width, height = height, dpi = 300, limitsize = FALSE,
         bg = "white")
  ggsave(file.path(out_dir, paste0(filename_stub, ".pdf")), p,
         width = width, height = height, limitsize = FALSE,
         bg = "white")
  p
}

for (method_name in methods) {
  method_stub <- method_name
  method_label <- stringr::str_to_title(method_name)
  p_all <- plot_heatmap(
    method_name, all_outcomes,
    paste0("heatmap_", method_stub, "_all_rates_transfers_ages"),
    paste0("SoilR rates, transfer coefficients and ages vs geochemistry: ",
           method_label)
  )
  p_ages <- plot_heatmap(
    method_name, age_outcomes,
    paste0("heatmap_", method_stub, "_ages_transit_only"),
    paste0("SoilR pool ages and system transit time vs geochemistry: ",
           method_label)
  )
  print(p_all)
  print(p_ages)
}

# Descriptive shortlist only; this is NOT a significance filter.
correlation_summary %>%
  filter(method == "spearman", abs(correlation_median) >= 0.5) %>%
  arrange(desc(abs(correlation_median))) %>%
  readr::write_csv(
    file.path(out_dir, "spearman_median_absrho_at_least_0.5.csv")
  )

saveRDS(
  list(
    persistence_draws = persistence_draws,
    geochemistry_summary = geochem_summary,
    correlation_draws = correlation_draws,
    correlation_summary = correlation_summary,
    outcomes = all_outcomes,
    predictors = valid_predictors
  ),
  file.path(out_dir, "soilR_rates_transfer_ages_correlation_analysis.rds")
)

message(
  "\nFinished posterior-draw correlation analysis.",
  "\nOutputs saved to: ", out_dir,
  "\nIncludes combined heatmaps (rates + transfer coefficients + ages/transit time),",
  "\nage/transit-only heatmaps, summary/draw CSVs, sampled persistence draws,",
  "\nand an RDS backup.",
  "\nInterpret cautiously: only five independent treatment-level observations."
)
