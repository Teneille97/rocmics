# =============================================================================
# MALLE ERW: Bayesian treatment comparisons of integrated dataset
#
# Purpose:
#   Compare Control, Lime, HB 50, EG 50 and BD 50 using Bayesian regression.
#   Generate Hühnerberg dose-response line plots across application rates.
#
# Included response categories:
#   A. Soil chemistry
#   B. CEC and exchangeable cations
#   C. Oxalate/DCB-extractable pools and ratios
#   D. Tessier carbonate fractions (Acetic acid)
#   E. Tessier exchangeable fractions (NH4OAc)
#   F. Tessier oxidizable fractions (H2O2)
#   G. Tessier reducible fractions (Hydroxylamine)
#   H. Aggregate stability indices
#   I. FTIR spectral changes
#   J. Carbon stocks
#   K. Radiocarbon signatures
#
# Comparisons:
#   1. All pairwise treatment comparisons (grouped by parameter type)
#   2. Primary treatment-vs-Control comparisons (Lime, HB 50, EG 50, BD 50 minus Control)
#
# Model:
#   Standardized response ~ 0 + Treatment_PCA
#   Student-t likelihood, df = 4
#   Regularizing priors
#
# Outputs:
#   Posterior treatment means
#   All pairwise treatment contrasts (grouped by parameter type)
#   Primary Control-vs contrasts with CSV significance flag (* if 95% CI excludes 0 & prob > 0.75)
#   Standardized forest plots (no stars)
#   Hühnerberg dose-response line plots across application rates
# =============================================================================


# ---- Code section 24: Bayesian treatment comparisons ------------------------

# =============================================================================
# 24.1 Packages
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(cmdstanr)
library(posterior)
library(readr)
library(tibble)
library(here)


# =============================================================================
# 24.2 Output directories
# =============================================================================
csv_dir <- here::here("csv_files")
outputs_dir <- here::here("outputs")

bayesian_treatment_out_dir <- file.path(
  outputs_dir,
  "Bayesian_treatment_comparisons"
)

bayesian_treatment_models_dir <- file.path(
  bayesian_treatment_out_dir,
  "models"
)

bayesian_treatment_plots_dir <- file.path(
  bayesian_treatment_out_dir,
  "plots"
)

bayesian_treatment_diagnostics_dir <- file.path(
  bayesian_treatment_out_dir,
  "diagnostics"
)


dir.create(
  bayesian_treatment_out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  bayesian_treatment_models_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  bayesian_treatment_plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  bayesian_treatment_diagnostics_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 24.0 Load integrated Malle ERW dataset
# =============================================================================

integrated <- read_csv(
  here::here(
    "outputs",
    "integrated_geochemistry_for_bayesian.csv"
  ),
  na = c("", "NA", "#VALUE!"),
  show_col_types = FALSE
)

# Check dataset
glimpse(integrated)

# Check treatment labels
print(unique(integrated$Treatment))


# =============================================================================
# 24.3 Define treatment selection
# =============================================================================

bayesian_treatment_levels <- c(
  "Control",
  "Lime",
  "HB 50",
  "EG 50",
  "BD 50"
)


bayesian_treatment_data <- integrated %>%
  filter(
    Treatment == "Control" |
      grepl("^Lime", Treatment) |
      Treatment %in% c(
        "Huhnerberg_50",
        "Eifelgold_50",
        "Bolsdorfer_50"
      )
  ) %>%
  mutate(
    Treatment_PCA = case_when(
      Treatment == "Control" ~ "Control",
      grepl("^Lime", Treatment) ~ "Lime",
      Treatment == "Huhnerberg_50" ~ "HB 50",
      Treatment == "Eifelgold_50" ~ "EG 50",
      Treatment == "Bolsdorfer_50" ~ "BD 50",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Treatment_PCA)) %>%
  mutate(
    Treatment_PCA = factor(
      Treatment_PCA,
      levels = bayesian_treatment_levels
    )
  )


# Check treatment replicate numbers
bayesian_treatment_replicates <- bayesian_treatment_data %>%
  count(Treatment_PCA, name = "n")

print(bayesian_treatment_replicates)

write_csv(
  bayesian_treatment_replicates,
  file.path(bayesian_treatment_out_dir, "treatment_replicate_counts.csv")
)


# =============================================================================
# 24.4 Define response variables and categorized parameter lists
# =============================================================================

parameter_categories <- list(
  "Soil_Chemistry" = c(
    "soil_pH",
    "soil_EC"
  ),
  "CEC_Cations" = c(
    "CEC_CEC",
    "CEC_Base_saturation",
    "CEC_H",
    "CEC_Ca",
    "CEC_K",
    "CEC_Mg",
    "CEC_Na",
    "CEC_Al",
    "CEC_Fe",
    "CEC_Mn"
  ),
  "Oxalate_Pools_Ratios" = c(
    "Fe_ox_over_DCB",
    "Al_ox_over_DCB",
    "Mn_ox_over_DCB",
    "Si_ox_over_DCB",
    "Fe_ox_over_Si_ox",
    "Al_ox_over_Si_ox"
  ),
  "Tessier_Carbonate" = c(
    "Tessier_Acetic_Ca_mg_g",
    "Tessier_Acetic_K_mg_g",
    "Tessier_Acetic_Na_mg_g",
    "Tessier_Acetic_Mg_mg_g"
  ),
  "Tessier_Exchangeable" = c(
    "Tessier_NH4OAc_Ca_mg_g",
    "Tessier_NH4OAc_K_mg_g",
    "Tessier_NH4OAc_Na_mg_g",
    "Tessier_NH4OAc_Mg_mg_g"
  ),
  "Tessier_Oxidizable" = c(
    "Tessier_H2O2_Fe_mg_g",
    "Tessier_H2O2_Al_mg_g",
    "Tessier_H2O2_Mn_mg_g",
    "Tessier_H2O2_Si_mg_g"
  ),
  "Tessier_Reducible" = c(
    "Tessier_Hydroxylamine_Fe_mg_g",
    "Tessier_Hydroxylamine_Al_mg_g",
    "Tessier_Hydroxylamine_Mn_mg_g",
    "Tessier_Hydroxylamine_Si_mg_g"
  ),
  "Aggregate_Stability" = c(
    "GMD_h2o",
    "D_h2o",
    "MWD_h2o",
    "GMD_oxalate_over_h2o",
    "D_oxalate_over_h2o",
    "MWD_oxalate_over_h2o"
  ),
  "FTIR_Spectral" = c(
    "FTIR_delta_aromatic",
    "FTIR_delta_Carbonate_bend",
    "FTIR_delta_carbonate_stretch",
    "FTIR_delta_carbonyl",
    "FTIR_delta_Aliphatic_total",
    "FTIR_delta_COO_total",
    "FTIR_delta_polysaccharide",
    "FTIR_delta_poorly_crystalline",
    "FTIR_delta_Silicate",
    "FTIR_delta_total_area"
  ),
  "Carbon_Stocks" = c(
    "C_Bulk",
    "C_MAOM",
    "C_oPOM",
    "C_fPOM"
  ),
  "Radiocarbon_Signatures" = c(
    "d14C_2026_Bulk",
    "d14C_2026_MAOM",
    "d14C_2026_oPOM",
    "d14C_2026_fPOM"
  )
)

bayesian_treatment_parameters <- unname(unlist(parameter_categories))

# Retain only variables present in integrated dataset
bayesian_treatment_parameters_available <- intersect(
  bayesian_treatment_parameters,
  names(bayesian_treatment_data)
)

bayesian_treatment_parameters_missing <- setdiff(
  bayesian_treatment_parameters,
  names(bayesian_treatment_data)
)

if (length(bayesian_treatment_parameters_missing) > 0) {
  warning(
    "Variables not found in integrated: ",
    paste(bayesian_treatment_parameters_missing, collapse = ", ")
  )
}

cat("\nNumber of response variables available:", length(bayesian_treatment_parameters_available), "\n")


# =============================================================================
# 24.5 Define pairwise and primary comparisons (Strictly Versus Control)
# =============================================================================

# All 10 pairwise comparisons
bayesian_treatment_all_pairs <- combn(
  bayesian_treatment_levels,
  2,
  simplify = FALSE
) %>%
  map_dfr(
    function(x) {
      tibble(
        Treatment_1 = x[1],
        Treatment_2 = x[2]
      )
    }
  ) %>%
  mutate(
    Comparison = paste(Treatment_1, "minus", Treatment_2)
  )

# Primary comparisons: Each treatment minus Control (Lime, HB 50, EG 50, BD 50 minus Control)
bayesian_treatment_primary_pairs <- tibble(
  Treatment_1 = c("Lime", "HB 50", "EG 50", "BD 50"),
  Treatment_2 = c("Control", "Control", "Control", "Control")
) %>%
  mutate(
    Comparison = paste(Treatment_1, "minus", Treatment_2)
  )

print(bayesian_treatment_all_pairs)
print(bayesian_treatment_primary_pairs)


# =============================================================================
# 24.6 CmdStan model settings and compilation
# =============================================================================

bayesian_treatment_chains <- 4
bayesian_treatment_iterations <- 4000
bayesian_treatment_warmup <- 2000
bayesian_treatment_cores <- min(4, parallel::detectCores())
bayesian_treatment_seed <- 20261109

cat("CmdStan version:", as.character(cmdstanr::cmdstan_version()), "\n")
cat("CmdStan path:", cmdstanr::cmdstan_path(), "\n")

bayesian_treatment_stan_file <- file.path(
  bayesian_treatment_out_dir,
  "bayesian_treatment_comparison.stan"
)

bayesian_treatment_stan_code <- '
data {
  int<lower=1> N;
  int<lower=2> K;
  array[N] int<lower=1, upper=K> treatment;
  vector[N] y;
}

parameters {
  vector[K] mu;
  real<lower=0> sigma;
  real<lower=2> nu;
}

model {
  mu ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);
  nu ~ gamma(2, 0.1);

  for (n in 1:N) {
    y[n] ~ student_t(nu, mu[treatment[n]], sigma);
  }
}

generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;

  for (n in 1:N) {
    y_rep[n] = student_t_rng(nu, mu[treatment[n]], sigma);
    log_lik[n] = student_t_lpdf(y[n] | nu, mu[treatment[n]], sigma);
  }
}
'

writeLines(bayesian_treatment_stan_code, bayesian_treatment_stan_file)

bayesian_treatment_stan_model <- cmdstanr::cmdstan_model(
  stan_file = bayesian_treatment_stan_file,
  quiet = TRUE
)

cat("\nStan model compiled successfully.\n")


# =============================================================================
# 24.7 Helper function: Fit one Bayesian model & calculate contrasts with Sig Flag
# =============================================================================

fit_bayesian_treatment_model <- function(response_variable, data) {
  
  model_data <- data %>%
    select(Treatment_PCA, all_of(response_variable)) %>%
    rename(response = all_of(response_variable)) %>%
    drop_na(response, Treatment_PCA) %>%
    mutate(
      Treatment_PCA = factor(
        Treatment_PCA,
        levels = bayesian_treatment_levels
      )
    )
  
  treatment_counts <- model_data %>%
    count(Treatment_PCA, name = "n")
  
  if (n_distinct(model_data$Treatment_PCA) != 5) {
    return(list(status = "Skipped", reason = "Not all five treatment groups represented", response = response_variable))
  }
  
  if (any(treatment_counts$n < 2)) {
    return(list(status = "Skipped", reason = "At least one treatment has < 2 observations", response = response_variable))
  }
  
  response_mean <- mean(model_data$response)
  response_sd <- sd(model_data$response)
  
  if (!is.finite(response_sd) || response_sd == 0) {
    return(list(status = "Skipped", reason = "Response has zero or invalid variance", response = response_variable))
  }
  
  model_data <- model_data %>%
    mutate(response_z = (response - response_mean) / response_sd)
  
  stan_data <- list(
    N = nrow(model_data),
    K = length(bayesian_treatment_levels),
    treatment = as.integer(model_data$Treatment_PCA),
    y = as.numeric(model_data$response_z)
  )
  
  model_fit <- bayesian_treatment_stan_model$sample(
    data = stan_data,
    chains = bayesian_treatment_chains,
    parallel_chains = bayesian_treatment_cores,
    iter_warmup = bayesian_treatment_warmup,
    iter_sampling = (bayesian_treatment_iterations - bayesian_treatment_warmup),
    seed = bayesian_treatment_seed,
    adapt_delta = 0.99,
    max_treedepth = 12,
    refresh = 0,
    output_dir = bayesian_treatment_models_dir,
    output_basename = paste0(response_variable, "_chain")
  )
  
  posterior_mu_draws <- model_fit$draws(variables = "mu", format = "matrix")
  
  posterior_means <- do.call(
    cbind,
    lapply(
      seq_along(bayesian_treatment_levels),
      function(k) {
        posterior_mu_draws[, paste0("mu[", k, "]")]
      }
    )
  )
  colnames(posterior_means) <- bayesian_treatment_levels
  
  # Posterior treatment means summary
  posterior_treatment_means <- map_dfr(
    bayesian_treatment_levels,
    function(treatment) {
      draws_z <- posterior_means[, treatment]
      draws_original <- (draws_z * response_sd) + response_mean
      
      tibble(
        Response = response_variable,
        Treatment = treatment,
        Mean_standardized = mean(draws_z),
        Median_standardized = median(draws_z),
        Lower_95_standardized = quantile(draws_z, 0.025),
        Upper_95_standardized = quantile(draws_z, 0.975),
        Mean_original = mean(draws_original),
        Median_original = median(draws_original),
        Lower_95_original = quantile(draws_original, 0.025),
        Upper_95_original = quantile(draws_original, 0.975)
      )
    }
  )
  
  # All pairwise contrasts with significance flag (CI excludes 0 AND Prob > 0.75)
  posterior_pairwise_contrasts <- map_dfr(
    seq_len(nrow(bayesian_treatment_all_pairs)),
    function(i) {
      t1 <- bayesian_treatment_all_pairs$Treatment_1[i]
      t2 <- bayesian_treatment_all_pairs$Treatment_2[i]
      
      draws_diff_z <- posterior_means[, t1] - posterior_means[, t2]
      draws_diff_orig <- draws_diff_z * response_sd
      
      p_higher <- mean(draws_diff_z > 0)
      p_lower  <- mean(draws_diff_z < 0)
      lower_q  <- quantile(draws_diff_z, 0.025)
      upper_q  <- quantile(draws_diff_z, 0.975)
      
      # Flag if Credible Interval excludes zero AND probability > 0.75 for being higher or lower
      is_sig <- (lower_q > 0 | upper_q < 0) & (p_higher > 0.75 | p_lower > 0.75)
      sig_flag <- if_else(is_sig, "*", "")
      
      tibble(
        Response = response_variable,
        Treatment_1 = t1,
        Treatment_2 = t2,
        Comparison = paste(t1, "minus", t2),
        Difference_standardized = mean(draws_diff_z),
        Median_difference_standardized = median(draws_diff_z),
        Lower_95_standardized = lower_q,
        Upper_95_standardized = upper_q,
        Difference_original = mean(draws_diff_orig),
        Median_difference_original = median(draws_diff_orig),
        Lower_95_original = quantile(draws_diff_orig, 0.025),
        Upper_95_original = quantile(draws_diff_orig, 0.975),
        Probability_Treatment1_higher = p_higher,
        Probability_Treatment1_lower = p_lower,
        Sig_Flag = sig_flag
      )
    }
  )
  
  # Primary contrasts (versus Control) with significance flag
  posterior_primary_contrasts <- map_dfr(
    seq_len(nrow(bayesian_treatment_primary_pairs)),
    function(i) {
      t1 <- bayesian_treatment_primary_pairs$Treatment_1[i]
      t2 <- bayesian_treatment_primary_pairs$Treatment_2[i]
      
      draws_diff_z <- posterior_means[, t1] - posterior_means[, t2]
      draws_diff_orig <- draws_diff_z * response_sd
      
      p_higher <- mean(draws_diff_z > 0)
      p_lower  <- mean(draws_diff_z < 0)
      lower_q  <- quantile(draws_diff_z, 0.025)
      upper_q  <- quantile(draws_diff_z, 0.975)
      
      is_sig <- (lower_q > 0 | upper_q < 0) & (p_higher > 0.75 | p_lower > 0.75)
      sig_flag <- if_else(is_sig, "*", "")
      
      tibble(
        Response = response_variable,
        Treatment_1 = t1,
        Treatment_2 = t2,
        Comparison = paste(t1, "minus", t2),
        Difference_standardized = mean(draws_diff_z),
        Median_difference_standardized = median(draws_diff_z),
        Lower_95_standardized = lower_q,
        Upper_95_standardized = upper_q,
        Difference_original = mean(draws_diff_orig),
        Median_difference_original = median(draws_diff_orig),
        Lower_95_original = quantile(draws_diff_orig, 0.025),
        Upper_95_original = quantile(draws_diff_orig, 0.975),
        Probability_Treatment1_higher = p_higher,
        Probability_Treatment1_lower = p_lower,
        Sig_Flag = sig_flag
      )
    }
  )
  
  diagnostic_draws <- model_fit$draws(variables = c("mu", "sigma", "nu"))
  diagnostic_summary <- posterior::summarise_draws(
    diagnostic_draws,
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail
  )
  sampler_diagnostics <- model_fit$diagnostic_summary(quiet = TRUE)
  
  model_diagnostics <- tibble(
    Response = response_variable,
    Maximum_Rhat = max(diagnostic_summary$rhat, na.rm = TRUE),
    Minimum_bulk_ESS = min(diagnostic_summary$ess_bulk, na.rm = TRUE),
    Minimum_tail_ESS = min(diagnostic_summary$ess_tail, na.rm = TRUE),
    Divergent_transitions = sum(sampler_diagnostics$num_divergent),
    Observations = nrow(model_data)
  )
  
  list(
    status = "Success",
    response = response_variable,
    fit = model_fit,
    model_data = model_data,
    treatment_means = posterior_treatment_means,
    pairwise_contrasts = posterior_pairwise_contrasts,
    primary_contrasts = posterior_primary_contrasts,
    diagnostics = model_diagnostics
  )
}


# =============================================================================
# 24.8 Fit models for all parameters
# =============================================================================

test_soil_pH <- fit_bayesian_treatment_model(
  response_variable = "soil_pH",
  data = bayesian_treatment_data
)

if (test_soil_pH$status != "Success") {
  stop("The soil_pH test model did not succeed. Resolve this before fitting all models.")
}

bayesian_treatment_results <- setNames(
  vector("list", length(bayesian_treatment_parameters_available)),
  bayesian_treatment_parameters_available
)

bayesian_treatment_results[["soil_pH"]] <- test_soil_pH

for (response_variable in setdiff(bayesian_treatment_parameters_available, "soil_pH")) {
  cat("\nFitting Bayesian model for:", response_variable, "\n")
  
  result <- tryCatch(
    fit_bayesian_treatment_model(
      response_variable = response_variable,
      data = bayesian_treatment_data
    ),
    error = function(e) {
      warning("Model failed for ", response_variable, ": ", conditionMessage(e))
      list(status = "Failed", response = response_variable, reason = conditionMessage(e))
    }
  )
  
  bayesian_treatment_results[[response_variable]] <- result
  
  if (identical(result$status, "Success")) {
    saveRDS(
      result$fit,
      file.path(
        bayesian_treatment_models_dir,
        paste0(response_variable, "_Bayesian_treatment_model.rds")
      )
    )
  }
}


# =============================================================================
# 24.9 Combine results & export CSVs
# =============================================================================

bayesian_treatment_successful_results <- keep(bayesian_treatment_results, ~ .x$status == "Success")
bayesian_treatment_skipped_results    <- keep(bayesian_treatment_results, ~ .x$status != "Success")

bayesian_treatment_means_all <- map_dfr(bayesian_treatment_successful_results, "treatment_means")
bayesian_treatment_contrasts_all <- map_dfr(bayesian_treatment_successful_results, "pairwise_contrasts")
bayesian_treatment_primary_contrasts_all <- map_dfr(bayesian_treatment_successful_results, "primary_contrasts")
bayesian_treatment_diagnostics_all <- map_dfr(bayesian_treatment_successful_results, "diagnostics")

write_csv(bayesian_treatment_means_all, file.path(bayesian_treatment_out_dir, "Bayesian_treatment_posterior_means.csv"))
write_csv(bayesian_treatment_contrasts_all, file.path(bayesian_treatment_out_dir, "Bayesian_all_pairwise_treatment_contrasts.csv"))
write_csv(bayesian_treatment_primary_contrasts_all, file.path(bayesian_treatment_out_dir, "Bayesian_primary_basalt_vs_Control_Lime.csv"))
write_csv(bayesian_treatment_diagnostics_all, file.path(bayesian_treatment_diagnostics_dir, "Bayesian_model_diagnostics.csv"))


# =============================================================================
# 24.10 Generate Parameter-Grouped Plots (Forest, Violin, & Dose-Response)
# =============================================================================

# Helper function to attach parameter category metadata
add_category_metadata <- function(df) {
  cat_map <- stack(parameter_categories) %>%
    as_tibble() %>%
    rename(Response = values, Category = ind)
  
  df %>% left_join(cat_map, by = "Response")
}

primary_contrasts_categorized <- add_category_metadata(bayesian_treatment_primary_contrasts_all)
pairwise_contrasts_categorized <- add_category_metadata(bayesian_treatment_contrasts_all)
raw_data_categorized <- bayesian_treatment_data %>%
  pivot_longer(
    cols = any_of(bayesian_treatment_parameters_available),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  drop_na(Value) %>%
  add_category_metadata()

# -----------------------------------------------------------------------------
# 1. Primary Forest Plots in Standardized Units (No Stars on Plots)
# -----------------------------------------------------------------------------
for (cat_name in names(parameter_categories)) {
  cat_data <- primary_contrasts_categorized %>% filter(Category == cat_name)
  if (nrow(cat_data) == 0) next
  
  p_forest <- ggplot(cat_data, aes(x = Difference_standardized, y = Response, color = Comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = Lower_95_standardized, xmax = Upper_95_standardized), height = 0.2, position = position_dodge(width = 0.5)) +
    labs(
      title = paste("Primary Treatment Contrasts (vs Control) -", cat_name),
      subtitle = "Differences in Standardized Units (95% Credible Intervals)",
      x = "Standardized Difference minus Control (SD Units)",
      y = "Parameter"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(
    filename = file.path(bayesian_treatment_plots_dir, paste0("Primary_Forest_Plot_Standardized_", cat_name, ".png")),
    plot = p_forest,
    width = 9, height = max(4, length(unique(cat_data$Response)) * 0.6), dpi = 300
  )
}

# -----------------------------------------------------------------------------
# 2. Violin Plots Grouped by Parameter Category
# -----------------------------------------------------------------------------
for (cat_name in names(parameter_categories)) {
  cat_raw <- raw_data_categorized %>% filter(Category == cat_name)
  if (nrow(cat_raw) == 0) next
  
  p_violin <- ggplot(cat_raw, aes(x = Treatment_PCA, y = Value, fill = Treatment_PCA)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.8, outlier.shape = NA) +
    facet_wrap(~ Response, scales = "free_y") +
    labs(
      title = paste("Treatment Distribution (Violin Plots) -", cat_name),
      x = "Treatment",
      y = "Measured Value"
    ) +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(
    filename = file.path(bayesian_treatment_plots_dir, paste0("Violin_Plots_", cat_name, ".png")),
    plot = p_violin,
    width = 10, height = 8, dpi = 300
  )
}

# -----------------------------------------------------------------------------
# 3. Hühnerberg Dose-Response Line Plots Across Application Rates
# -----------------------------------------------------------------------------
# Extract all Hühnerberg treatments and Control across application rates
hb_dose_data <- integrated %>%
  filter(
    Treatment == "Control" | 
      grepl("^Huhnerberg", Treatment) | 
      grepl("^HB_", Treatment) | 
      grepl("^HB[0-9]", Treatment)
  ) %>%
  mutate(
    # Extract numerical dose from treatment string (Control = 0)
    Dose_t_ha = case_when(
      Treatment == "Control" ~ 0,
      TRUE ~ as.numeric(gsub("[^0-9]", "", Treatment))
    )
  ) %>%
  filter(!is.na(Dose_t_ha)) %>%
  pivot_longer(
    cols = any_of(bayesian_treatment_parameters_available),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  drop_na(Value) %>%
  add_category_metadata()

for (cat_name in names(parameter_categories)) {
  cat_hb <- hb_dose_data %>% filter(Category == cat_name)
  if (nrow(cat_hb) == 0) next
  
  p_dose_hb <- ggplot(cat_hb, aes(x = Dose_t_ha, y = Value)) +
    geom_line(color = "darkblue", size = 0.8) +
    geom_point(color = "darkblue", size = 2.5) +
    facet_wrap(~ Response, scales = "free_y") +
    labs(
      title = paste("Hühnerberg Dose-Response Trends -", cat_name),
      subtitle = "Application Rates (0, 10, 20, 50, 100 t/ha)",
      x = "Hühnerberg Application Rate (t/ha)",
      y = "Measured Parameter Value"
    ) +
    theme_bw() +
    scale_x_continuous(breaks = unique(cat_hb$Dose_t_ha))
  
  ggsave(
    filename = file.path(bayesian_treatment_plots_dir, paste0("Huhnerberg_Dose_Response_", cat_name, ".png")),
    plot = p_dose_hb,
    width = 10, height = 8, dpi = 300
  )
}

cat("\nAll modeling, standardized forest plots, and Hühnerberg dose-response plots successfully completed!\n")