# =============================================================================
# Malle ERW: Bayesian SoilR persistence vs treatment-level geochemistry
#
# Run from the RStudio project root.
#
# IMPORTANT:
# - Only five independent SoilR treatment outcomes exist.
# - Plot-level geochemistry is summarised within treatment.
# - Never replicate one SoilR outcome across three plots and count those plots
#   as independent response observations.
# - All regressions are exploratory; four models each use only five rows.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(SoilR)
  library(brms)
  library(posterior)
  library(dplyr)
})

set.seed(20261026)

integrated_file <- here::here(
  "outputs", "integrated_geochemistry_for_bayesian.csv"
)
mcmc_dir <- here::here("mod_runs")
out_dir <- here::here("outputs", "bayesian_persistence_withmodfit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

burnin <- 20000L
n_posterior_draws <- 1000L
age_grid <- seq(0, 500, by = 1)
mean_C_inputs <- 100

# Edit filenames if your saved MCMC files differ.
mcmc_files <- c(
  Control = "MCMC_3p_Control.Rdata",
  Lime = "MCMC_3p_Lime.Rdata",
  Bolsdorfer = "MCMC_3p_Bolsdorfer.Rdata",
  Huhnerberg = "MCMC_3p_Huhnerberg.Rdata",
  Eifelgold = "MCMC_3p_Eifelgold.Rdata"
)

model_treatments <- c(
  "Control", "Lime", "Bolsdorfer", "Huhnerberg", "Eifelgold"
)

# ---------------------------------------------------------------------------
# 1. Treatment labels and integrated geochemistry
# ---------------------------------------------------------------------------

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

if (!file.exists(integrated_file)) {
  stop(
    "Cannot find ", integrated_file,
    ". Export the integrated dataset to ",
    "outputs/integrated_geochemistry_for_bayesian.csv."
  )
}

integrated <- readr::read_csv(
  integrated_file,
  show_col_types = FALSE
)

if (!all(c("sample", "Treatment") %in% names(integrated))) {
  stop("The integrated CSV must contain columns named sample and Treatment.")
}

# Map labels once, then retain only the five treatments represented by the
# five SoilR MCMC files. This excludes other Huhnerberg application rates.
integrated <- integrated %>%
  mutate(Treatment_model = normalise_treatment(Treatment)) %>%
  filter(
    Treatment %in% c(
      "Control",
      "Lime_2",
      "Bolsdorfer_50",
      "Huhnerberg_50",
      "Eifelgold_50"
    )
  )

if (anyNA(integrated$Treatment_model)) {
  stop(
    "Unmapped treatment labels: ",
    paste(
      unique(integrated$Treatment[is.na(integrated$Treatment_model)]),
      collapse = ", "
    )
  )
}

treatment_check <- integrated %>%
  count(Treatment, Treatment_model, name = "n_plots")
print(treatment_check, n = Inf)

if (nrow(treatment_check) != 5L || any(treatment_check$n_plots != 3L)) {
  stop(
    "Expected exactly five modeled treatments with three plots each. ",
    "Check the treatment labels and integrated CSV."
  )
}

# ---------------------------------------------------------------------------
# 2. Predictor specifications
# ---------------------------------------------------------------------------

model_specs <- list(
  MAOM_age = c(
    "soil_pH",
    "Fe_ox_over_DCB",
    "CEC_CEC"
  ),
  oPOM_age = c(
    "soil_EC",
    "CEC_Mg"
  ),
  fPOM_age = c(
    "soil_pH",
    "CEC_K"
  ),
  system_transit_time = c(
    "Fe_ox_over_DCB",
    "CEC_CEC"
  )
)
all_predictors <- unique(unlist(model_specs))

missing_predictors <- setdiff(all_predictors, names(integrated))
if (length(missing_predictors) > 0L) {
  stop(
    "Missing predictor columns in integrated CSV: ",
    paste(missing_predictors, collapse = ", "),
    ". Check exact names and edit model_specs."
  )
}

for (v in all_predictors) {
  integrated[[v]] <- suppressWarnings(as.numeric(integrated[[v]]))
}

# Treatment means, between-plot SDs, sample sizes and standard errors.
geochem_summary <- integrated %>%
  group_by(Treatment_model) %>%
  summarise(
    across(
      all_of(all_predictors),
      list(
        mean = ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE),
        sd = ~ if (sum(is.finite(.x)) < 2L) NA_real_ else sd(.x, na.rm = TRUE),
        n = ~ sum(is.finite(.x))
      ),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  )

if (nrow(geochem_summary) != 5L ||
    !setequal(geochem_summary$Treatment_model, model_treatments)) {
  stop("The geochemistry summary does not contain the expected five treatments.")
}

readr::write_csv(
  geochem_summary,
  file.path(out_dir, "geochemistry_treatment_summary.csv")
)

# ---------------------------------------------------------------------------
# 3. Reconstruct SoilR persistence metrics from posterior parameter draws
# ---------------------------------------------------------------------------

# Parameter order assumed to match the supplied three-pool SoilR script:
# kf, ki, ks, alpha_fi, alpha_is.
build_A_u <- function(pars) {
  pars <- as.numeric(pars)
  if (length(pars) < 5L) {
    stop("Expected five SoilR parameters: kf, ki, ks, alpha_fi, alpha_is.")
  }
  
  kf <- pars[1]
  ki <- pars[2]
  ks <- pars[3]
  alpha_fi <- pars[4]
  alpha_is <- pars[5]
  
  A <- diag(-c(kf, ki, ks))
  A[2, 1] <- kf * alpha_fi
  A[3, 2] <- ki * alpha_is
  
  list(
    A = A,
    u = matrix(c(mean_C_inputs, 0, 0), ncol = 1)
  )
}

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

one_persistence_draw <- function(pars) {
  au <- build_A_u(pars)
  
  sa <- SoilR::systemAge(
    A = au$A,
    u = au$u,
    a = age_grid
  )
  tt <- SoilR::transitTime(
    A = au$A,
    u = au$u,
    a = age_grid
  )
  
  if (is.null(sa$poolAgeDensity) ||
      ncol(sa$poolAgeDensity) < 3L ||
      is.null(tt$transitTimeDensity)) {
    stop(
      "Unexpected output from SoilR::systemAge() or transitTime(). ",
      "Check the installed SoilR version and the returned object structure."
    )
  }
  
  tibble(
    MAOM_age = mean_age(age_grid, sa$poolAgeDensity[, 3]),
    oPOM_age = mean_age(age_grid, sa$poolAgeDensity[, 2]),
    fPOM_age = mean_age(age_grid, sa$poolAgeDensity[, 1]),
    system_transit_time = mean_age(age_grid, tt$transitTimeDensity)
  )
}

get_par_draws <- function(path) {
  if (!file.exists(path)) {
    stop("MCMC file not found: ", path)
  }
  
  e <- new.env(parent = emptyenv())
  obj_names <- load(path, envir = e)
  
  if (!"MCMC" %in% obj_names || is.null(e$MCMC$pars)) {
    stop(path, " must contain an object named MCMC with MCMC$pars.")
  }
  
  pars <- as.data.frame(e$MCMC$pars)
  wanted <- c("kf", "ki", "ks", "alpha_fi", "alpha_is")
  
  if (all(wanted %in% names(pars))) {
    pars <- pars[, wanted, drop = FALSE]
  } else {
    if (ncol(pars) < 5L) {
      stop("Fewer than five parameter columns in ", path)
    }
    warning(
      "Using the first five MCMC columns in assumed parameter order for ",
      basename(path), ": ",
      paste(names(pars)[1:5], collapse = ", ")
    )
    pars <- pars[, 1:5, drop = FALSE]
  }
  
  if (nrow(pars) <= burnin) {
    stop(
      "burnin = ", burnin, " leaves no draws in ", basename(path),
      ". Check chain length and adjust burnin."
    )
  }
  
  pars <- pars[(burnin + 1L):nrow(pars), , drop = FALSE]
  pars <- pars[
    complete.cases(pars) &
      apply(pars, 1, function(x) all(is.finite(x))),
    ,
    drop = FALSE
  ]
  
  if (nrow(pars) == 0L) {
    stop("No finite posterior draws in ", path)
  }
  
  pars[
    sample(seq_len(nrow(pars)), min(n_posterior_draws, nrow(pars))),
    ,
    drop = FALSE
  ]
}

# Run full SoilR persistence calculation across posterior draws
persistence_draws_list <- list()

for (tr in names(mcmc_files)) {
  mcmc_path <- file.path(mcmc_dir, mcmc_files[[tr]])
  message("Processing MCMC file for treatment: ", tr, " (", mcmc_path, ")")
  
  pars <- get_par_draws(mcmc_path)
  
  draw_res <- purrr::map_dfr(seq_len(nrow(pars)), function(i) {
    one_persistence_draw(pars[i, ])
  }) %>%
    mutate(Treatment_model = tr, draw = row_number())
  
  persistence_draws_list[[tr]] <- draw_res
}

persistence_draws <- bind_rows(persistence_draws_list)

# ---------------------------------------------------------------------------
# 4. Assemble analysis dataset
# ---------------------------------------------------------------------------

soilr_summary <- persistence_draws %>%
  pivot_longer(
    cols = c(MAOM_age, oPOM_age, fPOM_age, system_transit_time),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  group_by(Treatment_model, outcome) %>%
  summarise(
    y_mean = mean(value, na.rm = TRUE),
    y_se   = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  soilr_summary,
  file.path(out_dir, "soilR_persistence_summary.csv")
)

# Reshape geochem summary for merging with outcomes
geochem_long <- geochem_summary %>%
  pivot_longer(
    cols = -Treatment_model,
    names_to = c("predictor", "stat"),
    names_pattern = "^(.*)__(mean|sd|n)$"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(se = sd / sqrt(n))

geochem_wide_means <- geochem_long %>%
  select(Treatment_model, predictor, mean) %>%
  pivot_wider(names_from = predictor, values_from = mean)

geochem_wide_ses <- geochem_long %>%
  select(Treatment_model, predictor, se) %>%
  pivot_wider(names_from = predictor, values_from = se, names_prefix = "") %>%
  rename_with(~ paste0(.x, "_se"), -Treatment_model)

geochem_prepped <- left_join(geochem_wide_means, geochem_wide_ses, by = "Treatment_model")

model_data_long <- soilr_summary %>%
  left_join(geochem_prepped, by = "Treatment_model")

readr::write_csv(
  model_data_long,
  file.path(out_dir, "treatment_level_model_data_long.csv")
)

# Optional session backup
dir.create("outputs/bayesian_persistence/session_backup",
           recursive = TRUE, showWarnings = FALSE)

saveRDS(
  model_data_long,
  "outputs/bayesian_persistence/session_backup/model_data_long.rds"
)

saveRDS(
  geochem_summary,
  "outputs/bayesian_persistence/session_backup/geochem_summary.rds"
)

# Commented out reading saved data for clean first-time run
# model_data_long <- readRDS(
#   here::here(
#     "outputs", "bayesian_persistence",
#     "session_backup", "model_data_long.rds"
#   )
# )

# ============================================================================
# 5. FME regressions: geochemistry versus SoilR persistence
#
# Uses the already-saved SoilR MCMC files and persistence_draws.
# Does NOT refit the SoilR decomposition models.
#
# IMPORTANT:
# - Five independent treatment-level observations per regression.
# - Predictor means are used; predictor measurement error is not propagated
#   in this FME version.
# - y_se is the SD of the SoilR-derived persistence draws for each treatment.
# - FME modMCMC is likelihood-based; it is not identical to a Stan model
#   with explicit Bayesian priors.
# ============================================================================

library(FME)
library(tidyverse)
library(here)

set.seed(20261027)

# Confirm inputs from sections 1-4 exist.
stopifnot(
  exists("model_data_long"),
  exists("model_specs"),
  exists("out_dir")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# MCMC settings.
niter_regression <- 50000L
burnin_regression <- 10000L
covariance_scale <- 0.01

# ---------------------------------------------------------------------------
# Helper: weighted regression cost
# ---------------------------------------------------------------------------

make_regression_cost <- function(dat, predictors) {
  
  # Standardise predictors and response.
  # Predictor SEs are intentionally not used in this FME implementation.
  x_means <- vapply(
    predictors,
    function(p) mean(dat[[p]]),
    numeric(1)
  )
  
  x_sds <- vapply(
    predictors,
    function(p) sd(dat[[p]]),
    numeric(1)
  )
  
  if (any(!is.finite(x_sds)) || any(x_sds <= 0)) {
    stop("One or more predictors have no between-treatment variation.")
  }
  
  y_center <- mean(dat$y_mean)
  y_scale <- sd(dat$y_mean)
  
  if (!is.finite(y_scale) || y_scale <= 0) {
    stop("The response has no between-treatment variation.")
  }
  
  x_std <- as.data.frame(
    sweep(
      sweep(
        as.matrix(dat[predictors]),
        2, x_means, "-"
      ),
      2, x_sds, "/"
    )
  )
  
  names(x_std) <- predictors
  
  y_std <- (dat$y_mean - y_center) / y_scale
  y_se_std <- dat$y_se / y_scale
  
  if (any(!is.finite(y_std)) ||
      any(!is.finite(y_se_std)) ||
      any(y_se_std <= 0)) {
    stop("Non-finite response values or invalid response uncertainty.")
  }
  
  regression_data <- bind_cols(
    tibble(
      idx = seq_len(nrow(dat)),
      y = y_std,
      y_se = y_se_std
    ),
    x_std
  )
  
  # Parameter names are alpha and one beta per predictor.
  parameter_names <- c("alpha", predictors)
  
  # Predictions are on the standardised response scale.
  predict_standardised <- function(pars) {
    alpha <- pars[["alpha"]]
    betas <- pars[predictors]
    
    as.numeric(
      alpha +
        as.matrix(regression_data[predictors]) %*% betas
    )
  }
  
  # FME modCost weights residuals using the supplied y_se column.
  cost_function <- function(pars) {
    
    pred <- data.frame(
      idx = regression_data$idx,
      y = predict_standardised(pars)
    )
    
    obs <- data.frame(
      idx = regression_data$idx,
      y = regression_data$y,
      y_se = regression_data$y_se
    )
    
    FME::modCost(
      model = pred,
      obs = obs,
      x = "idx",
      err = "y_se"
    )
  }
  
  list(
    data = regression_data,
    predictors = predictors,
    parameter_names = parameter_names,
    cost_function = cost_function,
    predict_standardised = predict_standardised,
    x_means = x_means,
    x_sds = x_sds,
    y_center = y_center,
    y_scale = y_scale
  )
}

# ---------------------------------------------------------------------------
# Helper: summarise regression MCMC parameter draws
# ---------------------------------------------------------------------------

summarise_draws <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    lower_95 = unname(quantile(x, 0.025, na.rm = TRUE)),
    upper_95 = unname(quantile(x, 0.975, na.rm = TRUE)),
    Pr_positive = mean(x > 0, na.rm = TRUE),
    Pr_negative = mean(x < 0, na.rm = TRUE)
  )
}

# ---------------------------------------------------------------------------
# Fit one outcome using modFit + modMCMC
# ---------------------------------------------------------------------------

fit_fme_regression <- function(outcome, predictors) {
  
  dat <- model_data_long %>%
    filter(.data$outcome == .env$outcome) %>%
    arrange(match(Treatment_model, model_treatments))
  
  if (nrow(dat) != 5L ||
      dplyr::n_distinct(dat$Treatment_model) != 5L) {
    stop(outcome, ": expected one row for each of five treatments.")
  }
  
  required <- c("y_mean", "y_se", predictors)
  
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop(
      outcome, ": missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  
  if (!all(vapply(dat[required], is.numeric, logical(1))) ||
      !all(vapply(dat[required], function(x) all(is.finite(x)),
                  logical(1)))) {
    stop(outcome, ": non-numeric or non-finite model data.")
  }
  
  if (any(dat$y_se <= 0)) {
    stop(outcome, ": y_se must be positive.")
  }
  
  message(
    "\nFitting FME regression: ", outcome,
    "\nPredictors: ", paste(predictors, collapse = ", "),
    "\nIndependent treatment observations: ", nrow(dat)
  )
  
  setup <- make_regression_cost(dat, predictors)
  
  # Starting values on standardised scales.
  p0 <- setNames(
    rep(0, length(setup$parameter_names)),
    setup$parameter_names
  )
  
  # 1. Find the best-fitting coefficients.
  mFit <- FME::modFit(
    f = setup$cost_function,
    p = p0,
    method = "Marq"
  )
  
  bestpars <- mFit$par
  
  if (any(!is.finite(bestpars))) {
    stop(outcome, ": modFit returned non-finite coefficients.")
  }
  
  # Save the initial optimisation result.
  saveRDS(
    mFit,
    file.path(out_dir, paste0("modFit_", outcome, ".rds"))
  )
  
  # 2. Obtain the fitted coefficient covariance matrix.
  cov0 <- tryCatch(
    as.matrix(summary(mFit)$cov.scaled),
    error = function(e) NULL
  )
  
  p <- length(bestpars)
  
  if (is.null(cov0) ||
      !all(dim(cov0) == c(p, p)) ||
      any(!is.finite(cov0))) {
    
    warning(
      outcome,
      ": modFit covariance was unavailable or invalid. ",
      "Using a diagonal jump proposal as a fallback."
    )
    
    cov0 <- diag(rep(0.01, p))
    dimnames(cov0) <- list(names(bestpars), names(bestpars))
  }
  
  # Symmetrise and add a small diagonal term for numerical stability.
  cov0 <- (cov0 + t(cov0)) / 2
  
  diag(cov0) <- diag(cov0) +
    max(1e-8, max(abs(diag(cov0)), na.rm = TRUE) * 1e-6)
  
  # The jump argument is a covariance matrix. Scaling it by 0.01
  # scales proposal standard deviations by sqrt(0.01) = 0.1.
  jump_cov <- cov0 * covariance_scale
  
  # 3. Sample coefficient uncertainty using FME.
  MCMC <- FME::modMCMC(
    f = setup$cost_function,
    p = bestpars,
    niter = niter_regression,
    jump = jump_cov,
    var0 = mFit$var_ms_unweighted,
    wvar0 = 1,
    updatecov = 1000,
    burninlength = 1000
  )
  
  saveRDS(
    MCMC,
    file.path(out_dir, paste0("modMCMC_", outcome, ".rds"))
  )
  
  # Retain the post-burn-in draws for summaries.
  draws_all <- as.data.frame(MCMC$pars)
  
  if (!all(setup$parameter_names %in% names(draws_all))) {
    stop(
      outcome, ": MCMC parameter columns differ from expected names. ",
      "Found: ", paste(names(draws_all), collapse = ", ")
    )
  }
  
  draws_all <- draws_all[, setup$parameter_names, drop = FALSE]
  
  post_burn <- min(
    burnin_regression,
    max(0L, nrow(draws_all) - 1L)
  )
  
  draws <- draws_all[
    seq.int(post_burn + 1L, nrow(draws_all)),
    ,
    drop = FALSE
  ]
  
  if (nrow(draws) < 100L) {
    stop(outcome, ": too few post-burn-in MCMC draws.")
  }
  
  readr::write_csv(
    draws,
    file.path(out_dir, paste0("regression_draws_", outcome, ".csv"))
  )
  
  # -------------------------------------------------------------------------
  # Convert coefficients back to original units
  # -------------------------------------------------------------------------
  
  coefficient_rows <- list()
  
  alpha_original <- setup$y_center +
    setup$y_scale * draws$alpha
  
  for (j in seq_along(predictors)) {
    pred <- predictors[j]
    
    beta_std <- draws[[pred]]
    beta_original <- beta_std *
      setup$y_scale / setup$x_sds[pred]
    
    alpha_original <- alpha_original -
      beta_original * setup$x_means[pred]
    
    coefficient_rows[[pred]] <- summarise_draws(beta_std) %>%
      mutate(
        outcome = outcome,
        parameter = pred,
        scale = "standardised",
        .before = 1
      )
    
    coefficient_rows[[paste0(pred, "_original")]] <-
      summarise_draws(beta_original) %>%
      mutate(
        outcome = outcome,
        parameter = pred,
        scale = "original response units per predictor unit",
        .before = 1
      )
  }
  
  coefficient_rows[["intercept_standardised"]] <-
    summarise_draws(draws$alpha) %>%
    mutate(
      outcome = outcome,
      parameter = "Intercept",
      scale = "standardised",
      .before = 1
    )
  
  coefficient_rows[["intercept_original"]] <-
    summarise_draws(alpha_original) %>%
    mutate(
      outcome = outcome,
      parameter = "Intercept",
      scale = "original response units",
      .before = 1
    )
  
  coefficients <- bind_rows(coefficient_rows)
  
  readr::write_csv(
    coefficients,
    file.path(out_dir, paste0("coefficients_", outcome, ".csv"))
  )
  
  # -------------------------------------------------------------------------
  # Predictions at the five observed treatment-level predictor means
  # -------------------------------------------------------------------------
  
  prediction_draws <- t(vapply(
    seq_len(nrow(draws)),
    function(i) {
      pars_i <- as.numeric(draws[i, setup$parameter_names])
      names(pars_i) <- setup$parameter_names
      
      setup$predict_standardised(pars_i) * setup$y_scale +
        setup$y_center
    },
    numeric(nrow(dat))
  ))
  
  prediction_summary <- bind_rows(
    lapply(seq_len(nrow(dat)), function(i) {
      tibble(
        Treatment_model = dat$Treatment_model[i],
        observed = dat$y_mean[i],
        observed_se = dat$y_se[i],
        predicted_median = median(prediction_draws[, i]),
        predicted_lower = unname(quantile(prediction_draws[, i], 0.025)),
        predicted_upper = unname(quantile(prediction_draws[, i], 0.975))
      )
    })
  )
  
  readr::write_csv(
    prediction_summary,
    file.path(out_dir, paste0("predictions_", outcome, ".csv"))
  )
  
  # -------------------------------------------------------------------------
  # Observed-versus-predicted plot
  # -------------------------------------------------------------------------
  
  p_obs_pred <- ggplot(
    prediction_summary,
    aes(x = reorder(Treatment_model, observed))
  ) +
    geom_linerange(
      aes(ymin = predicted_lower, ymax = predicted_upper),
      colour = "#2878B5",
      linewidth = 1
    ) +
    geom_point(
      aes(y = predicted_median),
      colour = "#2878B5",
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = observed - observed_se,
        ymax = observed + observed_se
      ),
      width = 0.15
    ) +
    geom_point(
      aes(y = observed),
      size = 3
    ) +
    coord_flip() +
    labs(
      title = paste("FME regression:", outcome),
      subtitle = paste(
        "Black = observed SoilR estimate ± 1 posterior SD;",
        "blue = median fitted value and 95% coefficient interval"
      ),
      x = "Treatment",
      y = outcome
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    file.path(out_dir, paste0("observed_vs_predicted_", outcome, ".png")),
    p_obs_pred,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  # -------------------------------------------------------------------------
  # MCMC trace and density plots
  # -------------------------------------------------------------------------
  
  draws_long <- draws %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(
      cols = all_of(setup$parameter_names),
      names_to = "parameter",
      values_to = "value"
    )
  
  p_trace <- ggplot(
    draws_long,
    aes(x = iteration, y = value)
  ) +
    geom_line(alpha = 0.5, linewidth = 0.25) +
    facet_wrap(~ parameter, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste("FME MCMC traces:", outcome),
      x = "Post-burn-in iteration",
      y = "Parameter value"
    )
  
  ggsave(
    file.path(out_dir, paste0("trace_", outcome, ".png")),
    p_trace,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  p_density <- ggplot(
    draws_long,
    aes(x = value)
  ) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    facet_wrap(~ parameter, scales = "free") +
    theme_bw() +
    labs(
      title = paste("FME MCMC coefficient densities:", outcome),
      x = "Parameter value",
      y = "Density"
    )
  
  ggsave(
    file.path(out_dir, paste0("density_", outcome, ".png")),
    p_density,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  # Basic sampler summary. This is not equivalent to Stan's R-hat/ESS.
  sampler_summary <- tibble(
    outcome = outcome,
    n_treatments = nrow(dat),
    n_parameters = length(bestpars),
    n_draws = nrow(draws),
    n_accepted = MCMC$naccepted,
    acceptance_percent = 100 * MCMC$naccepted / niter_regression,
    best_fit_SS = mFit$ssr,
    best_MCMC_SS = min(MCMC$SS, na.rm = TRUE)
  )
  
  list(
    mFit = mFit,
    MCMC = MCMC,
    draws = draws,
    data = dat,
    setup = setup,
    coefficients = coefficients,
    predictions = prediction_summary,
    sampler_summary = sampler_summary
  )
}

# ---------------------------------------------------------------------------
# Fit all four Hypothesis B regressions
# ---------------------------------------------------------------------------

fits_fme <- list()

for (nm in names(model_specs)) {
  fits_fme[[nm]] <- fit_fme_regression(
    outcome = nm,
    predictors = model_specs[[nm]]
  )
}

saveRDS(
  fits_fme,
  file.path(out_dir, "all_FME_persistence_regressions.rds")
)

# Combine coefficient summaries, predictions and sampler summaries.
all_coefficients <- bind_rows(
  lapply(fits_fme, function(x) x$coefficients)
)

all_predictions <- bind_rows(
  lapply(names(fits_fme), function(nm) {
    fits_fme[[nm]]$predictions %>%
      mutate(outcome = nm, .before = 1)
  })
)

all_sampler_summaries <- bind_rows(
  lapply(fits_fme, function(x) x$sampler_summary)
)

readr::write_csv(
  all_coefficients,
  file.path(out_dir, "all_FME_coefficients.csv")
)

readr::write_csv(
  all_predictions,
  file.path(out_dir, "all_FME_predictions.csv")
)

readr::write_csv(
  all_sampler_summaries,
  file.path(out_dir, "all_FME_sampler_summaries.csv")
)

print(all_sampler_summaries, n = Inf)

message(
  "\nFinished FME persistence regressions.",
  "\nOutputs saved to: ", out_dir,
  "\nReview acceptance rates, trace plots, coefficient uncertainty, ",
  "and prediction plots before interpreting the results."
)
