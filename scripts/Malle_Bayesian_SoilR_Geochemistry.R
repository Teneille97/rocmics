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
})

set.seed(20261026)

integrated_file <- here::here(
  "outputs", "integrated_geochemistry_for_bayesian.csv"
)
mcmc_dir <- here::here("mod_runs")
out_dir <- here::here("outputs", "bayesian_persistence")
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
  MAOM_age = c("soil_pH", "Fe_ox_over_DCB", "Si_ox_over_DCB"),
  oPOM_age = c("soil_EC", "CEC_Ca", "GMD_h2o"),
  fPOM_age = c(
    "Tessier_H2O2_Na_mg_g",
    "Mn_ox_over_DCB",
    "Tessier_Acetic_Mg_mg_g"
  ),
  system_transit_time = c("soil_pH", "Fe_ox_over_DCB", "CEC_Ca")
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

# Check all expected MCMC files before starting the computational loop.
missing_mcmc_files <- file.path(mcmc_dir, unname(mcmc_files))
missing_mcmc_files <- missing_mcmc_files[!file.exists(missing_mcmc_files)]

if (length(missing_mcmc_files) > 0L) {
  stop(
    "The following MCMC files are missing:\n",
    paste(missing_mcmc_files, collapse = "\n")
  )
}

soilR_draws <- purrr::map_dfr(names(mcmc_files), function(trt) {
  message("Calculating persistence draws for ", trt)
  
  pars <- get_par_draws(file.path(mcmc_dir, mcmc_files[[trt]]))
  
  purrr::map_dfr(seq_len(nrow(pars)), function(i) {
    one_persistence_draw(as.numeric(pars[i, ]))
  }) %>%
    mutate(
      Treatment_model = trt,
      draw_id = row_number(),
      .before = 1
    )
}) %>%
  filter(
    if_all(
      c(MAOM_age, oPOM_age, fPOM_age, system_transit_time),
      is.finite
    )
  )

if (nrow(soilR_draws) == 0L) {
  stop(
    "No valid SoilR persistence draws were produced. ",
    "Check the MCMC parameter columns and SoilR density outputs."
  )
}

readr::write_csv(
  soilR_draws,
  file.path(out_dir, "soilR_persistence_posterior_draws.csv")
)

response_names <- c(
  "MAOM_age", "oPOM_age", "fPOM_age", "system_transit_time"
)

soilR_summary <- soilR_draws %>%
  pivot_longer(
    all_of(response_names),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  group_by(Treatment_model, outcome) %>%
  summarise(
    n_draws = n(),
    mean = mean(value),
    median = median(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975),
    .groups = "drop"
  )

readr::write_csv(
  soilR_summary,
  file.path(out_dir, "soilR_persistence_summary.csv")
)

# ---------------------------------------------------------------------------
# 4. Build the treatment-level regression dataset
# ---------------------------------------------------------------------------

# Posterior means and SDs summarise uncertainty in each SoilR-derived outcome.
outcome_summary <- soilR_draws %>%
  pivot_longer(
    all_of(response_names),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  group_by(Treatment_model, outcome) %>%
  summarise(
    y_mean = mean(value),
    y_se = sd(value),
    .groups = "drop"
  )

model_data_long <- geochem_summary %>%
  inner_join(outcome_summary, by = "Treatment_model")

# Create predictor means and standard errors expected by fit_one().
# Predictor SE = between-plot SD / sqrt(number of available plots).
for (p in all_predictors) {
  model_data_long[[p]] <- model_data_long[[paste0(p, "__mean")]]
  model_data_long[[paste0(p, "_se")]] <-
    model_data_long[[paste0(p, "__sd")]] /
    sqrt(model_data_long[[paste0(p, "__n")]])
}

readr::write_csv(
  model_data_long,
  file.path(out_dir, "treatment_level_model_data_long.csv")
)

# Validate response/treatment layout before fitting any models.
model_counts <- model_data_long %>%
  count(outcome, name = "n_treatments")
print(model_counts, n = Inf)

if (!setequal(model_counts$outcome, response_names) ||
    any(model_counts$n_treatments != 5L)) {
  stop("Expected exactly five treatment-level observations for every outcome.")
}

# Inspect the treatment-level predictor means and SEs used by the regressions.
print(
  model_data_long %>%
    select(
      Treatment_model,
      outcome,
      all_of(all_predictors),
      all_of(paste0(all_predictors, "_se"))
    ),
  n = Inf
)

# Create a backup folder
dir.create("outputs/bayesian_persistence/session_backup",
           recursive = TRUE, showWarnings = FALSE)

# Save all objects in the current R environment
save.image(
  file = "outputs/bayesian_persistence/session_backup/session_before_restart.RData"
)
saveRDS(
  model_data_long,
  "outputs/bayesian_persistence/session_backup/model_data_long.rds"
)

saveRDS(
  geochem_summary,
  "outputs/bayesian_persistence/session_backup/geochem_summary.rds"
)
load("outputs/bayesian_persistence/session_backup/session_before_restart.RData")
# ---------------------------------------------------------------------------
# 5. Bayesian measurement-error regressions
# ---------------------------------------------------------------------------

# Standardise each predictor and its SE on the same scale within each model.
# The response uncertainty is represented by the posterior SD of the SoilR
# persistence metric. These models are exploratory because n = 5 per model


priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 1), class = "Intercept"),
  prior(exponential(1), class = "sigma")
)

fit_one <- function(outcome, predictors) {
  dat <- model_data_long %>%
    filter(.data$outcome == .env$outcome)
  
  message(
    "Fitting ", outcome,
    ": ", nrow(dat), " rows; treatments = ",
    paste(dat$Treatment_model, collapse = ", ")
  )
  
  stopifnot(nrow(dat) == 5L)
  stopifnot(dplyr::n_distinct(dat$Treatment_model) == 5L)
  
  for (p in predictors) {
    x <- dat[[p]]
    x_sd <- sd(x, na.rm = TRUE)
    
    if (!is.finite(x_sd) || x_sd == 0) {
      stop(p, " has no variation across the analysis dataset for ", outcome, ".")
    }
    
    dat[[paste0(p, "_se")]] <- dat[[paste0(p, "_se")]] / x_sd
    dat[[p]] <- (x - mean(x, na.rm = TRUE)) / x_sd
  }
  
  needed <- c("y_mean", "y_se", predictors, paste0(predictors, "_se"))
  
  dat <- dat %>%
    select(Treatment_model, all_of(needed)) %>%
    filter(if_all(all_of(needed), ~ is.finite(.x)))
  
  if (nrow(dat) != 5L) {
    stop(
      outcome, " has ", nrow(dat),
      " complete treatment rows; expected five. Check missing values."
    )
  }
  
  if (any(dat$y_se <= 0)) {
    stop(outcome, " has non-positive response SDs. Check SoilR posterior draws.")
  }
  
  predictor_se_matrix <- as.matrix(dat[paste0(predictors, "_se")])
  if (any(predictor_se_matrix < 0)) {
    stop(outcome, " has negative predictor standard errors.")
  }
  
  f <- as.formula(paste(
    "y_mean | se(y_se, sigma = TRUE) ~",
    paste(
      vapply(
        predictors,
        function(p) paste0("me(", p, ", ", p, "_se)"),
        character(1)
      ),
      collapse = " + "
    )
  ))
  
  message("\nFitting ", outcome, ": ", deparse(f))
  
  n_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
  
  fit <- brms::brm(
    formula = f,
    data = dat,
    family = gaussian(),
    prior = priors,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    cores = min(4L, n_cores),
    seed = 20261026,
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    refresh = 500,
    save_pars = brms::save_pars(all = TRUE)
  )
  
  saveRDS(fit, file.path(out_dir, paste0("brms_", outcome, ".rds")))
  capture.output(
    print(summary(fit)),
    file = file.path(out_dir, paste0("brms_", outcome, "_summary.txt"))
  )
  
  d <- posterior::as_draws_df(fit)
  bn <- grep("^b_", names(d), value = TRUE)
  
  coef <- purrr::map_dfr(bn, function(nm) {
    x <- d[[nm]]
    tibble(
      outcome = outcome,
      parameter = nm,
      mean = mean(x),
      median = median(x),
      sd = sd(x),
      q2.5 = unname(quantile(x, 0.025)),
      q97.5 = unname(quantile(x, 0.975)),
      Pr_positive = mean(x > 0),
      Pr_negative = mean(x < 0)
    )
  })
  
  readr::write_csv(
    coef,
    file.path(out_dir, paste0("coefficients_", outcome, ".csv"))
  )
  
  pp <- brms::pp_check(fit, type = "dens_overlay", ndraws = 100)
  ggsave(
    file.path(out_dir, paste0("pp_check_", outcome, ".png")),
    pp, width = 7, height = 5, dpi = 300
  )
  
  tr <- brms::mcmc_plot(fit, type = "trace")
  ggsave(
    file.path(out_dir, paste0("trace_", outcome, ".png")),
    tr, width = 9, height = 7, dpi = 300
  )
  
  list(fit = fit, data = dat, coefficients = coef)
}

fits <- list()
for (nm in names(model_specs)) {
  fits[[nm]] <- fit_one(nm, model_specs[[nm]])
}

# ---------------------------------------------------------------------------
# 6. Optional projpred sensitivity analysis (off by default)
# ---------------------------------------------------------------------------

RUN_PROJPRED <- FALSE

if (RUN_PROJPRED) {
  if (!requireNamespace("projpred", quietly = TRUE)) {
    stop("Install projpred with install.packages('projpred').")
  }
  
  for (nm in names(fits)) {
    tryCatch({
      vs <- projpred::varsel(
        fits[[nm]]$fit,
        method = "forward",
        ndraws = 400,
        nterms_max = 3,
        refit_prj = FALSE
      )
      
      saveRDS(vs, file.path(out_dir, paste0("projpred_", nm, ".rds")))
      capture.output(
        print(summary(vs)),
        file = file.path(out_dir, paste0("projpred_", nm, "_summary.txt"))
      )
      ggsave(
        file.path(out_dir, paste0("projpred_", nm, ".png")),
        plot(vs), width = 7, height = 5, dpi = 300
      )
    }, error = function(e) {
      warning("projpred failed for ", nm, ": ", conditionMessage(e))
    })
  }
}

# ---------------------------------------------------------------------------
# 7. Sampling diagnostics
# ---------------------------------------------------------------------------

diagnostics <- purrr::map_dfr(names(fits), function(nm) {
  fit <- fits[[nm]]$fit
  s <- posterior::summarise_draws(posterior::as_draws_array(fit))
  np <- brms::nuts_params(fit)
  
  tibble(
    model = nm,
    n_treatments = nrow(fits[[nm]]$data),
    max_Rhat = if (all(is.na(s$rhat))) NA_real_ else max(s$rhat, na.rm = TRUE),
    min_bulk_ESS = if (all(is.na(s$ess_bulk))) NA_real_ else min(s$ess_bulk, na.rm = TRUE),
    min_tail_ESS = if (all(is.na(s$ess_tail))) NA_real_ else min(s$ess_tail, na.rm = TRUE),
    divergences = sum(
      np$Parameter == "divergent__" & np$Value == 1,
      na.rm = TRUE
    ),
    treedepth_hits = sum(
      np$Parameter == "treedepth__" & np$Value >= 12,
      na.rm = TRUE
    )
  )
})

readr::write_csv(diagnostics, file.path(out_dir, "model_diagnostics.csv"))
print(diagnostics, n = Inf)

message(
  "\nFinished. Outputs saved to: ", out_dir,
  "\nInspect soilR_persistence_summary.csv, geochemistry_treatment_summary.csv, ",
  "treatment_level_model_data_long.csv, coefficients_*.csv, model_diagnostics.csv, ",
  "pp_check_*.png, and trace_*.png.",
  "\nAge densities are integrated over 0-500 years. Extend age_grid if the ",
  "density tails are not negligible at 500 years.",
  "\nReminder: these are highly exploratory regressions with five treatment-level ",
  "observations per model."
)
