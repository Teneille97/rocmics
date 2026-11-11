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

# =============================================================================
# PLOT-ONLY MODE: LOAD SAVED POSTERIOR SUMMARIES
# =============================================================================
# This script deliberately does NOT compile Stan or fit/sample any models.
# Required inputs are the saved CSV summaries from the full modelling script
# and the integrated dataset used for the raw-data violin/dose plots.

primary_csv <- file.path(
  bayesian_treatment_out_dir,
  "Bayesian_primary_basalt_vs_Control_Lime.csv"
)
pairwise_csv <- file.path(
  bayesian_treatment_out_dir,
  "Bayesian_all_pairwise_treatment_contrasts.csv"
)
means_csv <- file.path(
  bayesian_treatment_out_dir,
  "Bayesian_treatment_posterior_means.csv"
)

required_saved_files <- c(primary_csv, pairwise_csv, means_csv)
missing_saved_files <- required_saved_files[!file.exists(required_saved_files)]
if (length(missing_saved_files) > 0) {
  stop(
    "Missing saved Bayesian output file(s):\n",
    paste(missing_saved_files, collapse = "\n"),
    "\n\nRun the full Bayesian modelling script once to create these files."
  )
}

primary_contrasts <- readr::read_csv(primary_csv, show_col_types = FALSE)
pairwise_contrasts <- readr::read_csv(pairwise_csv, show_col_types = FALSE)
bayesian_treatment_means_all <- readr::read_csv(means_csv, show_col_types = FALSE)

cat("\nLoaded saved posterior summaries. No models will be fitted.\n")


# =============================================================================
# 5. DEFINE PLOT GROUPS
# =============================================================================
#
# Each group produces:
#   - one primary forest plot
#   - one all-pairwise forest plot
#   - one violin plot
#   - one dose-response plot, where dose-response data are available
#
# Carbon stocks and radiocarbon signatures are combined so that each SOM
# fraction can be compared across C stocks and d14C.
#
# The four Tessier extractants are kept separate.
# =============================================================================

plot_groups <- list(
  
  Soil_Chemistry = c(
    "soil_pH",
    "soil_EC"
  ),
  
  Aggregate_Stability = c(
    "GMD_h2o",
    "D_h2o",
    "MWD_h2o",
    "GMD_oxalate_over_h2o",
    "D_oxalate_over_h2o",
    "MWD_oxalate_over_h2o"
  ),
  
  CEC_Cations = c(
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
  
  Oxalate_DCB_Ratios = c(
    "Fe_ox_over_DCB",
    "Al_ox_over_DCB",
    "Mn_ox_over_DCB",
    "Si_ox_over_DCB",
    "Fe_ox_over_Si_ox",
    "Al_ox_over_Si_ox"
  ),
  
  Tessier_Carbonate = c(
    "Tessier_Acetic_Ca_mg_g",
    "Tessier_Acetic_K_mg_g",
    "Tessier_Acetic_Na_mg_g",
    "Tessier_Acetic_Mg_mg_g"
  ),
  
  Tessier_Exchangeable = c(
    "Tessier_NH4OAc_Ca_mg_g",
    "Tessier_NH4OAc_K_mg_g",
    "Tessier_NH4OAc_Na_mg_g",
    "Tessier_NH4OAc_Mg_mg_g"
  ),
  
  Tessier_Oxidizable = c(
    "Tessier_H2O2_Fe_mg_g",
    "Tessier_H2O2_Al_mg_g",
    "Tessier_H2O2_Mn_mg_g",
    "Tessier_H2O2_Si_mg_g"
  ),
  
  Tessier_Reducible = c(
    "Tessier_Hydroxylamine_Fe_mg_g",
    "Tessier_Hydroxylamine_Al_mg_g",
    "Tessier_Hydroxylamine_Mn_mg_g",
    "Tessier_Hydroxylamine_Si_mg_g"
  ),
  
  FTIR_Spectral = c(
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
  
  Carbon_and_Radiocarbon = c(
    "C_Bulk",
    "C_MAOM",
    "C_oPOM",
    "C_fPOM",
    "d14C_2026_Bulk",
    "d14C_2026_MAOM",
    "d14C_2026_oPOM",
    "d14C_2026_fPOM"
  )
)


# Retain only response variables actually present in the integrated dataset.

plot_groups_available <- lapply(
  plot_groups,
  function(x) {
    intersect(x, names(bayesian_treatment_data))
  }
)

# Remove empty groups.

plot_groups_available <- plot_groups_available[
  lengths(plot_groups_available) > 0
]


# Report missing variables.

all_requested_responses <- unique(
  unlist(plot_groups, use.names = FALSE)
)

all_available_responses <- unique(
  unlist(plot_groups_available, use.names = FALSE)
)

missing_responses <- setdiff(
  all_requested_responses,
  names(bayesian_treatment_data)
)

if (length(missing_responses) > 0) {
  cat("\nWARNING: The following requested responses are absent:\n")
  print(missing_responses)
}


cat("\nAvailable plot groups:\n")

for (g in names(plot_groups_available)) {
  cat(
    "\n",
    g,
    ": ",
    length(plot_groups_available[[g]]),
    " responses\n",
    sep = ""
  )
  print(plot_groups_available[[g]])
}



# =============================================================================
# PLOT LABELS, UNITS AND TREATMENT COLOURS
# =============================================================================

# Control is neutral grey; the four amendments use the requested colours.
treatment_colours <- c(
  "Control" = "#777777",
  "Lime" = "#78A641",
  "HB 50" = "#C65A9E",
  "EG 50" = "#D99A16",
  "BD 50" = "#2878B5"
)

treatment_display <- c(
  "Control" = "Control",
  "Lime" = "Lime",
  "HB 50" = "Hühnerberg",
  "EG 50" = "Eifelgold",
  "BD 50" = "Bolsdorfer"
)

group_display <- c(
  Soil_Chemistry = "Soil chemistry",
  Aggregate_Stability = "Aggregate stability",
  CEC_Cations = "CEC and exchangeable cations",
  Oxalate_DCB_Ratios = "Oxalate and DCB ratios",
  Tessier_Carbonate = "Tessier carbonate fraction",
  Tessier_Exchangeable = "Tessier exchangeable fraction",
  Tessier_Oxidizable = "Tessier oxidizable fraction",
  Tessier_Reducible = "Tessier reducible fraction",
  FTIR_Spectral = "FTIR spectral changes",
  Carbon_and_Radiocarbon = "Carbon and radiocarbon"
)

# Units are shown in each facet heading because several panels in a grouped
# figure have different units. This avoids assigning one misleading unit to
# the entire shared axis.
response_display <- c(
  soil_pH = "Soil pH (no units)",
  soil_EC = "Soil EC (µS cm⁻¹)",
  GMD_h2o = "GMD-H₂O (mm)",
  D_h2o = "D-H₂O (no units)",
  MWD_h2o = "MWD-H₂O (mm)",
  GMD_oxalate_over_h2o = "GMD-Oxalate / GMD-H₂O (no units)",
  D_oxalate_over_h2o = "D-Oxalate / D-H₂O (no units)",
  MWD_oxalate_over_h2o = "MWD-Oxalate / MWD-H₂O (no units)",
  CEC_CEC = "CEC (cmolc kg⁻¹)",
  CEC_Base_saturation = "Base saturation (%)",
  CEC_H = "Exchangeable H (cmolc kg⁻¹)",
  CEC_Ca = "Exchangeable Ca (cmolc kg⁻¹)",
  CEC_K = "Exchangeable K (cmolc kg⁻¹)",
  CEC_Mg = "Exchangeable Mg (cmolc kg⁻¹)",
  CEC_Na = "Exchangeable Na (cmolc kg⁻¹)",
  CEC_Al = "Exchangeable Al (cmolc kg⁻¹)",
  CEC_Fe = "Exchangeable Fe (cmolc kg⁻¹)",
  CEC_Mn = "Exchangeable Mn (cmolc kg⁻¹)",
  Fe_ox_over_DCB = "Fe-ox / Fe-DCB (no units)",
  Al_ox_over_DCB = "Al-ox / Al-DCB (no units)",
  Mn_ox_over_DCB = "Mn-ox / Mn-DCB (no units)",
  Si_ox_over_DCB = "Si-ox / Si-DCB (no units)",
  Fe_ox_over_Si_ox = "Fe-ox / Si-ox (no units)",
  Al_ox_over_Si_ox = "Al-ox / Si-ox (no units)",
  Tessier_Acetic_Ca_mg_g = "Ca-Acetic (mg g⁻¹)",
  Tessier_Acetic_K_mg_g = "K-Acetic (mg g⁻¹)",
  Tessier_Acetic_Na_mg_g = "Na-Acetic (mg g⁻¹)",
  Tessier_Acetic_Mg_mg_g = "Mg-Acetic (mg g⁻¹)",
  Tessier_NH4OAc_Ca_mg_g = "Ca-NH₄OAc (mg g⁻¹)",
  Tessier_NH4OAc_K_mg_g = "K-NH₄OAc (mg g⁻¹)",
  Tessier_NH4OAc_Na_mg_g = "Na-NH₄OAc (mg g⁻¹)",
  Tessier_NH4OAc_Mg_mg_g = "Mg-NH₄OAc (mg g⁻¹)",
  Tessier_H2O2_Fe_mg_g = "Fe-H₂O₂ (mg g⁻¹)",
  Tessier_H2O2_Al_mg_g = "Al-H₂O₂ (mg g⁻¹)",
  Tessier_H2O2_Mn_mg_g = "Mn-H₂O₂ (mg g⁻¹)",
  Tessier_H2O2_Si_mg_g = "Si-H₂O₂ (mg g⁻¹)",
  Tessier_Hydroxylamine_Fe_mg_g = "Fe-Hydroxylamine (mg g⁻¹)",
  Tessier_Hydroxylamine_Al_mg_g = "Al-Hydroxylamine (mg g⁻¹)",
  Tessier_Hydroxylamine_Mn_mg_g = "Mn-Hydroxylamine (mg g⁻¹)",
  Tessier_Hydroxylamine_Si_mg_g = "Si-Hydroxylamine (mg g⁻¹)",
  FTIR_delta_aromatic = "Aromatic (no units)",
  FTIR_delta_Carbonate_bend = "Carbonate bend (no units)",
  FTIR_delta_carbonate_stretch = "Carbonate stretch (no units)",
  FTIR_delta_carbonyl = "Carbonyl (no units)",
  FTIR_delta_Aliphatic_total = "Aliphatic (no units)",
  FTIR_delta_COO_total = "Carboxylic (no units)",
  FTIR_delta_polysaccharide = "Polysaccharide (no units)",
  FTIR_delta_poorly_crystalline = "Amorphous (no units)",
  FTIR_delta_Silicate = "Silicate (no units)",
  FTIR_delta_total_area = "Total (no units)",
  C_Bulk = "Bulk carbon (%)",
  C_MAOM = "MAOM carbon (%)",
  C_oPOM = "oPOM carbon (%)",
  C_fPOM = "fPOM carbon (%)",
  d14C_2026_Bulk = "Bulk Δ¹⁴C (‰)",
  d14C_2026_MAOM = "MAOM Δ¹⁴C (‰)",
  d14C_2026_oPOM = "oPOM Δ¹⁴C (‰)",
  d14C_2026_fPOM = "fPOM Δ¹⁴C (‰)"
)

display_response <- function(x) {
  out <- unname(response_display[x])
  missing <- is.na(out)
  out[missing] <- gsub("_", " ", x[missing])
  out
}

display_group <- function(x) {
  out <- unname(group_display[x])
  ifelse(is.na(out), gsub("_", " ", x), out)
}

display_comparison <- function(x) {
  x <- gsub("HB 50", "Hühnerberg", x, fixed = TRUE)
  x <- gsub("EG 50", "Eifelgold", x, fixed = TRUE)
  x <- gsub("BD 50", "Bolsdorfer", x, fixed = TRUE)
  x <- gsub("minus", "−", x, fixed = TRUE)
  x
}


# =============================================================================
# 6. CHECK AND PREPARE POSTERIOR CONTRAST TABLES
# =============================================================================

required_primary_columns <- c(
  "Response",
  "Treatment_1",
  "Treatment_2",
  "Difference_original",
  "Lower_95_original",
  "Upper_95_original"
)

required_pairwise_columns <- required_primary_columns


missing_primary_columns <- setdiff(
  required_primary_columns,
  names(primary_contrasts)
)

missing_pairwise_columns <- setdiff(
  required_pairwise_columns,
  names(pairwise_contrasts)
)

if (length(missing_primary_columns) > 0) {
  stop(
    "Missing columns in primary contrast CSV: ",
    paste(missing_primary_columns, collapse = ", ")
  )
}

if (length(missing_pairwise_columns) > 0) {
  stop(
    "Missing columns in pairwise contrast CSV: ",
    paste(missing_pairwise_columns, collapse = ", ")
  )
}


# Primary comparisons must all be treatment minus Control.

primary_contrasts <- primary_contrasts %>%
  filter(
    Treatment_2 == "Control",
    Treatment_1 %in% c("Lime", "HB 50", "EG 50", "BD 50")
  ) %>%
  mutate(
    Treatment_1 = factor(
      Treatment_1,
      levels = c("Lime", "HB 50", "EG 50", "BD 50")
    ),
    Comparison = paste(as.character(Treatment_1), "minus Control")
  )


# Keep the original comparison labels for the all-pairwise results.

pairwise_contrasts <- pairwise_contrasts %>%
  mutate(
    Comparison = if ("Comparison" %in% names(.)) {
      Comparison
    } else {
      paste(Treatment_1, "minus", Treatment_2)
    }
  )


# =============================================================================
# 7. CONVERT RAW DATA TO LONG FORMAT FOR VIOLIN PLOTS
# =============================================================================

# IMPORTANT:
#   Build a plain character vector of existing response columns.
#   This avoids the pivot_longer()/any_of() error.

plot_responses <- unique(
  unlist(plot_groups_available, use.names = FALSE)
)

plot_responses <- intersect(
  plot_responses,
  names(bayesian_treatment_data)
)

if (length(plot_responses) == 0) {
  stop("No plot response variables were found in bayesian_treatment_data.")
}


raw_data_plot <- bayesian_treatment_data %>%
  pivot_longer(
    cols = all_of(plot_responses),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  drop_na(Value) %>%
  mutate(
    Treatment_PCA = factor(
      Treatment_PCA,
      levels = bayesian_treatment_levels
    ),
    Response = as.character(Response)
  )


# =============================================================================
# HELPER: DISPLAY GROUPED FOREST PLOTS IN ORIGINAL UNITS
# =============================================================================

make_grouped_forest_plot <- function(
    data,
    responses,
    group_name,
    comparison_type = "Primary") {
  
  plot_data <- data %>%
    filter(Response %in% responses) %>%
    filter(
      is.finite(Difference_original),
      is.finite(Lower_95_original),
      is.finite(Upper_95_original)
    ) %>%
    mutate(
      Response = factor(
        Response,
        levels = responses,
        labels = display_response(responses)
      )
    )
  
  if (nrow(plot_data) == 0) {
    message(
      "Skipping forest plot: ",
      group_name,
      " (no original-scale contrast data)"
    )
    return(NULL)
  }
  
  # Order comparisons consistently.
  
  if (comparison_type == "Primary") {
    
    comparison_levels <- c(
      "Lime minus Control",
      "HB 50 minus Control",
      "EG 50 minus Control",
      "BD 50 minus Control"
    )
    
  } else {
    
    comparison_levels <- unique(plot_data$Comparison)
    
  }
  
  plot_data <- plot_data %>%
    mutate(
      Comparison = factor(
        display_comparison(as.character(Comparison)),
        levels = rev(display_comparison(comparison_levels))
      ),
      Treatment_colour_key = if_else(
        as.character(Treatment_1) == "Control",
        as.character(Treatment_2),
        as.character(Treatment_1)
      )
    )
  
  p <- ggplot(
    plot_data,
    aes(
      x = Difference_original,
      y = Comparison,
      colour = Treatment_colour_key
    )
  ) +
    
    # Zero reference line.
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.5
    ) +
    
    # 95% credible intervals.
    geom_errorbar(
      aes(
        xmin = Lower_95_original,
        xmax = Upper_95_original
      ),
      orientation = "y",
      width = 0.18,
      linewidth = 0.6
    ) +
    
    # Posterior mean differences.
    geom_point(size = 2.4) +
    scale_colour_manual(
      values = treatment_colours,
      breaks = c("Lime", "HB 50", "EG 50", "BD 50"),
      labels = treatment_display[c("Lime", "HB 50", "EG 50", "BD 50")],
      name = "Treatment"
    ) +
    
    # Independent horizontal scale for each response.
    facet_wrap(
      ~ Response,
      ncol = 2,
      scales = "free_x"
    ) +
    
    labs(
      title = paste(
        ifelse(comparison_type == "Primary",
               "Primary treatment contrasts:",
               "All pairwise treatment contrasts:"),
        display_group(group_name)
      ),
      subtitle = "Points show posterior mean differences; lines show 95% credible intervals",
      x = "Treatment difference (units specified in each panel heading)",
      y = NULL,
      colour = "Treatment"
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey40"
      ),
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9),
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      plot.subtitle = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(1, "lines")
    )
  
  n_responses <- length(unique(plot_data$Response))
  
  plot_height <- max(
    4.5,
    ceiling(n_responses / 2) * 2.8
  )
  
  plot_height <- max(4.5, ceiling(n_responses / 2) * 2.8)
  plot_prefix <- if (comparison_type == "Primary") "Primary_Forest_" else "All_Pairwise_Forest_"
  filename <- file.path(
    bayesian_treatment_plots_dir,
    paste0(plot_prefix, group_name, ".png")
  )
  ggsave(filename, p, width = 13, height = plot_height, dpi = 300, bg = "white")
  ggsave(sub("\\.png$", ".pdf", filename), p, width = 13, height = plot_height, bg = "white")
  message("Saved: ", filename)

  print(p)
  
  return(p)
}


# =============================================================================
# 9. HELPER: SAVE GROUPED VIOLIN PLOTS
# =============================================================================

make_grouped_violin_plot <- function(
    data,
    responses,
    group_name,
    output_dir = bayesian_treatment_plots_dir) {
  
  plot_data <- data %>%
    filter(Response %in% responses) %>%
    filter(is.finite(Value)) %>%
    mutate(
      Response = factor(
        Response,
        levels = responses,
        labels = display_response(responses)
      ),
      Treatment_PCA = factor(
        Treatment_PCA,
        levels = bayesian_treatment_levels
      )
    )
  
  if (nrow(plot_data) == 0) {
    message(
      "Skipping violin plot: ",
      group_name,
      " (no raw data)"
    )
    return(NULL)
  }
  
  p <- ggplot(
    plot_data,
    aes(
      x = Treatment_PCA,
      y = Value,
      fill = Treatment_PCA
    )
  ) +
    
    geom_violin(
      alpha = 0.65,
      trim = FALSE,
      scale = "width",
      na.rm = TRUE
    ) +
    
    geom_boxplot(
      width = 0.12,
      colour = "black",
      alpha = 0.8,
      outlier.shape = NA,
      linewidth = 0.35,
      na.rm = TRUE
    ) +
    
    facet_wrap(
      ~ Response,
      scales = "free_y",
      ncol = 2
    ) +
    
    scale_x_discrete(labels = treatment_display) +
    
    scale_fill_manual(
      values = treatment_colours,
      breaks = bayesian_treatment_levels,
      labels = treatment_display[bayesian_treatment_levels],
      drop = FALSE
    ) +


        labs(
      title = paste("Treatment distributions:", display_group(group_name)),
      x = NULL,
      y = "Measured value (units specified in each panel heading)"
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      legend.position = "none",
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey40"
      ),
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      axis.text.x = element_text(
        angle = 35,
        hjust = 1,
        size = 8
      ),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(1, "lines")
    )
  
  n_responses <- length(unique(plot_data$Response))
  
  plot_height <- max(
    5,
    ceiling(n_responses / 2) * 3
  )
  
  filename <- file.path(
    output_dir,
    paste0(
      "Violin_",
      group_name,
      ".png"
    )
  )
  
  ggsave(
    filename = filename,
    plot = p,
    width = 13,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )
  
  message("Saved: ", filename)
  
  return(p)
}


# =============================================================================
# 10. GENERATE PRIMARY AND PAIRWISE FOREST PLOTS
# =============================================================================

# Primary treatment versus Control forest plots.

for (group_name in names(plot_groups_available)) {
  
  make_grouped_forest_plot(
    data = primary_contrasts,
    responses = plot_groups_available[[group_name]],
    group_name = group_name,
    comparison_type = "Primary"
  )
}


# All-pairwise treatment forest plots.

for (group_name in names(plot_groups_available)) {
  
  make_grouped_forest_plot(
    data = pairwise_contrasts,
    responses = plot_groups_available[[group_name]],
    group_name = group_name,
    comparison_type = "All_Pairwise"
  )
}

# =============================================================================
# 12. GENERATE GROUPED VIOLIN PLOTS
# =============================================================================

cat("\n\nGenerating grouped violin plots...\n")

violin_plots <- list()

for (group_name in names(plot_groups_available)) {
  
  responses <- plot_groups_available[[group_name]]
  
  violin_plots[[group_name]] <-
    make_grouped_violin_plot(
      data = raw_data_plot,
      responses = responses,
      group_name = group_name
    )
}


# =============================================================================
# 13. HÜHNERBERG DOSE-RESPONSE DATA
# =============================================================================
#
# This section uses the original integrated dataset rather than the
# five-group dataset, so that the full Hühnerberg dose series is retained.
#
# Control is assigned a dose of zero.
# Dose is extracted from the treatment label.
# =============================================================================

cat("\n\nPreparing Hühnerberg dose-response data...\n")


hb_dose_data <- integrated %>%
  filter(
    Treatment == "Control" |
      grepl("^Huhnerberg", Treatment) |
      grepl("^HB_", Treatment) |
      grepl("^HB[0-9]", Treatment)
  ) %>%
  mutate(
    Dose_t_ha = case_when(
      Treatment == "Control" ~ 0,
      TRUE ~ readr::parse_number(Treatment)
    )
  ) %>%
  filter(!is.na(Dose_t_ha)) %>%
  pivot_longer(
    cols = all_of(
      intersect(
        plot_responses,
        names(.)
      )
    ),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  drop_na(Value) %>%
  mutate(
    Response = as.character(Response)
  )


# =============================================================================
# 14. HELPER: SAVE GROUPED DOSE-RESPONSE PLOTS
# =============================================================================

make_grouped_dose_plot <- function(
    data,
    responses,
    group_name,
    output_dir = bayesian_treatment_plots_dir) {
  
  plot_data <- data %>%
    filter(Response %in% responses) %>%
    filter(
      is.finite(Dose_t_ha),
      is.finite(Value)
    ) %>%
    mutate(
      Response = factor(
        Response,
        levels = responses,
        labels = display_response(responses)
      )
    )
  
  if (nrow(plot_data) == 0) {
    message(
      "Skipping dose-response plot: ",
      group_name,
      " (no observations)"
    )
    return(NULL)
  }
  
  # Summarise replicate values at each dose.
  # Points and error bars represent mean +/- SD.
  
  dose_summary <- plot_data %>%
    group_by(Response, Dose_t_ha) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      n = sum(!is.na(Value)),
      SE = SD / sqrt(n),
      .groups = "drop"
    )
  
  p <- ggplot(
    dose_summary,
    aes(
      x = Dose_t_ha,
      y = Mean
    )
  ) +
    
    geom_line(
      linewidth = 0.7,
      colour = "darkblue"
    ) +
    
    geom_point(
      size = 2.2,
      colour = "darkblue"
    ) +
    
    geom_errorbar(
      aes(
        ymin = Mean - SD,
        ymax = Mean + SD
      ),
      width = 0.5,
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    
    facet_wrap(
      ~ Response,
      scales = "free_y",
      ncol = 2
    ) +
    
    scale_x_continuous(
      breaks = sort(unique(dose_summary$Dose_t_ha))
    ) +
    
    labs(
      title = paste("Hühnerberg dose response:", display_group(group_name)),
      subtitle = "Mean ± SD across replicate observations",
      x = "Hühnerberg application rate (t ha⁻¹)",
      y = "Measured value (units specified in each panel heading)"
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey40"
      ),
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(1, "lines")
    )
  
  n_responses <- length(unique(plot_data$Response))
  
  plot_height <- max(
    5,
    ceiling(n_responses / 2) * 3
  )
  
  filename <- file.path(
    output_dir,
    paste0(
      "Huhnerberg_Dose_Response_",
      group_name,
      ".png"
    )
  )
  
  ggsave(
    filename = filename,
    plot = p,
    width = 13,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )
  
  message("Saved: ", filename)
  
  return(p)
}


# =============================================================================
# 15. GENERATE GROUPED DOSE-RESPONSE PLOTS
# =============================================================================

cat("\n\nGenerating grouped Hühnerberg dose-response plots...\n")

dose_response_plots <- list()

for (group_name in names(plot_groups_available)) {
  
  responses <- plot_groups_available[[group_name]]
  
  dose_response_plots[[group_name]] <-
    make_grouped_dose_plot(
      data = hb_dose_data,
      responses = responses,
      group_name = group_name
    )
}



cat("\nPLOT-ONLY RUN COMPLETE\n")
cat("Plots saved in: ", bayesian_treatment_plots_dir, "\n", sep = "")
