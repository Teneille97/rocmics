
# =============================================================================
# MALLE ERW: BAYESIAN GEOCHEMISTRY TREATMENT COMPARISONS
# PLOT-ONLY SCRIPT
#
# PURPOSE:
#   Load saved Bayesian geochemistry treatment comparisons and generate:
#
#   1. Primary treatment-vs-Control forest plots
#   2. All-pairwise treatment forest plots
#   3. Raw-data treatment violin plots
#   4. Hühnerberg dose-response plots
#   5. Posterior treatment density distributions (if posterior draws saved)
#   6. Posterior treatment-contrast density distributions (if contrast draws saved)
#
# IMPORTANT:
#   This script does NOT fit models, compile Stan or generate posterior draws.
#
#   Posterior density plots are generated only from actual saved posterior
#   draws. Summary statistics are not used to reconstruct distributions.
#
# TREATMENTS:
#   Control, Lime, HB 50, EG 50, BD 50
#
# =============================================================================


# =============================================================================
# 1. PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  
  library(tidyverse)
  library(here)
  
})


# =============================================================================
# 2. PATHS
# =============================================================================

outputs_dir <- here::here("outputs")

bayesian_treatment_out_dir <- file.path(
  outputs_dir,
  "Bayesian_treatment_comparisons"
)

bayesian_treatment_plots_dir <- file.path(
  bayesian_treatment_out_dir,
  "plots"
)

dir.create(
  bayesian_treatment_plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 3. LOAD INTEGRATED DATASET
# =============================================================================

integrated <- read_csv(
  file.path(
    outputs_dir,
    "integrated_geochemistry_for_bayesian.csv"
  ),
  na = c("", "NA", "#VALUE!"),
  show_col_types = FALSE
)


# =============================================================================
# 4. TREATMENT SELECTION
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


# =============================================================================
# 5. RESPONSE GROUPS
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


plot_groups_available <- lapply(
  plot_groups,
  function(x) intersect(x, names(bayesian_treatment_data))
)

plot_groups_available <- plot_groups_available[
  lengths(plot_groups_available) > 0
]


plot_responses <- unique(
  unlist(plot_groups_available, use.names = FALSE)
)


# =============================================================================
# 6. COLOURS AND DISPLAY LABELS
# =============================================================================

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
# 7. LOAD SAVED POSTERIOR SUMMARIES
# =============================================================================

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


required_saved_files <- c(
  primary_csv,
  pairwise_csv,
  means_csv
)


missing_saved_files <- required_saved_files[
  !file.exists(required_saved_files)
]


if (length(missing_saved_files) > 0) {
  
  stop(
    "Missing saved Bayesian output file(s):\n",
    paste(missing_saved_files, collapse = "\n")
  )
  
}


primary_contrasts <- read_csv(
  primary_csv,
  show_col_types = FALSE
)

pairwise_contrasts <- read_csv(
  pairwise_csv,
  show_col_types = FALSE
)

bayesian_treatment_means_all <- read_csv(
  means_csv,
  show_col_types = FALSE
)


cat("\nLoaded saved Bayesian summaries.\n")


# =============================================================================
# 8. PREPARE CONTRAST TABLES
# =============================================================================

required_columns <- c(
  "Response",
  "Treatment_1",
  "Treatment_2",
  "Difference_original",
  "Lower_95_original",
  "Upper_95_original"
)


for (column in required_columns) {
  
  if (!column %in% names(primary_contrasts)) {
    
    stop("Missing primary contrast column: ", column)
    
  }
  
  if (!column %in% names(pairwise_contrasts)) {
    
    stop("Missing pairwise contrast column: ", column)
    
  }
  
}


primary_contrasts <- primary_contrasts %>%
  
  filter(
    Treatment_2 == "Control",
    Treatment_1 %in% c("Lime", "HB 50", "EG 50", "BD 50")
  ) %>%
  
  mutate(
    Comparison = paste(Treatment_1, "minus Control")
  )


pairwise_contrasts <- pairwise_contrasts %>%
  
  mutate(
    
    Comparison = if ("Comparison" %in% names(.)) {
      Comparison
    } else {
      paste(Treatment_1, "minus", Treatment_2)
    }
    
  )


# =============================================================================
# 9. RAW DATA LONG FORMAT
# =============================================================================

raw_data_plot <- bayesian_treatment_data %>%
  
  pivot_longer(
    cols = all_of(plot_responses),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  
  filter(is.finite(Value)) %>%
  
  mutate(
    Response = as.character(Response),
    Treatment_PCA = factor(
      Treatment_PCA,
      levels = bayesian_treatment_levels
    )
  )


# =============================================================================
# 10. GROUPED FOREST PLOTS
# =============================================================================

make_grouped_forest_plot <- function(
    data,
    responses,
    group_name,
    comparison_type = "Primary"
) {
  
  plot_data <- data %>%
    
    filter(
      Response %in% responses,
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
  
  
  if (nrow(plot_data) == 0) return(NULL)
  
  
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
        display_comparison(Comparison),
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
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.5
    ) +
    
    geom_errorbar(
      aes(
        xmin = Lower_95_original,
        xmax = Upper_95_original
      ),
      orientation = "y",
      width = 0.18,
      linewidth = 0.6
    ) +
    
    geom_point(size = 2.4) +
    
    scale_colour_manual(
      values = treatment_colours,
      breaks = c("Lime", "HB 50", "EG 50", "BD 50"),
      labels = treatment_display[c("Lime", "HB 50", "EG 50", "BD 50")],
      name = "Treatment"
    ) +
    
    facet_wrap(
      ~ Response,
      ncol = 2,
      scales = "free_x"
    ) +
    
    labs(
      title = paste(
        ifelse(
          comparison_type == "Primary",
          "Primary treatment contrasts:",
          "All pairwise treatment contrasts:"
        ),
        display_group(group_name)
      ),
      subtitle = "Posterior mean differences and 95% credible intervals",
      x = "Treatment difference (units specified in each panel heading)",
      y = NULL
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey40"
      ),
      strip.text = element_text(face = "bold", size = 9),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(1, "lines")
    )
  
  
  height <- max(
    4.5,
    ceiling(length(unique(plot_data$Response)) / 2) * 2.8
  )
  
  
  prefix <- if (comparison_type == "Primary") {
    "Primary_Forest_"
  } else {
    "All_Pairwise_Forest_"
  }
  
  
  filename <- file.path(
    bayesian_treatment_plots_dir,
    paste0(prefix, group_name, ".png")
  )
  
  
  ggsave(filename, p, width = 13, height = height, dpi = 300, bg = "white")
  
  ggsave(
    sub("\\.png$", ".pdf", filename),
    p,
    width = 13,
    height = height,
    bg = "white"
  )
  
  
  print(p)
  
  p
  
}


# =============================================================================
# 11. RAW DATA VIOLIN PLOTS
# =============================================================================

make_grouped_violin_plot <- function(
    data,
    responses,
    group_name
) {
  
  plot_data <- data %>%
    
    filter(
      Response %in% responses,
      is.finite(Value)
    ) %>%
    
    mutate(
      Response = factor(
        Response,
        levels = responses,
        labels = display_response(responses)
      )
    )
  
  
  if (nrow(plot_data) == 0) return(NULL)
  
  
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
      scale = "width"
    ) +
    
    geom_boxplot(
      width = 0.12,
      colour = "black",
      alpha = 0.8,
      outlier.shape = NA,
      linewidth = 0.35
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
      title = paste(
        "Treatment distributions:",
        display_group(group_name)
      ),
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
      strip.text = element_text(face = "bold", size = 9),
      axis.text.x = element_text(
        angle = 35,
        hjust = 1,
        size = 8
      ),
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  
  height <- max(
    5,
    ceiling(length(responses) / 2) * 3
  )
  
  
  filename <- file.path(
    bayesian_treatment_plots_dir,
    paste0("Violin_", group_name, ".png")
  )
  
  
  ggsave(
    filename,
    p,
    width = 13,
    height = height,
    dpi = 300,
    bg = "white"
  )
  
  
  print(p)
  
  p
  
}


# =============================================================================
# 12. POSTERIOR TREATMENT DENSITY DISTRIBUTIONS
# =============================================================================
#
# These require actual posterior treatment draws.
#
# The treatment posterior means CSV is NOT a substitute for posterior draws.
#


posterior_draw_candidates <- c(
  
  file.path(
    bayesian_treatment_out_dir,
    "Bayesian_treatment_posterior_draws.csv"
  ),
  
  file.path(
    bayesian_treatment_out_dir,
    "posterior_treatment_draws.csv"
  ),
  
  file.path(
    bayesian_treatment_out_dir,
    "Bayesian_posterior_treatment_draws.csv"
  ),
  
  file.path(
    bayesian_treatment_out_dir,
    "posterior_draws.csv"
  )
  
)


posterior_draw_file <- posterior_draw_candidates[
  file.exists(posterior_draw_candidates)
][1]


make_posterior_treatment_density <- function(
    data,
    responses,
    group_name
) {
  
  plot_data <- data %>%
    
    filter(
      Response %in% responses,
      is.finite(Value)
    ) %>%
    
    mutate(
      Response = factor(
        Response,
        levels = responses,
        labels = display_response(responses)
      ),
      Treatment = factor(
        Treatment,
        levels = bayesian_treatment_levels
      )
    )
  
  
  if (nrow(plot_data) == 0) return(NULL)
  
  
  p <- ggplot(
    plot_data,
    aes(
      x = Value,
      fill = Treatment,
      colour = Treatment
    )
  ) +
    
    geom_density(
      alpha = 0.30,
      linewidth = 0.65,
      na.rm = TRUE
    ) +
    
    facet_wrap(
      ~ Response,
      scales = "free",
      ncol = 2
    ) +
    
    scale_fill_manual(
      values = treatment_colours,
      breaks = bayesian_treatment_levels,
      labels = treatment_display[bayesian_treatment_levels]
    ) +
    
    scale_colour_manual(
      values = treatment_colours,
      breaks = bayesian_treatment_levels,
      labels = treatment_display[bayesian_treatment_levels]
    ) +
    
    labs(
      title = paste(
        "Posterior treatment distributions:",
        display_group(group_name)
      ),
      subtitle = "Full Bayesian posterior distributions for each treatment",
      x = "Posterior treatment value (units specified in each panel heading)",
      y = "Density",
      fill = "Treatment",
      colour = "Treatment"
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey40"
      ),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  
  height <- max(
    5,
    ceiling(length(responses) / 2) * 3
  )
  
  
  filename <- file.path(
    bayesian_treatment_plots_dir,
    paste0(
      "Posterior_Treatment_Distributions_",
      group_name,
      ".png"
    )
  )
  
  
  ggsave(
    filename,
    p,
    width = 13,
    height = height,
    dpi = 300,
    bg = "white"
  )
  
  
  ggsave(
    sub("\\.png$", ".pdf", filename),
    p,
    width = 13,
    height = height,
    bg = "white"
  )
  
  
  print(p)
  
  p
  
}


posterior_treatment_density_plots <- list()


if (
  length(posterior_draw_file) > 0 &&
  !is.na(posterior_draw_file)
) {
  
  cat("\nLoading posterior treatment draws:\n")
  cat(posterior_draw_file, "\n")
  
  posterior_treatment_draws <- read_csv(
    posterior_draw_file,
    show_col_types = FALSE
  )
  
  
  required_draw_columns <- c(
    "Response",
    "Treatment",
    "Value"
  )
  
  
  missing_draw_columns <- setdiff(
    required_draw_columns,
    names(posterior_treatment_draws)
  )
  
  
  if (length(missing_draw_columns) == 0) {
    
    for (group_name in names(plot_groups_available)) {
      
      posterior_treatment_density_plots[[group_name]] <-
        
        make_posterior_treatment_density(
          data = posterior_treatment_draws,
          responses = plot_groups_available[[group_name]],
          group_name = group_name
        )
      
    }
    
  } else {
    
    warning(
      "Posterior treatment draws were found, but the expected columns ",
      "Response, Treatment and Value are missing. ",
      "Treatment density plots were skipped."
    )
    
  }
  
} else {
  
  cat(
    "\nNo saved posterior treatment draws found. ",
    "Posterior treatment density plots skipped.\n"
  )
  
}


# =============================================================================
# 13. POSTERIOR TREATMENT-CONTRAST DENSITY DISTRIBUTIONS
# =============================================================================
#
# These show treatment minus Control.
#
# A vertical zero line indicates no treatment difference.
#


contrast_draw_candidates <- c(
  
  file.path(
    bayesian_treatment_out_dir,
    "Bayesian_primary_contrast_draws.csv"
  ),
  
  file.path(
    bayesian_treatment_out_dir,
    "Bayesian_primary_posterior_contrast_draws.csv"
  ),
  
  file.path(
    bayesian_treatment_out_dir,
    "primary_contrast_draws.csv"
  ),
  
  file.path(
    bayesian_treatment_out_dir,
    "posterior_contrast_draws.csv"
  ),
  
  file.path(
    bayesian_treatment_out_dir,
    "Bayesian_contrast_draws.csv"
  )
  
)


contrast_draw_file <- contrast_draw_candidates[
  file.exists(contrast_draw_candidates)
][1]


make_contrast_density_plot <- function(
    data,
    responses,
    group_name
) {
  
  plot_data <- data %>%
    
    filter(
      Response %in% responses,
      is.finite(Difference)
    ) %>%
    
    mutate(
      Response = factor(
        Response,
        levels = responses,
        labels = display_response(responses)
      ),
      Treatment = factor(
        Treatment,
        levels = c("Lime", "HB 50", "EG 50", "BD 50")
      )
    )
  
  
  if (nrow(plot_data) == 0) return(NULL)
  
  
  p <- ggplot(
    plot_data,
    aes(
      x = Difference,
      fill = Treatment,
      colour = Treatment
    )
  ) +
    
    geom_density(
      alpha = 0.35,
      linewidth = 0.65,
      na.rm = TRUE
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.5
    ) +
    
    facet_wrap(
      ~ Response,
      scales = "free",
      ncol = 2
    ) +
    
    scale_fill_manual(
      values = treatment_colours[
        c("Lime", "HB 50", "EG 50", "BD 50")
      ],
      breaks = c("Lime", "HB 50", "EG 50", "BD 50"),
      labels = treatment_display[
        c("Lime", "HB 50", "EG 50", "BD 50")
      ]
    ) +
    
    scale_colour_manual(
      values = treatment_colours[
        c("Lime", "HB 50", "EG 50", "BD 50")
      ],
      breaks = c("Lime", "HB 50", "EG 50", "BD 50"),
      labels = treatment_display[
        c("Lime", "HB 50", "EG 50", "BD 50")
      ]
    ) +
    
    labs(
      title = paste(
        "Posterior treatment differences:",
        display_group(group_name)
      ),
      subtitle = "Treatment minus Control; zero indicates no difference",
      x = "Posterior treatment difference (units specified in each panel heading)",
      y = "Density",
      fill = "Treatment",
      colour = "Treatment"
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey40"
      ),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  
  height <- max(
    5,
    ceiling(length(responses) / 2) * 3
  )
  
  
  filename <- file.path(
    bayesian_treatment_plots_dir,
    paste0(
      "Posterior_Contrast_Distributions_",
      group_name,
      ".png"
    )
  )
  
  
  ggsave(
    filename,
    p,
    width = 13,
    height = height,
    dpi = 300,
    bg = "white"
  )
  
  
  ggsave(
    sub("\\.png$", ".pdf", filename),
    p,
    width = 13,
    height = height,
    bg = "white"
  )
  
  
  print(p)
  
  p
  
}


contrast_density_plots <- list()


if (
  length(contrast_draw_file) > 0 &&
  !is.na(contrast_draw_file)
) {
  
  cat("\nLoading posterior contrast draws:\n")
  cat(contrast_draw_file, "\n")
  
  posterior_contrast_draws <- read_csv(
    contrast_draw_file,
    show_col_types = FALSE
  )
  
  
  required_contrast_columns <- c(
    "Response",
    "Treatment",
    "Difference"
  )
  
  
  missing_contrast_columns <- setdiff(
    required_contrast_columns,
    names(posterior_contrast_draws)
  )
  
  
  if (length(missing_contrast_columns) == 0) {
    
    for (group_name in names(plot_groups_available)) {
      
      contrast_density_plots[[group_name]] <-
        
        make_contrast_density_plot(
          data = posterior_contrast_draws,
          responses = plot_groups_available[[group_name]],
          group_name = group_name
        )
      
    }
    
  } else {
    
    warning(
      "Posterior contrast draws were found, but the expected columns ",
      "Response, Treatment and Difference are missing. ",
      "Contrast density plots were skipped."
    )
    
  }
  
} else {
  
  cat(
    "\nNo saved posterior contrast draws found. ",
    "Posterior contrast density plots skipped.\n"
  )
  
}


# =============================================================================
# 14. GENERATE PRIMARY AND PAIRWISE FOREST PLOTS
# =============================================================================

cat("\nGenerating primary forest plots...\n")

primary_forest_plots <- list()


for (group_name in names(plot_groups_available)) {
  
  primary_forest_plots[[group_name]] <-
    
    make_grouped_forest_plot(
      data = primary_contrasts,
      responses = plot_groups_available[[group_name]],
      group_name = group_name,
      comparison_type = "Primary"
    )
  
}


cat("\nGenerating all-pairwise forest plots...\n")

pairwise_forest_plots <- list()


for (group_name in names(plot_groups_available)) {
  
  pairwise_forest_plots[[group_name]] <-
    
    make_grouped_forest_plot(
      data = pairwise_contrasts,
      responses = plot_groups_available[[group_name]],
      group_name = group_name,
      comparison_type = "All_Pairwise"
    )
  
}


# =============================================================================
# 15. GENERATE RAW-DATA VIOLIN PLOTS
# =============================================================================

cat("\nGenerating raw-data violin plots...\n")

violin_plots <- list()


for (group_name in names(plot_groups_available)) {
  
  violin_plots[[group_name]] <-
    
    make_grouped_violin_plot(
      data = raw_data_plot,
      responses = plot_groups_available[[group_name]],
      group_name = group_name
    )
  
}


# =============================================================================
# 16. HÜHNERBERG DOSE-RESPONSE DATA
# =============================================================================

cat("\nPreparing Hühnerberg dose-response data...\n")


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
      
      TRUE ~ parse_number(Treatment)
      
    )
    
  ) %>%
  
  filter(!is.na(Dose_t_ha)) %>%
  
  pivot_longer(
    cols = any_of(plot_responses),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  
  filter(is.finite(Value))


# =============================================================================
# 17. HÜHNERBERG DOSE-RESPONSE PLOTS
# =============================================================================

make_grouped_dose_plot <- function(
    data,
    responses,
    group_name
) {
  
  plot_data <- data %>%
    
    filter(
      Response %in% responses,
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
  
  
  if (nrow(plot_data) == 0) return(NULL)
  
  
  dose_summary <- plot_data %>%
    
    group_by(Response, Dose_t_ha) %>%
    
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      n = sum(!is.na(Value)),
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
      title = paste(
        "Hühnerberg dose response:",
        display_group(group_name)
      ),
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
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  
  height <- max(
    5,
    ceiling(length(responses) / 2) * 3
  )
  
  
  filename <- file.path(
    bayesian_treatment_plots_dir,
    paste0(
      "Huhnerberg_Dose_Response_",
      group_name,
      ".png"
    )
  )
  
  
  ggsave(
    filename,
    p,
    width = 13,
    height = height,
    dpi = 300,
    bg = "white"
  )
  
  
  print(p)
  
  p
  
}


cat("\nGenerating Hühnerberg dose-response plots...\n")

dose_response_plots <- list()


for (group_name in names(plot_groups_available)) {
  
  dose_response_plots[[group_name]] <-
    
    make_grouped_dose_plot(
      data = hb_dose_data,
      responses = plot_groups_available[[group_name]],
      group_name = group_name
    )
  
}


# =============================================================================
# 18. COMPLETION
# =============================================================================

cat("\n\nGEOCHEMISTRY PLOT-ONLY RUN COMPLETE\n")

cat(
  "\nPlots saved in:\n",
  bayesian_treatment_plots_dir,
  "\n",
  sep = ""
)

cat("\nGenerated plot categories:\n")

cat("- Primary treatment forest plots\n")
cat("- All-pairwise treatment forest plots\n")
cat("- Raw-data treatment violin plots\n")
cat("- Hühnerberg dose-response plots\n")

if (length(posterior_treatment_density_plots) > 0) {
  cat("- Posterior treatment density distributions\n")
}

if (length(contrast_density_plots) > 0) {
  cat("- Posterior treatment-contrast density distributions\n")
}

cat("\nNo Bayesian models were fitted or sampled.\n")
