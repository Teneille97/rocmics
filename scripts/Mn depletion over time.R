
# =============================================================================
# MALLE ERW: INITIAL Mn, THREE-YEAR Mn CHANGE AND fPOM RADIOCARBON
# PLOT-ONLY SCRIPT
#
# PURPOSE:
#   Investigate whether initial exchangeable Mn and its subsequent
#   redistribution/depletion are associated with final fPOM Δ14C.
#
# VARIABLES:
#   Initial Mn:       Mn_Oct2023
#   Three-year Mn:    CEC_Mn
#   Final radiocarbon:d14C_2026_fPOM
#
# IMPORTANT:
#   Verify that initial and final Mn measurements represent comparable
#   operational pools, extraction methods and units BEFORE calculating
#   their difference.
#
#   Mn_change = Mn_2026 - Mn_initial
#
#   Negative values = lower measured Mn at the final sampling.
#
#   No Bayesian models are fitted.
#
# OUTPUT:
#   outputs/Bayesian_treatment_comparisons/Mn_radiocarbon_analysis/
#
# =============================================================================


# =============================================================================
# 1. PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  
  library(tidyverse)
  library(here)
  library(ggtext)
  
})


# =============================================================================
# 2. PATHS AND SETTINGS
# =============================================================================

outputs_dir <- here::here("outputs")

input_file <- file.path(
  outputs_dir,
  "integrated_geochemistry_for_bayesian.csv"
)

output_dir <- file.path(
  outputs_dir,
  "Bayesian_treatment_comparisons",
  "Mn_radiocarbon_analysis"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 3. VARIABLE DEFINITIONS
# =============================================================================

# Confirm these assignments against your original dataset.

initial_mn_col <- "Mn_Oct2023"

final_mn_col <- "CEC_Mn"

radiocarbon_col <- "d14C_2026_fPOM"


# IMPORTANT:
#
# Set TRUE only if initial_mn_col and final_mn_col:
#   1. Represent the same operational Mn pool.
#   2. Were measured using comparable extraction methods.
#   3. Have the same units.
#
# If this is FALSE, the script will still generate the association plots
# but will not calculate Mn change or generate change-score plots.

calculate_mn_change <- TRUE


# Display labels.

initial_mn_label <- "Initial exchangeable Mn"

final_mn_label <- "Exchangeable Mn after three years"

radiocarbon_label <- "fPOM Δ¹⁴C after three years (‰)"

change_label <- expression(
  Delta*" exchangeable Mn (final - initial)"
)


# =============================================================================
# 4. LOAD DATA
# =============================================================================

if (!file.exists(input_file)) {
  
  stop("Integrated dataset not found: ", input_file)
  
}


integrated <- read_csv(
  input_file,
  na = c("", "NA", "#VALUE!"),
  show_col_types = FALSE
)


required_columns <- c(
  initial_mn_col,
  final_mn_col,
  radiocarbon_col,
  "Treatment"
)


missing_columns <- setdiff(
  required_columns,
  names(integrated)
)


if (length(missing_columns) > 0) {
  
  stop(
    "The following required columns are missing:\n",
    paste(missing_columns, collapse = "\n")
  )
  
}


cat("\nLoaded integrated dataset.\n")

cat("Initial Mn:", initial_mn_col, "\n")
cat("Final Mn:", final_mn_col, "\n")
cat("Radiocarbon:", radiocarbon_col, "\n")


# =============================================================================
# 5. PREPARE ANALYSIS DATA
# =============================================================================

mn_data <- integrated %>%
  
  mutate(
    
    Initial_Mn = as.numeric(.data[[initial_mn_col]]),
    
    Final_Mn = as.numeric(.data[[final_mn_col]]),
    
    fPOM_d14C = as.numeric(.data[[radiocarbon_col]]),
    
    Treatment_display = case_when(
      
      Treatment == "Control" ~ "Control",
      
      grepl("^Lime", Treatment) ~ "Lime",
      
      grepl("^Huhnerberg", Treatment) ~ "Hühnerberg",
      
      grepl("^Eifelgold", Treatment) ~ "Eifelgold",
      
      grepl("^Bolsdorfer", Treatment) ~ "Bolsdorfer",
      
      TRUE ~ as.character(Treatment)
      
    )
    
  ) %>%
  
  mutate(
    
    Treatment_display = factor(
      
      Treatment_display,
      
      levels = c(
        "Control",
        "Lime",
        "Hühnerberg",
        "Eifelgold",
        "Bolsdorfer"
      )
      
    )
    
  )


# Remove observations without any of the three measurements.

mn_data <- mn_data %>%
  
  filter(
    
    is.finite(Initial_Mn) |
      is.finite(Final_Mn) |
      is.finite(fPOM_d14C)
    
  )


# =============================================================================
# 6. CHECK DATA AND SAMPLE MATCHING
# =============================================================================

cat("\n============================================\n")
cat("DATA CHECKS\n")
cat("============================================\n")

cat("\nTotal observations:", nrow(mn_data), "\n")

cat(
  "Initial Mn observations:",
  sum(is.finite(mn_data$Initial_Mn)),
  "\n"
)

cat(
  "Final Mn observations:",
  sum(is.finite(mn_data$Final_Mn)),
  "\n"
)

cat(
  "fPOM radiocarbon observations:",
  sum(is.finite(mn_data$fPOM_d14C)),
  "\n"
)


if ("sample" %in% names(mn_data)) {
  
  cat(
    "\nNumber of unique sample identifiers:",
    n_distinct(mn_data$sample),
    "\n"
  )
  
}


if ("Sample" %in% names(mn_data)) {
  
  cat(
    "\nNumber of unique Sample identifiers:",
    n_distinct(mn_data$Sample),
    "\n"
  )
  
}


# =============================================================================
# 7. CALCULATE Mn CHANGE
# =============================================================================

if (calculate_mn_change) {
  
  mn_data <- mn_data %>%
    
    mutate(
      
      Mn_change = Final_Mn - Initial_Mn,
      
      Mn_depletion = Initial_Mn - Final_Mn
      
    )
  
  cat("\nMn change calculated:\n")
  
  cat("Mn_change = Final Mn - Initial Mn\n")
  
  cat("Negative Mn_change = lower final Mn\n")
  
  cat("Positive Mn_depletion = greater Mn depletion\n")
  
} else {
  
  mn_data$Mn_change <- NA_real_
  
  mn_data$Mn_depletion <- NA_real_
  
  cat(
    "\nMn change calculation disabled because measurement ",
    "comparability has not been confirmed.\n"
  )
  
}


# Save analysis data.

write_csv(
  
  mn_data,
  
  file.path(
    output_dir,
    "Initial_final_Mn_fPOM_radiocarbon_data.csv"
  )
  
)


# =============================================================================
# 8. TREATMENT COLOURS
# =============================================================================

treatment_colours <- c(
  
  "Control" = "#777777",
  
  "Lime" = "#78A641",
  
  "Hühnerberg" = "#C65A9E",
  
  "Eifelgold" = "#D99A16",
  
  "Bolsdorfer" = "#2878B5"
  
)


# =============================================================================
# 9. COMMON PLOT THEME
# =============================================================================

theme_mn <- function() {
  
  theme_bw(base_size = 12) +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      
      plot.subtitle = element_text(
        size = 10,
        colour = "grey30"
      ),
      
      axis.title = element_text(
        face = "bold",
        size = 11
      ),
      
      axis.text = element_text(
        colour = "black"
      ),
      
      legend.title = element_text(
        face = "bold"
      ),
      
      panel.grid.minor = element_blank(),
      
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.3
      ),
      
      legend.position = "bottom"
      
    )
  
}


# =============================================================================
# 10. CORRELATION FUNCTION
# =============================================================================

calculate_correlations <- function(
    data,
    x,
    y,
    comparison_name
) {
  
  complete_data <- data %>%
    
    filter(
      
      is.finite(.data[[x]]),
      
      is.finite(.data[[y]])
      
    )
  
  n <- nrow(complete_data)
  
  if (n < 3) {
    
    return(
      
      tibble(
        
        Comparison = comparison_name,
        
        n = n,
        
        Pearson_r = NA_real_,
        
        Pearson_p = NA_real_,
        
        Spearman_rho = NA_real_,
        
        Spearman_p = NA_real_
        
      )
      
    )
    
  }
  
  
  pearson <- cor.test(
    
    complete_data[[x]],
    
    complete_data[[y]],
    
    method = "pearson"
    
  )
  
  
  spearman <- suppressWarnings(
    
    cor.test(
      
      complete_data[[x]],
      
      complete_data[[y]],
      
      method = "spearman",
      
      exact = FALSE
      
    )
    
  )
  
  
  tibble(
    
    Comparison = comparison_name,
    
    n = n,
    
    Pearson_r = unname(pearson$estimate),
    
    Pearson_p = pearson$p.value,
    
    Spearman_rho = unname(spearman$estimate),
    
    Spearman_p = spearman$p.value
    
  )
  
}


# =============================================================================
# 11. PLOT SAVING FUNCTION
# =============================================================================

save_mn_plot <- function(
    plot,
    filename,
    width = 8,
    height = 6
) {
  
  ggsave(
    
    file.path(
      output_dir,
      paste0(filename, ".png")
    ),
    
    plot,
    
    width = width,
    
    height = height,
    
    dpi = 350,
    
    bg = "white"
    
  )
  
  
  ggsave(
    
    file.path(
      output_dir,
      paste0(filename, ".pdf")
    ),
    
    plot,
    
    width = width,
    
    height = height,
    
    bg = "white"
    
  )
  
  
  print(plot)
  
}


# =============================================================================
# 12. PLOT 1: INITIAL Mn VS FINAL Mn
# =============================================================================
#
# Question:
# Does initially exchangeable Mn correspond to residual exchangeable Mn
# after three years?
#
# Points below the 1:1 line indicate lower final Mn than initial Mn,
# provided both measurements have comparable units and extraction.
#


plot1_data <- mn_data %>%
  
  filter(
    
    is.finite(Initial_Mn),
    
    is.finite(Final_Mn)
    
  )


if (nrow(plot1_data) > 0) {
  
  plot1 <- ggplot(
    
    plot1_data,
    
    aes(
      
      x = Initial_Mn,
      
      y = Final_Mn,
      
      colour = Treatment_display
      
    )
    
  ) +
    
    geom_abline(
      
      slope = 1,
      
      intercept = 0,
      
      linetype = "dashed",
      
      colour = "black",
      
      linewidth = 0.7
      
    ) +
    
    geom_smooth(
      
      aes(group = 1),
      
      method = "lm",
      
      se = TRUE,
      
      colour = "black",
      
      fill = "grey75",
      
      linewidth = 0.8,
      
      inherit.aes = TRUE
      
    ) +
    
    geom_point(
      
      size = 3,
      
      alpha = 0.85
      
    ) +
    
    scale_colour_manual(
      
      values = treatment_colours,
      
      drop = FALSE
      
    ) +
    
    labs(
      
      title = "Initial versus residual exchangeable Mn",
      
      subtitle = paste(
        "Dashed line represents equal initial and final Mn",
        "| n =", nrow(plot1_data)
      ),
      
      x = paste0(
        initial_mn_label,
        " (",
        initial_mn_col,
        ")"
      ),
      
      y = paste0(
        final_mn_label,
        " (",
        final_mn_col,
        ")"
      ),
      
      colour = "Treatment"
      
    ) +
    
    theme_mn()
  
  
  save_mn_plot(
    
    plot1,
    
    "01_Initial_vs_Final_Mn",
    
    width = 8,
    
    height = 7
    
  )
  
}


# =============================================================================
# 13. PLOT 2: CHANGE IN Mn VS FINAL fPOM Δ14C
# =============================================================================
#
# Central hypothesis:
#
# Greater depletion of exchangeable Mn is associated with higher final
# fPOM Δ14C.
#
# Because Mn_change = final - initial:
#
#   Negative x values = Mn depletion.
#
#   Positive y direction = higher final fPOM Δ14C.
#
# Therefore, the hypothesised relationship is NEGATIVE:
# more negative Mn_change -> higher fPOM Δ14C.
#


if (calculate_mn_change) {
  
  plot2_data <- mn_data %>%
    
    filter(
      
      is.finite(Mn_change),
      
      is.finite(fPOM_d14C)
      
    )
  
  
  if (nrow(plot2_data) >= 3) {
    
    plot2 <- ggplot(
      
      plot2_data,
      
      aes(
        
        x = Mn_change,
        
        y = fPOM_d14C,
        
        colour = Treatment_display
        
      )
      
    ) +
      
      geom_hline(
        
        yintercept = 0,
        
        linetype = "dotted",
        
        colour = "grey50",
        
        linewidth = 0.5
        
      ) +
      
      geom_vline(
        
        xintercept = 0,
        
        linetype = "dashed",
        
        colour = "black",
        
        linewidth = 0.6
        
      ) +
      
      geom_smooth(
        
        aes(group = 1),
        
        method = "lm",
        
        se = TRUE,
        
        colour = "black",
        
        fill = "grey75",
        
        linewidth = 0.9
        
      ) +
      
      geom_point(
        
        size = 3.2,
        
        alpha = 0.85
        
      ) +
      
      scale_colour_manual(
        
        values = treatment_colours,
        
        drop = FALSE
        
      ) +
      
      labs(
        
        title = "Mn change versus final fPOM radiocarbon",
        
        subtitle = paste(
          "Hypothesis: greater Mn depletion is associated",
          "with higher fPOM Δ¹⁴C | n =", nrow(plot2_data)
        ),
        
        x = expression(
          Delta*" exchangeable Mn (final - initial)"
        ),
        
        y = radiocarbon_label,
        
        colour = "Treatment"
        
      ) +
      
      theme_mn()
    
    
    save_mn_plot(
      
      plot2,
      
      "02_Mn_Change_vs_fPOM_d14C",
      
      width = 8,
      
      height = 7
      
    )
    
  }
  
}


# =============================================================================
# 14. PLOT 3: INITIAL Mn VS FINAL fPOM Δ14C
# =============================================================================
#
# Question:
# Does initial exchangeable Mn remain positively associated with final
# fPOM Δ14C?
#
# This reproduces the initial Mn relationship independently of Mn change.
#


plot3_data <- mn_data %>%
  
  filter(
    
    is.finite(Initial_Mn),
    
    is.finite(fPOM_d14C)
    
  )


if (nrow(plot3_data) >= 3) {
  
  plot3 <- ggplot(
    
    plot3_data,
    
    aes(
      
      x = Initial_Mn,
      
      y = fPOM_d14C,
      
      colour = Treatment_display
      
    )
    
  ) +
    
    geom_hline(
      
      yintercept = 0,
      
      linetype = "dotted",
      
      colour = "grey50",
      
      linewidth = 0.5
      
    ) +
    
    geom_smooth(
      
      aes(group = 1),
      
      method = "lm",
      
      se = TRUE,
      
      colour = "black",
      
      fill = "grey75",
      
      linewidth = 0.9
      
    ) +
    
    geom_point(
      
      size = 3.2,
      
      alpha = 0.85
      
    ) +
    
    scale_colour_manual(
      
      values = treatment_colours,
      
      drop = FALSE
      
    ) +
    
    labs(
      
      title = "Initial exchangeable Mn versus final fPOM radiocarbon",
      
      subtitle = paste(
        "Relationship between initial Mn availability and",
        "fPOM Δ¹⁴C after three years | n =", nrow(plot3_data)
      ),
      
      x = initial_mn_label,
      
      y = radiocarbon_label,
      
      colour = "Treatment"
      
    ) +
    
    theme_mn()
  
  
  save_mn_plot(
    
    plot3,
    
    "03_Initial_Mn_vs_fPOM_d14C",
    
    width = 8,
    
    height = 7
    
  )
  
}


# =============================================================================
# 15. PLOT 4: INITIAL Mn VS Mn CHANGE
# =============================================================================
#
# Question:
# Do samples with initially greater Mn availability show greater subsequent
# Mn depletion?
#
# Points below zero indicate a decline in exchangeable Mn.
#
# Colour represents final fPOM Δ14C.
#
# This is particularly useful for examining the proposed combination:
#
# Initially high Mn + subsequent Mn depletion -> higher final fPOM Δ14C.
#


if (calculate_mn_change) {
  
  plot4_data <- mn_data %>%
    
    filter(
      
      is.finite(Initial_Mn),
      
      is.finite(Mn_change),
      
      is.finite(fPOM_d14C)
      
    )
  
  
  if (nrow(plot4_data) >= 3) {
    
    plot4 <- ggplot(
      
      plot4_data,
      
      aes(
        
        x = Initial_Mn,
        
        y = Mn_change,
        
        colour = fPOM_d14C
        
      )
      
    ) +
      
      geom_hline(
        
        yintercept = 0,
        
        linetype = "dashed",
        
        colour = "black",
        
        linewidth = 0.6
        
      ) +
      
      geom_smooth(
        
        aes(group = 1),
        
        method = "lm",
        
        se = TRUE,
        
        colour = "black",
        
        fill = "grey80",
        
        linewidth = 0.8
        
      ) +
      
      geom_point(
        
        size = 3.5,
        
        alpha = 0.9
        
      ) +
      
      scale_colour_gradient2(
        
        low = "#2878B5",
        
        mid = "grey90",
        
        high = "#C65A9E",
        
        midpoint = median(
          plot4_data$fPOM_d14C,
          na.rm = TRUE
        ),
        
        name = "fPOM Δ¹⁴C\n(‰)"
        
      ) +
      
      labs(
        
        title = "Initial Mn availability and subsequent Mn change",
        
        subtitle = paste(
          "Point colour indicates final fPOM Δ¹⁴C",
          "| n =", nrow(plot4_data)
        ),
        
        x = initial_mn_label,
        
        y = expression(
          Delta*" exchangeable Mn (final - initial)"
        )
        
      ) +
      
      theme_mn() +
      
      theme(
        
        legend.position = "right"
        
      )
    
    
    save_mn_plot(
      
      plot4,
      
      "04_Initial_Mn_vs_Mn_Change_Coloured_by_d14C",
      
      width = 9,
      
      height = 7
      
    )
    
  }
  
}


# =============================================================================
# 16. CORRELATION ANALYSIS
# =============================================================================

correlation_results <- bind_rows(
  
  calculate_correlations(
    
    mn_data,
    
    "Initial_Mn",
    
    "fPOM_d14C",
    
    "Initial Mn vs final fPOM Δ14C"
    
  ),
  
  calculate_correlations(
    
    mn_data,
    
    "Final_Mn",
    
    "fPOM_d14C",
    
    "Final Mn vs final fPOM Δ14C"
    
  )
  
)


if (calculate_mn_change) {
  
  correlation_results <- bind_rows(
    
    correlation_results,
    
    calculate_correlations(
      
      mn_data,
      
      "Mn_change",
      
      "fPOM_d14C",
      
      "Mn change vs final fPOM Δ14C"
      
    ),
    
    calculate_correlations(
      
      mn_data,
      
      "Initial_Mn",
      
      "Mn_change",
      
      "Initial Mn vs Mn change"
      
    ),
    
    calculate_correlations(
      
      mn_data,
      
      "Initial_Mn",
      
      "Final_Mn",
      
      "Initial Mn vs final Mn"
      
    )
    
  )
  
}


write_csv(
  
  correlation_results,
  
  file.path(
    
    output_dir,
    
    "Mn_radiocarbon_correlations.csv"
    
  )
  
)


cat("\n============================================\n")
cat("CORRELATION RESULTS\n")
cat("============================================\n")

print(
  
  correlation_results %>%
    
    mutate(
      
      across(
        
        where(is.numeric),
        
        ~ round(.x, 4)
        
      )
      
    )
  
)


# =============================================================================
# 17. OPTIONAL: TREATMENT-SPECIFIC CORRELATIONS
# =============================================================================

treatment_correlations <- mn_data %>%
  
  group_by(Treatment) %>%
  
  group_modify(
    
    ~ {
      
      bind_rows(
        
        calculate_correlations(
          
          .x,
          
          "Initial_Mn",
          
          "fPOM_d14C",
          
          "Initial Mn vs fPOM Δ14C"
          
        ),
        
        calculate_correlations(
          
          .x,
          
          "Final_Mn",
          
          "fPOM_d14C",
          
          "Final Mn vs fPOM Δ14C"
          
        ),
        
        if (calculate_mn_change) {
          
          calculate_correlations(
            
            .x,
            
            "Mn_change",
            
            "fPOM_d14C",
            
            "Mn change vs fPOM Δ14C"
            
          )
          
        }
        
      )
      
    }
    
  ) %>%
  
  ungroup()


write_csv(
  
  treatment_correlations,
  
  file.path(
    
    output_dir,
    
    "Mn_radiocarbon_correlations_by_treatment.csv"
    
  )
  
)


# =============================================================================
# 18. SUMMARY OF Mn CHANGE
# =============================================================================

if (calculate_mn_change) {
  
  mn_summary <- mn_data %>%
    
    group_by(Treatment) %>%
    
    summarise(
      
      n = sum(is.finite(Mn_change)),
      
      Initial_Mn_mean = mean(
        Initial_Mn,
        na.rm = TRUE
      ),
      
      Final_Mn_mean = mean(
        Final_Mn,
        na.rm = TRUE
      ),
      
      Mn_change_mean = mean(
        Mn_change,
        na.rm = TRUE
      ),
      
      Mn_change_SD = sd(
        Mn_change,
        na.rm = TRUE
      ),
      
      fPOM_d14C_mean = mean(
        fPOM_d14C,
        na.rm = TRUE
      ),
      
      .groups = "drop"
      
    )
  
  
  write_csv(
    
    mn_summary,
    
    file.path(
      
      output_dir,
      
      "Mn_change_summary_by_treatment.csv"
      
    )
    
  )
  
  
  cat("\n============================================\n")
  cat("Mn CHANGE BY TREATMENT\n")
  cat("============================================\n")
  
  print(mn_summary)
  
}


# =============================================================================
# 19. COMPLETION
# =============================================================================

cat("\n\nMN-RADIOCARBON ANALYSIS COMPLETE\n")

cat("\nPlots and results saved to:\n")

cat(output_dir, "\n")

cat("\nGenerated plots:\n")

cat("01. Initial versus final exchangeable Mn\n")

if (calculate_mn_change) {
  
  cat("02. Mn change versus final fPOM Δ14C\n")
  
}

cat("03. Initial Mn versus final fPOM Δ14C\n")

if (calculate_mn_change) {
  
  cat("04. Initial Mn versus Mn change, coloured by fPOM Δ14C\n")
  
}

cat("\nSaved correlation tables and analysis data.\n")

cat("\nNo Bayesian models were fitted.\n")

# =============================================================================
# END
# =============================================================================
