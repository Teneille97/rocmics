library(tidyverse)
library(here)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(stringr)

radiocarbon_all_rawdata <- read.csv(here("csv_files", "14C_all_rawdata.csv"), header = TRUE)
df_proportions <- read.csv(here("outputs", "df_proportions.csv"), header = TRUE)

# =========================================================
# Prepare dataframe for radiocarbon modelling
# =========================================================

# ---------------------------------------------------------
# 1. Fraction data
# ---------------------------------------------------------

fraction_data <- df_proportions %>%
  transmute(
    sample = sample.ID,
    Plot = Plot,
    Sampling_year = Sampling_year,
    SOM_fraction = as.character(SOM_fraction),
    Tmt = Tmt,
    App_rate = App_rate,
    fraction_prop = fraction_prop,
    Ctotal = Ctotal,
    Ntotal = Ntotal,
    CN = CN,
    Cprop = Cprop,
    Nprop = Nprop
  )


# ---------------------------------------------------------
# 2. Create bulk dataframe
# ---------------------------------------------------------

bulk_data <- df_proportions %>%
  distinct(
    core_id,
    Plot,
    Sampling_year,
    Tmt,
    App_rate,
    C_bulk,
    N_bulk,
    CN_bulk
  ) %>%
  mutate(
    sample = paste0(core_id, "_Bulk"),
    SOM_fraction = "Bulk",
    fraction_prop = 100,
    Ctotal = C_bulk,
    Ntotal = N_bulk,
    CN = CN_bulk
  ) %>%
  transmute(
    sample = sample,
    Plot = Plot,
    Sampling_year = Sampling_year,
    SOM_fraction = SOM_fraction,
    Tmt = Tmt,
    App_rate = App_rate,
    fraction_prop = fraction_prop,
    Ctotal = Ctotal,
    Ntotal = Ntotal,
    CN = CN,
    
    # Bulk is 100% of bulk C and N
    Cprop = 100,
    Nprop = 100
  )

# ---------------------------------------------------------
# 3. Combine fractions + bulk
# ---------------------------------------------------------

radiocarbon_samples <- bind_rows(
  fraction_data,
  bulk_data
) %>%
  arrange(Sampling_year, Plot, SOM_fraction)


radiocarbon_samples <- radiocarbon_samples %>%
  left_join(
    radiocarbon_all_rawdata,
    by = c("sample" = "Sample")
  )


write.csv(
  radiocarbon_samples,
  here("outputs", "radiocarbon_samples.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Radiocarbon: main experiment
# ---------------------------------------------------------
main_14C <- radiocarbon_samples %>%
  filter(
    Tmt %in% c("Control", "Lime", "Eifelgold", "Bolsdorfer", "Huhnerberg")
  ) %>%
  # Filter to retain non-Huhnerberg treatments OR Huhnerberg only where Dose is 50
  filter(Tmt != "Huhnerberg" | App_rate == 50) %>%
  filter(!is.na(X.14C.....)) %>%
  mutate(
    SOM_fraction = factor(
      SOM_fraction,
      levels = c("Bulk", "fPOM", "oPOM", "MAOM")
    ),
    Tmt = factor(
      Tmt,
      levels = c("Control", "Lime", "Eifelgold", "Bolsdorfer", "Huhnerberg")
    )
  )

ggplot(
  main_14C,
  aes(
    x = SOM_fraction,
    y = X.14C.....
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65
  ) +
  geom_jitter(
    aes(colour = factor(Sampling_year)),  # Apply color mapping ONLY to jitter points
    width = 0.12,
    alpha = 0.7,
    size = 2
  ) +
  facet_wrap(~ Tmt) +
  theme_bw() +
  labs(
    x = "SOM fraction",
    y = expression(Delta^14*C~("\u2030")),
    title = expression(paste(Delta^14, "C of bulk soil and SOM fractions")),
    colour = "Sampling Year"
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3)))

# ---------------------------------------------------------
# Radiocarbon: Huhnerberg dose-response
# ---------------------------------------------------------

huhnerberg_14C <- radiocarbon_samples %>%
  filter(Tmt == "Huhnerberg") %>%
  filter(!is.na(X.14C.....)) %>%
  mutate(
    SOM_fraction = factor(
      SOM_fraction,
      levels = c("Bulk", "fPOM", "oPOM", "MAOM")
    )
  )

ggplot(
  huhnerberg_14C,
  aes(
    x = App_rate,
    y = X.14C.....
  )
) +
  geom_boxplot(
    aes(group = App_rate), # Groups numeric application rates for the boxplot
    outlier.shape = NA,
    width = 0.65
  ) +
  geom_jitter(
    aes(colour = factor(Sampling_year)), # Colored discrete points
    width = 0.12,
    alpha = 0.8,
    size = 3
  ) +
  theme_bw() +
  facet_wrap(~ SOM_fraction) +
  labs(
    x = "Huhnerberg application rate (t ha\u207B\u00B9)",
    y = expression(Delta^14*C~("\u2030")),
    colour = "Sampling Year",
    title = expression(paste(Delta^14, "C response to Huhnerberg application rate"))
  )

# =========================================================
# Radiocarbon statistics - MAIN EXPERIMENT
# =========================================================
#
# Experimental design
#
# 2023 = pre-treatment reference samples
#        Only a subset of representative plots/samples
#        were measured.
#        NOT plot-matched.
#
# 2025 = post-treatment, replicated experimental plots
# 2026 = post-treatment, same plots as 2025
#
# Main treatments:
#   Control
#   Lime
#   Eifelgold
#   Bolsdorfer
#   Huhnerberg at 50 t ha-1
#
# Huhnerberg dose-response treatments are excluded.
#
# Formal statistical analysis:
#   2025 + 2026 only
#
# Model:
#   Delta14C ~ Treatment * Year * SOM_fraction +
#             (1 | Plot)
#
# 2023 is summarized separately as a descriptive
# pre-treatment reference.
# =========================================================


# =========================================================
# 0. Packages
# =========================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(openxlsx)


# =========================================================
# 1. Prepare main experimental dataset
# =========================================================

main_14C <- radiocarbon_samples %>%
  
  # -------------------------------------------------------
# Keep only treatments in the main experiment
# -------------------------------------------------------

filter(
  Tmt %in% c(
    "Control",
    "Lime",
    "Eifelgold",
    "Bolsdorfer",
    "Huhnerberg"
  )
) %>%
  
  # -------------------------------------------------------
# Huhnerberg:
# retain ONLY the 50 t ha-1 treatment
# -------------------------------------------------------

filter(
  Tmt != "Huhnerberg" | App_rate == 50
) %>%
  
  # -------------------------------------------------------
# Formal analysis uses 2025 and 2026
# -------------------------------------------------------

filter(
  Sampling_year %in% c(2025, 2026)
) %>%
  
  # -------------------------------------------------------
# Remove missing radiocarbon measurements
# -------------------------------------------------------

filter(
  !is.na(X.14C.....)
) %>%
  
  mutate(
    
    # -----------------------------------------------------
    # Year
    # -----------------------------------------------------
    
    Sampling_year = factor(
      Sampling_year,
      levels = c(2025, 2026)
    ),
    
    # -----------------------------------------------------
    # SOM fraction
    # -----------------------------------------------------
    
    SOM_fraction = factor(
      SOM_fraction,
      levels = c(
        "Bulk",
        "fPOM",
        "oPOM",
        "MAOM"
      )
    ),
    
    # -----------------------------------------------------
    # Treatment
    # -----------------------------------------------------
    
    Tmt = factor(
      Tmt,
      levels = c(
        "Control",
        "Lime",
        "Eifelgold",
        "Bolsdorfer",
        "Huhnerberg"
      )
    ),
    
    # -----------------------------------------------------
    # Plot = repeated-measures unit
    # -----------------------------------------------------
    
    Plot = factor(Plot)
  )


# =========================================================
# 2. Check experimental structure
# =========================================================

cat("\n========================================\n")
cat("OBSERVATIONS BY YEAR / FRACTION / TREATMENT\n")
cat("========================================\n")

print(
  table(
    main_14C$Sampling_year,
    main_14C$SOM_fraction,
    main_14C$Tmt
  )
)


cat("\n========================================\n")
cat("NUMBER OF PLOTS BY YEAR / TREATMENT\n")
cat("========================================\n")

plot_counts <- main_14C %>%
  distinct(
    Plot,
    Sampling_year,
    Tmt
  ) %>%
  count(
    Sampling_year,
    Tmt
  )

print(plot_counts)


# =========================================================
# 3. Check that the same plots occur in both years
# =========================================================

plot_year_check <- main_14C %>%
  distinct(
    Plot,
    Sampling_year
  ) %>%
  count(Plot) %>%
  count(n)

cat("\n========================================\n")
cat("NUMBER OF YEARS AVAILABLE PER PLOT\n")
cat("========================================\n")

print(plot_year_check)


# Ideally, most/all plots should have n = 2.


# =========================================================
# 4. Main mixed-effects model
# =========================================================
#
# This is the primary statistical model.
#
# Fixed effects:
#   Tmt
#   Sampling_year
#   SOM_fraction
#
# Interactions:
#   Tmt × Year
#   Tmt × SOM fraction
#   Year × SOM fraction
#   Tmt × Year × SOM fraction
#
# Random effect:
#   Plot
#
# The random Plot effect accounts for the fact that the
# same experimental plots were measured in 2025 and 2026.
# =========================================================

model_14C <- lmer(
  X.14C..... ~
    Tmt *
    Sampling_year *
    SOM_fraction +
    (1 | Plot),
  data = main_14C,
  REML = FALSE
)


# =========================================================
# 5. Model summary
# =========================================================

cat("\n========================================\n")
cat("MODEL SUMMARY\n")
cat("========================================\n")

print(
  summary(model_14C)
)


# =========================================================
# 6. Type III ANOVA
# =========================================================

anova_14C <- anova(
  model_14C,
  type = 3
)

cat("\n========================================\n")
cat("TYPE III ANOVA\n")
cat("========================================\n")

print(anova_14C)


# Convert to dataframe for export
anova_df <- as.data.frame(anova_14C) %>%
  rownames_to_column("Term") %>%
  mutate(
    Significance = case_when(
      `Pr(>F)` < 0.001 ~ "***",
      `Pr(>F)` < 0.01  ~ "**",
      `Pr(>F)` < 0.05  ~ "*",
      `Pr(>F)` < 0.10  ~ ".",
      TRUE             ~ ""
    )
  )


# =========================================================
# 7. Treatment vs Control in 2025
# =========================================================
#
# Formal post-treatment treatment comparison.
#
# Comparisons are made separately within each SOM fraction.
# Holm adjustment controls for multiple treatment comparisons.
# =========================================================

emm_2025 <- emmeans(
  model_14C,
  ~ Tmt | SOM_fraction,
  at = list(
    Sampling_year = "2025"
  )
)


treatment_vs_control_2025 <- contrast(
  emm_2025,
  method = "trt.vs.ctrl",
  ref = 1,
  adjust = "holm"
)


treatment_vs_control_2025_df <- as.data.frame(
  summary(
    treatment_vs_control_2025,
    infer = c(TRUE, TRUE)
  )
)


cat("\n========================================\n")
cat("2025: TREATMENT VS CONTROL\n")
cat("========================================\n")

print(treatment_vs_control_2025_df)


# =========================================================
# 8. Treatment vs Control in 2026
# =========================================================

emm_2026 <- emmeans(
  model_14C,
  ~ Tmt | SOM_fraction,
  at = list(
    Sampling_year = "2026"
  )
)


treatment_vs_control_2026 <- contrast(
  emm_2026,
  method = "trt.vs.ctrl",
  ref = 1,
  adjust = "holm"
)


treatment_vs_control_2026_df <- as.data.frame(
  summary(
    treatment_vs_control_2026,
    infer = c(TRUE, TRUE)
  )
)


cat("\n========================================\n")
cat("2026: TREATMENT VS CONTROL\n")
cat("========================================\n")

print(treatment_vs_control_2026_df)


# =========================================================
# 9. Change from 2025 -> 2026 within each treatment
# =========================================================
#
# This tells us whether each treatment increased/decreased
# between the two post-treatment sampling years.
#
# IMPORTANT:
# This is not by itself the treatment effect.
# The key comparison is the next section:
# change vs Control.
# =========================================================

emm_year <- emmeans(
  model_14C,
  ~ Sampling_year | Tmt * SOM_fraction
)


year_change <- contrast(
  emm_year,
  method = "revpairwise",
  adjust = "holm"
)


year_change_df <- as.data.frame(
  summary(
    year_change,
    infer = c(TRUE, TRUE)
  )
)


cat("\n========================================\n")
cat("CHANGE FROM 2025 TO 2026\n")
cat("========================================\n")

print(year_change_df)


# =========================================================
# 10. Change vs Control
# =========================================================
#
# THIS IS THE KEY FOLLOW-UP TEST FOR THE YEAR INTERACTION.
#
# It asks:
#
# Did the change from 2025 -> 2026 differ between each
# treatment and the Control?
#
# Example:
#
# (Eifelgold 2026 - Eifelgold 2025)
# -
# (Control 2026 - Control 2025)
#
# A significant result means the trajectory over time
# differed between that treatment and Control.
# =========================================================

emm_tmt_year <- emmeans(
  model_14C,
  ~ Tmt * Sampling_year | SOM_fraction
)


change_vs_control <- contrast(
  emm_tmt_year,
  interaction = "revpairwise",
  by = "SOM_fraction",
  adjust = "holm"
)


change_vs_control_df <- as.data.frame(
  summary(
    change_vs_control,
    infer = c(TRUE, TRUE)
  )
)


cat("\n========================================\n")
cat("2025 -> 2026 CHANGE VS CONTROL\n")
cat("========================================\n")

print(change_vs_control_df)


# =========================================================
# 11. Estimated marginal means
# =========================================================
#
# Model-estimated means and 95% confidence intervals.
# =========================================================

emm_all <- emmeans(
  model_14C,
  ~ Tmt * Sampling_year | SOM_fraction
)


emm_all_df <- as.data.frame(
  summary(
    emm_all,
    infer = c(TRUE, TRUE)
  )
)


cat("\n========================================\n")
cat("ESTIMATED MARGINAL MEANS\n")
cat("========================================\n")

print(emm_all_df)


# =========================================================
# 12. Significant 2025 treatment effects
# =========================================================

significant_2025 <- treatment_vs_control_2025_df %>%
  filter(
    p.value < 0.05
  ) %>%
  arrange(
    SOM_fraction,
    p.value
  )


# =========================================================
# 13. Significant 2026 treatment effects
# =========================================================

significant_2026 <- treatment_vs_control_2026_df %>%
  filter(
    p.value < 0.05
  ) %>%
  arrange(
    SOM_fraction,
    p.value
  )


# =========================================================
# 14. Significant 2025 -> 2026 changes
# =========================================================

significant_year_change <- year_change_df %>%
  filter(
    p.value < 0.05
  ) %>%
  arrange(
    SOM_fraction,
    p.value
  )


# =========================================================
# 15. Significant change vs Control
# =========================================================

significant_change_vs_control <- change_vs_control_df %>%
  filter(
    p.value < 0.05
  ) %>%
  arrange(
    SOM_fraction,
    p.value
  )


# =========================================================
# 16. 2023 PRE-TREATMENT REFERENCE
# =========================================================
#
# IMPORTANT:
#
# These are NOT used in the mixed model.
#
# Because 2023 samples were only representative samples
# rather than measurements from every experimental plot,
# they cannot be used as a plot-level baseline.
#
# We calculate descriptive statistics only.
# =========================================================

baseline_2023 <- radiocarbon_samples %>%
  
  filter(
    Sampling_year == 2023
  ) %>%
  
  filter(
    !is.na(X.14C.....)
  ) %>%
  
  mutate(
    SOM_fraction = factor(
      SOM_fraction,
      levels = c(
        "Bulk",
        "fPOM",
        "oPOM",
        "MAOM"
      )
    )
  ) %>%
  
  group_by(
    SOM_fraction
  ) %>%
  
  summarise(
    
    n = n(),
    
    mean_D14C = mean(
      X.14C.....,
      na.rm = TRUE
    ),
    
    sd_D14C = sd(
      X.14C.....,
      na.rm = TRUE
    ),
    
    se_D14C = sd_D14C / sqrt(n),
    
    lower_95CI = mean_D14C -
      qt(
        0.975,
        df = n - 1
      ) * se_D14C,
    
    upper_95CI = mean_D14C +
      qt(
        0.975,
        df = n - 1
      ) * se_D14C,
    
    .groups = "drop"
  )


cat("\n========================================\n")
cat("2023 PRE-TREATMENT REFERENCE\n")
cat("========================================\n")

print(baseline_2023)


# =========================================================
# 17. Optional: 2023 reference by actual available
#     sample categories
# =========================================================
#
# This can be useful if the 2023 samples have multiple
# representative locations/samples and you want to see
# exactly how many observations contribute to each fraction.
# =========================================================

baseline_2023_raw <- radiocarbon_samples %>%
  filter(
    Sampling_year == 2023,
    !is.na(X.14C.....)
  ) %>%
  select(
    sample,
    Plot,
    Sampling_year,
    SOM_fraction,
    X.14C.....
  )


# =========================================================
# 18. Model diagnostics
# =========================================================

# Residual vs fitted
plot(
  model_14C,
  which = 1
)

# Normal Q-Q
qqnorm(
  residuals(model_14C)
)

qqline(
  residuals(model_14C)
)


# =========================================================
# 19. Create Excel workbook
# =========================================================

output_file <- here(
  "outputs",
  "radiocarbon_main_experiment_statistics_2025_2026.xlsx"
)


wb <- createWorkbook()


# ---------------------------------------------------------
# Sheet 1: Model ANOVA
# ---------------------------------------------------------

addWorksheet(
  wb,
  "Model_ANOVA"
)

writeData(
  wb,
  "Model_ANOVA",
  anova_df
)


# ---------------------------------------------------------
# Sheet 2: 2025 treatment vs Control
# ---------------------------------------------------------

addWorksheet(
  wb,
  "2025_Tmt_vs_Control"
)

writeData(
  wb,
  "2025_Tmt_vs_Control",
  treatment_vs_control_2025_df
)


# ---------------------------------------------------------
# Sheet 3: 2026 treatment vs Control
# ---------------------------------------------------------

addWorksheet(
  wb,
  "2026_Tmt_vs_Control"
)

writeData(
  wb,
  "2026_Tmt_vs_Control",
  treatment_vs_control_2026_df
)


# ---------------------------------------------------------
# Sheet 4: 2025 -> 2026 change
# ---------------------------------------------------------

addWorksheet(
  wb,
  "Year_Change"
)

writeData(
  wb,
  "Year_Change",
  year_change_df
)


# ---------------------------------------------------------
# Sheet 5: Change vs Control
# ---------------------------------------------------------

addWorksheet(
  wb,
  "Change_vs_Control"
)

writeData(
  wb,
  "Change_vs_Control",
  change_vs_control_df
)


# ---------------------------------------------------------
# Sheet 6: Estimated means
# ---------------------------------------------------------

addWorksheet(
  wb,
  "Estimated_Means"
)

writeData(
  wb,
  "Estimated_Means",
  emm_all_df
)


# ---------------------------------------------------------
# Sheet 7: Significant 2025
# ---------------------------------------------------------

addWorksheet(
  wb,
  "Significant_2025"
)

writeData(
  wb,
  "Significant_2025",
  significant_2025
)


# ---------------------------------------------------------
# Sheet 8: Significant 2026
# ---------------------------------------------------------

addWorksheet(
  wb,
  "Significant_2026"
)

writeData(
  wb,
  "Significant_2026",
  significant_2026
)


# ---------------------------------------------------------
# Sheet 9: Significant year changes
# ---------------------------------------------------------

addWorksheet(
  wb,
  "Significant_Year_Change"
)

writeData(
  wb,
  "Significant_Year_Change",
  significant_year_change
)


# ---------------------------------------------------------
# Sheet 10: Significant change vs Control
# ---------------------------------------------------------

addWorksheet(
  wb,
  "Significant_Change_vs_Control"
)

writeData(
  wb,
  "Significant_Change_vs_Control",
  significant_change_vs_control
)


# ---------------------------------------------------------
# Sheet 11: 2023 baseline reference
# ---------------------------------------------------------

addWorksheet(
  wb,
  "2023_Reference"
)

writeData(
  wb,
  "2023_Reference",
  baseline_2023
)


# ---------------------------------------------------------
# Sheet 12: Raw 2023 data
# ---------------------------------------------------------

addWorksheet(
  wb,
  "2023_Raw_Data"
)

writeData(
  wb,
  "2023_Raw_Data",
  baseline_2023_raw
)


# ---------------------------------------------------------
# Sheet 13: Data used in model
# ---------------------------------------------------------

addWorksheet(
  wb,
  "Data_Used"
)

writeData(
  wb,
  "Data_Used",
  main_14C
)


# =========================================================
# 20. Formatting
# =========================================================

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center"
)


for (sheet in names(wb)) {
  
  setColWidths(
    wb,
    sheet,
    cols = 1:20,
    widths = "auto"
  )
  
  freezePane(
    wb,
    sheet,
    firstRow = TRUE
  )
  
  addStyle(
    wb,
    sheet,
    style = header_style,
    rows = 1,
    cols = 1:20,
    gridExpand = TRUE
  )
  
}


# =========================================================
# 21. Save workbook
# =========================================================

saveWorkbook(
  wb,
  output_file,
  overwrite = TRUE
)


cat(
  "\n========================================\n",
  "ANALYSIS COMPLETE\n",
  "========================================\n",
  "Formal analysis: 2025 + 2026\n",
  "2023: descriptive pre-treatment reference only\n",
  "\nResults saved to:\n",
  output_file,
  "\n========================================\n"
)
