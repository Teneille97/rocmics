# =========================================================
# C PROPORTIONAL ALLOCATION TO SOM POOLS
# MAIN EXPERIMENT STATISTICAL ANALYSIS
# =========================================================
#
# df_proportions contains:
#
# Cprop = proportion of bulk C allocated to each SOM pool
#
# Main experiment:
#
# Control
# Lime
# Eifelgold
# Bolsdorfer
# Huhnerberg at 50 t ha-1
#
# Excluded:
# Huhnerberg dose-response treatments other than 50
#
# Years:
#
# 2023 = pre-treatment reference
# representative samples only; NOT plot-matched
#
# 2025 = post-treatment
# 2026 = post-treatment
# same experimental plots followed through time
#
# Formal analysis:
#
# Cprop ~ Tmt * Year * SOM_fraction + (1 | Plot)
#
# 2023 is NOT included in the inferential model.
# It is retained as a descriptive pre-treatment reference.
#
# =========================================================


# =========================================================
# Packages
# =========================================================

library(tidyverse)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(openxlsx)
library(here)


# =========================================================
# Import already-prepared dataframe
# =========================================================

df_proportions <- read.csv(
  here(
    "outputs",
    "df_proportions.csv"
  ),
  header = TRUE
)


# =========================================================
# Check dataframe
# =========================================================

str(df_proportions)


# =========================================================
# Repair missing 2025 bulk C and N data
# =========================================================
#
# The existing df_proportions already contains all of the
# fractionation data and 2026 Cprop values.
#
# However, the saved dataframe has missing C_bulk and N_bulk
# for 2025 because the 2025 bulk CN data had not been
# successfully joined before Cprop was calculated.
#
# We therefore:
#
# 1. Import the 2025 bulk CN measurements
# 2. Average technical/measurement replicates by sample
# 3. Match those samples to core_id
# 4. Populate C_bulk, N_bulk and CN_bulk for 2025
# 5. Recalculate Cprop and Nprop
#
# No fractionation data are recreated here.
#
# =========================================================


CN_data_nov2025 <- read.csv(
  here(
    "csv_files",
    "CN_nov2025.csv"
  ),
  header = TRUE
)


# Keep the 63 relevant 2025 bulk measurements
CN_data_nov2025 <- CN_data_nov2025[1:63, ]


# Convert C and N to numeric
CN_data_nov2025 <- CN_data_nov2025 %>%
  mutate(
    C = as.numeric(X.C),
    N = as.numeric(X.N),
    CN = C / N
  )


# Average measurements belonging to the same sample
CN_data_nov2025 <- CN_data_nov2025 %>%
  group_by(sample) %>%
  summarise(
    
    C_mean = mean(
      C,
      na.rm = TRUE
    ),
    
    C_se = sd(
      C,
      na.rm = TRUE
    ) / sqrt(sum(!is.na(C))),
    
    N_mean = mean(
      N,
      na.rm = TRUE
    ),
    
    N_se = sd(
      N,
      na.rm = TRUE
    ) / sqrt(sum(!is.na(N))),
    
    CN_mean = mean(
      CN,
      na.rm = TRUE
    ),
    
    CN_se = sd(
      CN,
      na.rm = TRUE
    ) / sqrt(sum(!is.na(CN))),
    
    .groups = "drop"
  )


# Create core IDs matching df_proportions
#
# Example:
#
# sample = 1
# becomes:
# 1_2025
#
CN_2025 <- CN_data_nov2025 %>%
  mutate(
    core_id = paste0(
      sample,
      "_2025"
    )
  ) %>%
  dplyr::select(
    core_id,
    C_mean,
    N_mean,
    CN_mean
  )


# Check the 2025 bulk data
cat("\n========================================\n")
cat("2025 BULK CN DATA CHECK\n")
cat("========================================\n")

cat(
  "Number of 2025 bulk samples:",
  nrow(CN_2025),
  "\n"
)

cat(
  "Missing mean C values:",
  sum(is.na(CN_2025$C_mean)),
  "\n"
)

cat(
  "Missing mean N values:",
  sum(is.na(CN_2025$N_mean)),
  "\n"
)

cat(
  "Missing mean CN values:",
  sum(is.na(CN_2025$CN_mean)),
  "\n"
)


# Check whether all 2025 fractionation core IDs
# have matching CN data

matching_2025_ids <- df_proportions %>%
  filter(
    Sampling_year == 2025
  ) %>%
  distinct(
    core_id
  ) %>%
  pull(
    core_id
  ) %>%
  `%in%`(
    CN_2025$core_id
  )

cat(
  "Number of matching 2025 core IDs:",
  sum(matching_2025_ids),
  "\n"
)


# Populate missing 2025 bulk C/N values

df_proportions <- df_proportions %>%
  left_join(
    CN_2025,
    by = "core_id"
  ) %>%
  mutate(
    
    C_bulk = if_else(
      Sampling_year == 2025,
      C_mean,
      C_bulk
    ),
    
    N_bulk = if_else(
      Sampling_year == 2025,
      N_mean,
      N_bulk
    ),
    
    CN_bulk = if_else(
      Sampling_year == 2025,
      CN_mean,
      CN_bulk
    )
    
  ) %>%
  dplyr::select(
    -C_mean,
    -N_mean,
    -CN_mean
  )


# Recalculate Cprop and Nprop
#
# fraction_prop = % of total soil mass represented by the SOM pool
#
# Ctotal = % C in that SOM pool
#
# C_bulk = % C in bulk soil
#
# Therefore:
#
# Cprop = percentage of total bulk C contained in the pool
#
df_proportions <- df_proportions %>%
  mutate(
    
    Cprop =
      (
        Ctotal *
          fraction_prop
      ) /
      C_bulk,
    
    Nprop =
      (
        Ntotal *
          fraction_prop
      ) /
      N_bulk
    
  )


# Check that the repair worked

cat("\n========================================\n")
cat("Cprop DATA CHECK\n")
cat("========================================\n")

Cprop_check <- df_proportions %>%
  group_by(
    Sampling_year
  ) %>%
  summarise(
    
    n = n(),
    
    C_bulk_missing =
      sum(is.na(C_bulk)),
    
    N_bulk_missing =
      sum(is.na(N_bulk)),
    
    Cprop_missing =
      sum(is.na(Cprop)),
    
    Nprop_missing =
      sum(is.na(Nprop)),
    
    .groups = "drop"
  )

print(
  Cprop_check
)


# Save repaired dataframe

write.csv(
  df_proportions,
  here(
    "outputs",
    "df_proportions.csv"
  ),
  row.names = FALSE
)


# =========================================================
# Prepare main experiment
# =========================================================
#
# IMPORTANT:
#
# 2023 is excluded from the formal model because the
# 2023 samples were representative samples rather than
# measurements from every experimental plot.
#
# 2025 and 2026 are the repeated-plot experiment.
#
# =========================================================

main_Cprop <- df_proportions %>%
  
  # Keep main treatments
  filter(
    Tmt %in% c(
      "Control",
      "Lime",
      "Eifelgold",
      "Bolsdorfer",
      "Huhnerberg"
    )
  ) %>%
  
  # Huhnerberg only at 50 t ha-1
  filter(
    Tmt != "Huhnerberg" |
      App_rate == 50
  ) %>%
  
  # Only post-treatment years
  filter(
    Sampling_year %in% c(
      2025,
      2026
    )
  ) %>%
  
  # Keep observations with Cprop
  filter(
    !is.na(Cprop)
  ) %>%
  
  # Convert variables to factors
  mutate(
    
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
    
    Sampling_year = factor(
      Sampling_year,
      levels = c(
        2025,
        2026
      )
    ),
    
    SOM_fraction = factor(
      SOM_fraction,
      levels = c(
        "fPOM",
        "oPOM",
        "MAOM"
      )
    ),
    
    Plot = factor(
      Plot
    )
    
  ) %>%
  
  # Remove unused factor levels after filtering
  droplevels()


# =========================================================
# Check experimental structure
# =========================================================

cat("\n========================================\n")
cat("MAIN Cprop DATASET\n")
cat("========================================\n")

cat(
  "Number of observations:",
  nrow(main_Cprop),
  "\n"
)

cat("\nTreatments:\n")

print(
  table(
    main_Cprop$Tmt
  )
)

cat("\nYears:\n")

print(
  table(
    main_Cprop$Sampling_year
  )
)

cat("\nSOM fractions:\n")

print(
  table(
    main_Cprop$SOM_fraction
  )
)

cat("\nTreatment × Year:\n")

print(
  table(
    main_Cprop$Tmt,
    main_Cprop$Sampling_year
  )
)

cat("\nTreatment × Year × SOM fraction:\n")

print(
  with(
    main_Cprop,
    table(
      Tmt,
      Sampling_year,
      SOM_fraction
    )
  )
)


# =========================================================
# Check that plots are repeated between years
# =========================================================

plot_years <- main_Cprop %>%
  
  distinct(
    Plot,
    Sampling_year
  ) %>%
  
  group_by(
    Plot
  ) %>%
  
  summarise(
    
    n_years =
      n_distinct(
        Sampling_year
      ),
    
    years =
      paste(
        sort(
          unique(
            Sampling_year
          )
        ),
        collapse = ", "
      ),
    
    .groups = "drop"
  )


cat("\n========================================\n")
cat("PLOT REPEAT CHECK\n")
cat("========================================\n")

print(
  plot_years
)


# =========================================================
# Main mixed-effects model
# =========================================================
#
# Cprop is the response.
#
# Fixed effects:
#
# Tmt
# Sampling_year
# SOM_fraction
#
# including all interactions.
#
# Plot is a random effect because the same experimental
# plots were followed from 2025 to 2026.
#
# =========================================================

model_Cprop <- lmer(
  
  Cprop ~
    Tmt *
    Sampling_year *
    SOM_fraction +
    (1 | Plot),
  
  data = main_Cprop,
  
  REML = FALSE
)


# =========================================================
# Model summary
# =========================================================

cat("\n========================================\n")
cat("Cprop MIXED MODEL\n")
cat("========================================\n")

print(
  summary(
    model_Cprop
  )
)


# =========================================================
# Type III ANOVA
# =========================================================

anova_Cprop <- anova(
  model_Cprop,
  type = 3
)

cat("\n========================================\n")
cat("TYPE III ANOVA\n")
cat("========================================\n")

print(
  anova_Cprop
)


anova_Cprop_df <-
  as.data.frame(
    anova_Cprop
  ) %>%
  
  rownames_to_column(
    "Term"
  ) %>%
  
  mutate(
    
    Significance =
      case_when(
        
        `Pr(>F)` < 0.001 ~ "***",
        `Pr(>F)` < 0.01  ~ "**",
        `Pr(>F)` < 0.05  ~ "*",
        `Pr(>F)` < 0.10  ~ ".",
        TRUE ~ ""
      )
  )


# =========================================================
# Estimated marginal means
# =========================================================

# Obtain estimated means for every combination of:
# Treatment × Year × SOM fraction

emm_Cprop <- emmeans(
  model_Cprop,
  ~ Tmt * Sampling_year * SOM_fraction
)

emm_Cprop_df <- as.data.frame(
  summary(
    emm_Cprop,
    infer = c(TRUE, TRUE)
  )
)

# =========================================================
# Treatment vs Control — separately for each year
# and SOM fraction
# =========================================================

# This directly asks:
#
# Within each SOM fraction and year,
# does each treatment differ from Control?
#
# 2025 and 2026 are analysed separately.

emm_treatment_year <- emmeans(
  model_Cprop,
  ~ Tmt | Sampling_year * SOM_fraction
)

Cprop_treatment_vs_control <- contrast(
  emm_treatment_year,
  method = "trt.vs.ctrl",
  ref = 1,
  adjust = "holm"
)

Cprop_treatment_vs_control_df <- as.data.frame(
  summary(
    Cprop_treatment_vs_control,
    infer = c(TRUE, TRUE)
  )
)

# Add significance column explicitly

Cprop_treatment_vs_control_df <- Cprop_treatment_vs_control_df %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.10  ~ ".",
      TRUE ~ ""
    )
  )

cat("\n========================================\n")
cat("TREATMENT VS CONTROL BY YEAR AND SOM FRACTION\n")
cat("========================================\n")

print(
  Cprop_treatment_vs_control_df
)

# =========================================================
# Separate 2025 results
# =========================================================

Cprop_treatment_vs_control_2025_df <-
  Cprop_treatment_vs_control_df %>%
  filter(
    Sampling_year == "2025"
  ) %>%
  arrange(
    SOM_fraction,
    p.value
  )

cat("\n========================================\n")
cat("2025 TREATMENT VS CONTROL\n")
cat("========================================\n")

print(
  Cprop_treatment_vs_control_2025_df
)

# =========================================================
# Separate 2026 results
# =========================================================

Cprop_treatment_vs_control_2026_df <-
  Cprop_treatment_vs_control_df %>%
  filter(
    Sampling_year == "2026"
  ) %>%
  arrange(
    SOM_fraction,
    p.value
  )

cat("\n========================================\n")
cat("2026 TREATMENT VS CONTROL\n")
cat("========================================\n")

print(
  Cprop_treatment_vs_control_2026_df
)

# =========================================================
# 2025 -> 2026 change within each treatment
# =========================================================

# Question:
#
# Did C allocation to each SOM pool change between
# 2025 and 2026 within each treatment?

emm_Cprop_year <- emmeans(
  model_Cprop,
  ~ Sampling_year | Tmt * SOM_fraction
)

Cprop_year_change <- contrast(
  emm_Cprop_year,
  method = "revpairwise",
  adjust = "holm"
)

Cprop_year_change_df <- as.data.frame(
  summary(
    Cprop_year_change,
    infer = c(TRUE, TRUE)
  )
)

Cprop_year_change_df <- Cprop_year_change_df %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.10  ~ ".",
      TRUE ~ ""
    )
  ) %>%
  arrange(
    SOM_fraction,
    Tmt,
    p.value
  )

cat("\n========================================\n")
cat("2025 -> 2026 CHANGE WITHIN TREATMENT\n")
cat("========================================\n")

print(
  Cprop_year_change_df
)

# =========================================================
# Change from 2025 -> 2026 versus Control
# =========================================================

# This is the key longitudinal treatment test.
#
# It tests whether the change between 2025 and 2026
# differs between each treatment and Control.
#
# For example:
#
# (Eifelgold_2026 - Eifelgold_2025)
#
# versus
#
# (Control_2026 - Control_2025)
#
# separately for each SOM fraction.

emm_Cprop_change <- emmeans(
  model_Cprop,
  ~ Tmt * Sampling_year | SOM_fraction
)

Cprop_change_vs_control <- contrast(
  emm_Cprop_change,
  interaction = "revpairwise",
  by = "SOM_fraction",
  adjust = "holm"
)

Cprop_change_vs_control_df <- as.data.frame(
  summary(
    Cprop_change_vs_control,
    infer = c(TRUE, TRUE)
  )
)

Cprop_change_vs_control_df <- Cprop_change_vs_control_df %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.10  ~ ".",
      TRUE ~ ""
    )
  ) %>%
  arrange(
    SOM_fraction,
    p.value
  )

cat("\n========================================\n")
cat("CHANGE 2025 -> 2026 VS CONTROL\n")
cat("========================================\n")

print(
  Cprop_change_vs_control_df
)

# =========================================================
# Significant 2025 treatment effects
# =========================================================

Cprop_significant_2025 <-
  Cprop_treatment_vs_control_2025_df %>%
  filter(
    p.value < 0.05
  ) %>%
  arrange(
    SOM_fraction,
    p.value
  )

# =========================================================
# Significant 2026 treatment effects
# =========================================================

Cprop_significant_2026 <-
  Cprop_treatment_vs_control_2026_df %>%
  filter(
    p.value < 0.05
  ) %>%
  arrange(
    SOM_fraction,
    p.value
  )

# =========================================================
# Significant year changes
# =========================================================

Cprop_significant_year_change <-
  Cprop_year_change_df %>%
  filter(
    p.value < 0.05
  ) %>%
  arrange(
    SOM_fraction,
    Tmt,
    p.value
  )

# =========================================================
# Significant change vs Control
# =========================================================

Cprop_significant_change_vs_control <-
  Cprop_change_vs_control_df %>%
  filter(
    p.value < 0.05
  ) %>%
  arrange(
    SOM_fraction,
    p.value
  )

# =========================================================
# Print significant results
# =========================================================

cat("\n========================================\n")
cat("SIGNIFICANT 2025 TREATMENT EFFECTS\n")
cat("========================================\n")

print(
  Cprop_significant_2025
)

cat("\n========================================\n")
cat("SIGNIFICANT 2026 TREATMENT EFFECTS\n")
cat("========================================\n")

print(
  Cprop_significant_2026
)

cat("\n========================================\n")
cat("SIGNIFICANT 2025 -> 2026 CHANGES\n")
cat("========================================\n")

print(
  Cprop_significant_year_change
)

cat("\n========================================\n")
cat("SIGNIFICANT CHANGE VS CONTROL\n")
cat("========================================\n")

print(
  Cprop_significant_change_vs_control
)

# =========================================================
# 2023 descriptive pre-treatment reference
# =========================================================
#
# NOT a formal treatment comparison.
#
# These were representative samples rather than
# measurements from every experimental plot.
#
# =========================================================

Cprop_2023_reference <-
  df_proportions %>%
  
  filter(
    Sampling_year == 2023,
    !is.na(Cprop)
  ) %>%
  
  group_by(
    SOM_fraction
  ) %>%
  
  summarise(
    
    n = n(),
    
    mean_Cprop =
      mean(
        Cprop,
        na.rm = TRUE
      ),
    
    sd_Cprop =
      sd(
        Cprop,
        na.rm = TRUE
      ),
    
    se_Cprop =
      sd_Cprop /
      sqrt(n),
    
    .groups = "drop"
  )


# =========================================================
# Print significant results
# =========================================================

cat("\n========================================\n")
cat("SIGNIFICANT 2025 TREATMENT EFFECTS\n")
cat("========================================\n")

print(
  Cprop_significant_2025
)

cat("\n========================================\n")
cat("SIGNIFICANT 2026 TREATMENT EFFECTS\n")
cat("========================================\n")

print(
  Cprop_significant_2026
)

cat("\n========================================\n")
cat("SIGNIFICANT 2025 -> 2026 CHANGES\n")
cat("========================================\n")

print(
  Cprop_significant_year_change
)

cat("\n========================================\n")
cat("SIGNIFICANT CHANGE VS CONTROL\n")
cat("========================================\n")

print(
  Cprop_significant_change_vs_control
)


# =========================================================
# Excel export
# =========================================================

output_file <- here(
  "outputs",
  "Cprop_main_experiment_statistics.xlsx"
)

wb <- createWorkbook()

# Model ANOVA

addWorksheet(wb, "Model_ANOVA")
writeData(
  wb,
  "Model_ANOVA",
  anova_Cprop_df
)

# Estimated marginal means

addWorksheet(wb, "Estimated_Means")
writeData(
  wb,
  "Estimated_Means",
  emm_Cprop_df
)

# All treatment vs control comparisons

addWorksheet(wb, "Tmt_vs_Control_All")
writeData(
  wb,
  "Tmt_vs_Control_All",
  Cprop_treatment_vs_control_df
)

# 2025 treatment vs control

addWorksheet(wb, "2025_Tmt_vs_Control")
writeData(
  wb,
  "2025_Tmt_vs_Control",
  Cprop_treatment_vs_control_2025_df
)

# 2026 treatment vs control

addWorksheet(wb, "2026_Tmt_vs_Control")
writeData(
  wb,
  "2026_Tmt_vs_Control",
  Cprop_treatment_vs_control_2026_df
)

# 2025 -> 2026 change

addWorksheet(wb, "Year_Change")
writeData(
  wb,
  "Year_Change",
  Cprop_year_change_df
)

# Change vs Control

addWorksheet(wb, "Change_vs_Control")
writeData(
  wb,
  "Change_vs_Control",
  Cprop_change_vs_control_df
)

# Significant 2025

addWorksheet(wb, "Significant_2025")
writeData(
  wb,
  "Significant_2025",
  Cprop_significant_2025
)

# Significant 2026

addWorksheet(wb, "Significant_2026")
writeData(
  wb,
  "Significant_2026",
  Cprop_significant_2026
)

# Significant year changes

addWorksheet(wb, "Significant_Year_Change")
writeData(
  wb,
  "Significant_Year_Change",
  Cprop_significant_year_change
)

# Significant change vs Control

addWorksheet(wb, "Significant_Change_vs_Control")
writeData(
  wb,
  "Significant_Change_vs_Control",
  Cprop_significant_change_vs_control
)

# 2023 descriptive reference

addWorksheet(wb, "2023_Reference")
writeData(
  wb,
  "2023_Reference",
  Cprop_2023_reference
)

# Data used in model

addWorksheet(wb, "Data_Used")
writeData(
  wb,
  "Data_Used",
  main_Cprop
)

# Format

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

# Save

saveWorkbook(
  wb,
  output_file,
  overwrite = TRUE
)

cat(
  "\n========================================\n",
  "Cprop STATISTICAL ANALYSIS COMPLETE\n",
  "========================================\n",
  "\nResults saved to:\n",
  output_file,
  "\n========================================\n"
)