# =========================================================
# CEC AND EXCHANGEABLE CATIONS — FEBRUARY 2026
# ANOVA + TUKEY HSD + MEAN ± SE PLOTS
# =========================================================

# ---------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------

library(ggplot2)
library(gganimate)
library(ragg)
options(bitmapType = "cairo")
library(gifski)
library(viridis)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(car)
library(zoo)
library(ggthemr)
library(gridExtra)
library(here)
library(multcompView)

theme_update(
  legend.title = element_text(size = 10),
  legend.text  = element_text(size = 8)
)

# ---------------------------------------------------------
# 2. Theme
# ---------------------------------------------------------

ggthemr("flat dark")

# ---------------------------------------------------------
# 3. Import data
# ---------------------------------------------------------

CEC_Feb_26 <- read.csv(
  here("csv_files", "CEC_Feb_26.csv"),
  header = TRUE
)[1:252, ]

treatment_names <- read.csv(
  here("csv_files", "treatment_names.csv")
)

# ---------------------------------------------------------
# 4. Join treatment metadata
# ---------------------------------------------------------

CEC_Feb_26 <- CEC_Feb_26 %>%
  left_join(
    treatment_names %>%
      dplyr::rename(Sample = sample),
    by = "Sample"
  ) %>%
  na.omit()

# Make Treatment consistently character
CEC_Feb_26 <- CEC_Feb_26 %>%
  mutate(
    Treatment = as.character(Treatment)
  )

write.csv(CEC_Feb_26,file=here::here("outputs","CEC_Feb_26.csv"), row.names=FALSE)

# ---------------------------------------------------------
# 5. Define CEC properties
# ---------------------------------------------------------

cec_properties <- c(
  "H...cmolg.kg.1.",
  "CALCIUM..cmolg.kg.1.",
  "KALIUM..cmolg.kg.1.",
  "MAGNESIUM..cmolg.kg.1.",
  "NATRIUM..cmolg.kg.1.",
  "Aluminium..cmolg.kg.1.",
  "Ijzer..cmolg.kg.1.",
  "Mangaan..cmolg.kg.1.",
  "CEC..cmolg.kg.1.",
  "Basesaturation...."
)

# ---------------------------------------------------------
# 6. English facet labels
# ---------------------------------------------------------

property_labels <- c(
  "H...cmolg.kg.1."        = "H^'+'~(cmol[c]~kg^{-1})",
  "CALCIUM..cmolg.kg.1."   = "Calcium~(cmol[c]~kg^{-1})",
  "KALIUM..cmolg.kg.1."    = "Potassium~(cmol[c]~kg^{-1})",
  "MAGNESIUM..cmolg.kg.1." = "Magnesium~(cmol[c]~kg^{-1})",
  "NATRIUM..cmolg.kg.1."   = "Sodium~(cmol[c]~kg^{-1})",
  "Aluminium..cmolg.kg.1." = "Aluminium~(cmol[c]~kg^{-1})",
  "Ijzer..cmolg.kg.1."     = "Iron~(cmol[c]~kg^{-1})",
  "Mangaan..cmolg.kg.1."   = "Manganese~(cmol[c]~kg^{-1})",
  "CEC..cmolg.kg.1."       = "CEC~(cmol[c]~kg^{-1})",
  "Basesaturation...."      = "Base~saturation~plain('%')"
)

# ---------------------------------------------------------
# 7. ANOVA for every property
# ---------------------------------------------------------

anova_results <- list()
tukey_results <- list()

for (prop in cec_properties) {
  
  formula <- as.formula(
    paste0("`", prop, "` ~ Treatment")
  )
  
  model <- aov(
    formula,
    data = CEC_Feb_26
  )
  
  anova_results[[prop]] <- summary(model)
  
  # Only calculate Tukey when the overall ANOVA is significant
  p_value <- summary(model)[[1]][["Pr(>F)"]][1]
  
  if (!is.na(p_value) && p_value < 0.05) {
    tukey_results[[prop]] <- TukeyHSD(model)
  } else {
    tukey_results[[prop]] <- NULL
  }
}

# ---------------------------------------------------------
# 8. Print ANOVA p-values
# ---------------------------------------------------------

anova_pvalues <- data.frame(
  Property = cec_properties,
  p_value = sapply(
    anova_results,
    function(x) x[[1]][["Pr(>F)"]][1]
  )
)

print(anova_pvalues)

# ---------------------------------------------------------
# 9. Calculate Tukey letters
# ---------------------------------------------------------

tukey_letters <- list()

for (prop in cec_properties) {
  
  formula <- as.formula(
    paste0("`", prop, "` ~ Treatment")
  )
  
  model <- aov(
    formula,
    data = CEC_Feb_26
  )
  
  p_value <- summary(model)[[1]][["Pr(>F)"]][1]
  
  if (!is.na(p_value) && p_value < 0.05) {
    
    tukey <- TukeyHSD(model)
    
    letters <- multcompLetters4(
      model,
      tukey
    )
    
    tukey_letters[[prop]] <- data.frame(
      Treatment = names(letters$Treatment$Letters),
      Letters = as.character(letters$Treatment$Letters),
      Property = prop,
      stringsAsFactors = FALSE
    )
    
  } else {
    
    # No letters for non-significant ANOVAs
    tukey_letters[[prop]] <- data.frame(
      Treatment = unique(CEC_Feb_26$Treatment),
      Letters = "",
      Property = prop,
      stringsAsFactors = FALSE
    )
  }
}

tukey_letters_df <- bind_rows(tukey_letters) %>%
  mutate(
    Treatment = as.character(Treatment),
    Property = as.character(Property),
    Letters = as.character(Letters)
  )

# ---------------------------------------------------------
# 10. Calculate mean ± SE
# ---------------------------------------------------------

CEC_Feb_26_summary <- CEC_Feb_26 %>%
  group_by(Treatment) %>%
  summarise(
    across(
      all_of(cec_properties),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se = ~ sd(.x, na.rm = TRUE) /
          sqrt(sum(!is.na(.x)))
      ),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  )

# ---------------------------------------------------------
# 11. Convert summary data to long format
# ---------------------------------------------------------

CEC_long <- CEC_Feb_26_summary %>%
  pivot_longer(
    cols = -Treatment,
    names_to = c(".value", "Property"),
    names_pattern = "^(mean|se)_(.*)$"
  ) %>%
  mutate(
    Treatment = as.character(Treatment),
    Property = as.character(Property)
  )

# ---------------------------------------------------------
# 12. Keep only significant Tukey letters
# ---------------------------------------------------------

letter_positions <- tukey_letters_df %>%
  filter(
    !is.na(Letters),
    Letters != ""
  ) %>%
  left_join(
    CEC_long,
    by = c("Treatment", "Property")
  )

# ---------------------------------------------------------
# 13. Calculate positions for Tukey letters
# ---------------------------------------------------------

letter_positions <- letter_positions %>%
  group_by(Property) %>%
  mutate(
    y_max = max(
      mean + se,
      na.rm = TRUE
    ),
    
    y_min = min(
      mean - se,
      na.rm = TRUE
    ),
    
    y_range = y_max - y_min,
    
    # Prevent zero/invalid ranges
    y_range = if_else(
      !is.finite(y_range) | y_range == 0,
      1,
      y_range
    ),
    
    # Position letters slightly above each error bar
    y_position = mean + se + 0.05 * y_range
  ) %>%
  ungroup() %>%
  filter(
    !is.na(mean),
    !is.na(se),
    !is.na(y_position),
    is.finite(y_position)
  )

# ---------------------------------------------------------
# 14. Plot CEC and exchangeable cations
# ---------------------------------------------------------

CEC_plot <- ggplot(
  CEC_long,
  aes(
    x = Treatment,
    y = mean
  )
) +
  
  # Mean bars
  geom_col(
    width = 0.7
  ) +
  
  # Standard error
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.2,
    na.rm = TRUE
  ) +
  
  # Tukey letters
  geom_text(
    data = letter_positions,
    aes(
      x = Treatment,
      y = y_position,
      label = Letters
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4,
    vjust = 0
  ) +
  
  # Separate panel for each property
  facet_wrap(
    ~ Property,
    scales = "free_y",
    labeller = as_labeller(
      property_labels,
      default = label_parsed
    )
  ) +
  
  # Extra space above bars for Tukey letters
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.25)
    )
  ) +
  
  labs(
    x = "Treatment",
    y = "Mean ± SE"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

# ---------------------------------------------------------
# 15. Display plot
# ---------------------------------------------------------

print(CEC_plot)


# =========================================================
# CEC AND EXCHANGEABLE CATIONS — FEBRUARY 2026
# TWO PLOTS:
# 1. 50 t/ha treatments + Lime (ANOVA + Tukey HSD)
# 2. Huhnerberg dose-response (descriptive)
# Point symbols + mean ± SE
# =========================================================


# ---------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(dplyr)
library(here)
library(multcompView)
library(ggthemr)

options(bitmapType = "cairo")

ggthemr("flat dark")

theme_update(
  legend.title = element_text(size = 10),
  legend.text  = element_text(size = 8)
)


# ---------------------------------------------------------
# 2. Import data
# ---------------------------------------------------------

CEC_Feb_26 <- read.csv(
  here("csv_files", "CEC_Feb_26.csv"),
  header = TRUE
)[1:252, ]

treatment_names <- read.csv(
  here("csv_files", "treatment_names.csv")
)


# ---------------------------------------------------------
# 3. Join treatment metadata
# ---------------------------------------------------------

CEC_Feb_26 <- CEC_Feb_26 %>%
  left_join(
    treatment_names %>%
      dplyr::rename(Sample = sample),
    by = "Sample"
  ) %>%
  na.omit() %>%
  mutate(
    Treatment = as.character(Treatment)
  )


# Export joined data
write.csv(
  CEC_Feb_26,
  file = here("outputs", "CEC_Feb_26.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# 4. Define CEC properties
# ---------------------------------------------------------

cec_properties <- c(
  "H...cmolg.kg.1.",
  "CALCIUM..cmolg.kg.1.",
  "KALIUM..cmolg.kg.1.",
  "MAGNESIUM..cmolg.kg.1.",
  "NATRIUM..cmolg.kg.1.",
  "Aluminium..cmolg.kg.1.",
  "Ijzer..cmolg.kg.1.",
  "Mangaan..cmolg.kg.1.",
  "CEC..cmolg.kg.1.",
  "Basesaturation...."
)


# ---------------------------------------------------------
# 5. Facet labels
# ---------------------------------------------------------

property_labels <- c(
  "H...cmolg.kg.1."        = "H^'+'~(cmol[c]~kg^{-1})",
  "CALCIUM..cmolg.kg.1."   = "Calcium~(cmol[c]~kg^{-1})",
  "KALIUM..cmolg.kg.1."    = "Potassium~(cmol[c]~kg^{-1})",
  "MAGNESIUM..cmolg.kg.1." = "Magnesium~(cmol[c]~kg^{-1})",
  "NATRIUM..cmolg.kg.1."   = "Sodium~(cmol[c]~kg^{-1})",
  "Aluminium..cmolg.kg.1." = "Aluminium~(cmol[c]~kg^{-1})",
  "Ijzer..cmolg.kg.1."     = "Iron~(cmol[c]~kg^{-1})",
  "Mangaan..cmolg.kg.1."   = "Manganese~(cmol[c]~kg^{-1})",
  "CEC..cmolg.kg.1."       = "CEC~(cmol[c]~kg^{-1})",
  "Basesaturation...."     = "Base~saturation~plain('%')"
)


# ---------------------------------------------------------
# 6. Prepare treatment comparison data
# ---------------------------------------------------------

CEC_50 <- CEC_Feb_26 %>%
  filter(
    Treatment %in% c(
      "Control",
      "Bolsdorfer_50",
      "Eifelgold_50",
      "Huhnerberg_50",
      "Lime_2"
    )
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "Bolsdorfer_50",
        "Eifelgold_50",
        "Huhnerberg_50",
        "Lime_2"
      )
    )
  )


# Check treatment sample sizes
print(
  CEC_50 %>%
    group_by(Treatment) %>%
    summarise(
      n = n(),
      .groups = "drop"
    )
)


# =========================================================
# PLOT 1: 50 t/ha TREATMENTS + LIME
# =========================================================


# ---------------------------------------------------------
# 7. ANOVA + Tukey HSD
# ---------------------------------------------------------

anova_50_results <- list()
tukey_50_results <- list()
letters_50_results <- list()

for (prop in cec_properties) {
  
  formula <- as.formula(
    paste0("`", prop, "` ~ Treatment")
  )
  
  model <- aov(
    formula,
    data = CEC_50
  )
  
  anova_50_results[[prop]] <- summary(model)
  
  p_value <- summary(model)[[1]][["Pr(>F)"]][1]
  
  if (!is.na(p_value) && p_value < 0.05) {
    
    tukey <- TukeyHSD(model)
    
    tukey_50_results[[prop]] <- tukey
    
    letters <- multcompLetters4(
      model,
      tukey
    )
    
    letters_50_results[[prop]] <- data.frame(
      Treatment = names(letters$Treatment$Letters),
      Letters = as.character(letters$Treatment$Letters),
      Property = prop,
      stringsAsFactors = FALSE
    )
    
  } else {
    
    tukey_50_results[[prop]] <- NULL
    
    letters_50_results[[prop]] <- data.frame(
      Treatment = levels(CEC_50$Treatment),
      Letters = "",
      Property = prop,
      stringsAsFactors = FALSE
    )
  }
}


# ---------------------------------------------------------
# 8. Print ANOVA p-values
# ---------------------------------------------------------

anova_50_pvalues <- data.frame(
  Property = cec_properties,
  p_value = sapply(
    anova_50_results,
    function(x) x[[1]][["Pr(>F)"]][1]
  )
)

print(anova_50_pvalues)


# ---------------------------------------------------------
# 9. Calculate mean ± SE
# ---------------------------------------------------------

CEC_50_summary <- CEC_50 %>%
  group_by(Treatment) %>%
  summarise(
    across(
      all_of(cec_properties),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se = ~ sd(.x, na.rm = TRUE) /
          sqrt(sum(!is.na(.x)))
      ),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  )


CEC_50_long <- CEC_50_summary %>%
  pivot_longer(
    cols = -Treatment,
    names_to = c(".value", "Property"),
    names_pattern = "^(mean|se)_(.*)$"
  ) %>%
  mutate(
    Treatment = as.character(Treatment),
    Property = as.character(Property)
  )


# ---------------------------------------------------------
# 10. Prepare Tukey letters
# ---------------------------------------------------------

letters_50_df <- bind_rows(letters_50_results) %>%
  mutate(
    Treatment = as.character(Treatment),
    Property = as.character(Property),
    Letters = as.character(Letters)
  )


letter_positions_50 <- letters_50_df %>%
  filter(
    !is.na(Letters),
    Letters != ""
  ) %>%
  left_join(
    CEC_50_long,
    by = c("Treatment", "Property")
  ) %>%
  group_by(Property) %>%
  mutate(
    y_max = max(mean + se, na.rm = TRUE),
    y_min = min(mean - se, na.rm = TRUE),
    
    y_range = y_max - y_min,
    
    y_range = if_else(
      !is.finite(y_range) | y_range == 0,
      1,
      y_range
    ),
    
    y_position = mean + se + 0.08 * y_range
  ) %>%
  ungroup() %>%
  filter(
    !is.na(mean),
    !is.na(se),
    is.finite(y_position)
  )


# ---------------------------------------------------------
# 11. Plot 1 — 50 t/ha + Lime
#     Points + SE + Tukey letters
# ---------------------------------------------------------

CEC_plot_50 <- ggplot(
  CEC_50_long,
  aes(
    x = Treatment,
    y = mean
  )
) +
  
  geom_point(
    size = 3,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.2,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  
  geom_text(
    data = letter_positions_50,
    aes(
      x = Treatment,
      y = y_position,
      label = Letters
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4,
    vjust = 0
  ) +
  
  facet_wrap(
    ~ Property,
    scales = "free_y",
    labeller = as_labeller(
      property_labels,
      default = label_parsed
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.25)
    )
  ) +
  
  labs(
    x = "",
    y = "Mean ± SE",
    title = "CEC and Exchangeable Cations — 50 t ha-1 + Lime"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    plot.title = element_text(
      face = "bold",
      size = 13
    )
  )


print(CEC_plot_50)


# =========================================================
# PLOT 2: HUHNERBERG DOSE-RESPONSE
# =========================================================


# ---------------------------------------------------------
# 12. Filter Control + Huhnerberg dose-response
# ---------------------------------------------------------

CEC_dose <- CEC_Feb_26 %>%
  filter(
    Treatment == "Control" |
      Treatment %in% c(
        "Huhnerberg_2",
        "Huhnerberg_4",
        "Huhnerberg_8",
        "Huhnerberg_12",
        "Huhnerberg_20",
        "Huhnerberg_30",
        "Huhnerberg_50"
      )
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "Huhnerberg_2",
        "Huhnerberg_4",
        "Huhnerberg_8",
        "Huhnerberg_12",
        "Huhnerberg_20",
        "Huhnerberg_30",
        "Huhnerberg_50"
      ),
      labels = c(
        "Control",
        "2 t ha-1",
        "4 t ha-1",
        "8 t ha-1",
        "12 t ha-1",
        "20 t ha-1",
        "30 t ha-1",
        "50 t ha-1"
      )
    )
  )


# Check treatment sample sizes
print(
  CEC_dose %>%
    group_by(Treatment) %>%
    summarise(
      n = n(),
      .groups = "drop"
    )
)


# ---------------------------------------------------------
# 13. Calculate mean ± SE
# ---------------------------------------------------------

CEC_dose_summary <- CEC_dose %>%
  group_by(Treatment) %>%
  summarise(
    across(
      all_of(cec_properties),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se = ~ sd(.x, na.rm = TRUE) /
          sqrt(sum(!is.na(.x)))
      ),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  )


CEC_dose_long <- CEC_dose_summary %>%
  pivot_longer(
    cols = -Treatment,
    names_to = c(".value", "Property"),
    names_pattern = "^(mean|se)_(.*)$"
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "2 t ha-1",
        "4 t ha-1",
        "8 t ha-1",
        "12 t ha-1",
        "20 t ha-1",
        "30 t ha-1",
        "50 t ha-1"
      )
    ),
    Property = as.character(Property)
  )


# ---------------------------------------------------------
# 14. Plot 2 — Huhnerberg dose-response
#     Points + SE, no Tukey letters
# ---------------------------------------------------------

CEC_plot_dose <- ggplot(
  CEC_dose_long,
  aes(
    x = Treatment,
    y = mean,
    group = 1
  )
) +
  
  geom_point(
    size = 3,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.2,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ Property,
    scales = "free_y",
    labeller = as_labeller(
      property_labels,
      default = label_parsed
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  
  labs(
    x = "",
    y = "Mean ± SE",
    title = "CEC and Exchangeable Cations — Hühnerberg Dose-Response"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    plot.title = element_text(
      face = "bold",
      size = 13
    )
  )


print(CEC_plot_dose)
