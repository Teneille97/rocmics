# ==============================================================================
# TESSIER EXTRACTIONS — FEBRUARY 2026
# ABSOLUTE DIFFERENCE FROM CONTROL
# FOUR EXTRACTS: NH4OAc, ACETIC, HYDROXYLAMINE, H2O2
# ANOVA + TUKEY HSD + PLOTS
# ==============================================================================


# ==============================================================================
# 1. LIBRARIES & SETUP
# ==============================================================================

library(tidyverse)
library(here)
library(emmeans)
library(multcomp)
library(multcompView)
library(broom)
library(ggplot2)

theme_set(theme_bw())


# ==============================================================================
# 2. METADATA PROCESSING
# ==============================================================================

metadata <- read.csv(here("csv_files", "treatment_names.csv")) %>%
  mutate(
    analysis_group = case_when(
      Treatment == "Control" ~ "Control",
      Treatment == "Eifelgold_50" ~ "Eifelgold_50",
      Treatment == "Bolsdorfer_50" ~ "Bolsdorfer_50",
      Treatment == "Lime_2" ~ "Lime_2",
      Treatment == "Huhnerberg_50" ~ "Huhnerberg_50",
      TRUE ~ NA_character_
    )
  )


# ==============================================================================
# 3. IMPORT & CLEAN DATA
# ==============================================================================

tessier_feb26_nh4oac <- read.csv(
  here("csv_files", "tessier_feb2026_nh4oac.csv")
)

tessier_feb26_acetic <- read.csv(
  here("csv_files", "tessier_feb2026_acetic.csv")
)

tessier_feb26_hydroxylamine <- read.csv(
  here("csv_files", "tessier_feb2026_hydroxylamine.csv")
)

tessier_feb26_h2o2 <- read.csv(
  here("csv_files", "tessier_feb2026_h2o2.csv")
)


# ==============================================================================
# 4. PROCESS TESSIER EXTRACTS
# Convert concentrations to mg/g soil
# ==============================================================================

process_tessier_extract <- function(df, extract_name) {
  
  cbind(df, metadata) %>%
    separate(
      Treatment,
      into = c("Tmt", "App_rate"),
      sep = "_",
      fill = "right"
    ) %>%
    mutate(
      App_rate = replace_na(App_rate, "0"),
      
      Tmt = factor(Tmt),
      
      App_rate = factor(
        App_rate,
        levels = c("0", "2", "4", "8", "12", "20", "30", "50")
      ),
      
      Extract = extract_name,
      
      # Convert concentrations from mg/L to mg/g soil
      Fe = (Fe..mg.L. / (1000 / Corrected.extract.volume..mL.)) /
        Mass.soil..g.,
      
      Al = (Al..mg.L. / (1000 / Corrected.extract.volume..mL.)) /
        Mass.soil..g.,
      
      Mn = (Mn..mg.L. / (1000 / Corrected.extract.volume..mL.)) /
        Mass.soil..g.,
      
      Si = (Si..mg.L. / (1000 / Corrected.extract.volume..mL.)) /
        Mass.soil..g.,
      
      Ca = (Ca..mg.L. / (1000 / Corrected.extract.volume..mL.)) /
        Mass.soil..g.,
      
      K = (K..mg.L. / (1000 / Corrected.extract.volume..mL.)) /
        Mass.soil..g.,
      
      Na = (Na..mg.L. / (1000 / Corrected.extract.volume..mL.)) /
        Mass.soil..g.
    ) %>%
    dplyr::select(
      sample,
      analysis_group,
      Tmt,
      App_rate,
      Extract,
      Fe, Al, Mn, Si, Ca, K, Na
    ) %>%
    pivot_longer(
      cols = c(Fe, Al, Mn, Si, Ca, K, Na),
      names_to = "Element",
      values_to = "Concentration (mg g-1 soil)"
    )
}


# ==============================================================================
# 5. APPLY FUNCTION TO EACH EXTRACT
# ==============================================================================

tessier_feb26_nh4oac_long <- process_tessier_extract(
  tessier_feb26_nh4oac,
  "NH4OAc"
)

tessier_feb26_acetic_long <- process_tessier_extract(
  tessier_feb26_acetic,
  "Acetic"
)

tessier_feb26_hydroxylamine_long <- process_tessier_extract(
  tessier_feb26_hydroxylamine,
  "Hydroxylamine"
)

tessier_feb26_h2o2_long <- process_tessier_extract(
  tessier_feb26_h2o2,
  "H2O2"
)


# ==============================================================================
# 6. COMBINE ALL FOUR TESSIER EXTRACTS
# ==============================================================================

tessier_combined <- bind_rows(
  tessier_feb26_nh4oac_long,
  tessier_feb26_acetic_long,
  tessier_feb26_hydroxylamine_long,
  tessier_feb26_h2o2_long
) %>%
  mutate(
    Extract = factor(
      Extract,
      levels = c("NH4OAc", "Acetic", "Hydroxylamine", "H2O2")
    ),
    
    Element = factor(
      Element,
      levels = c("Fe", "Al", "Mn", "Si", "Ca", "K", "Na")
    )
  )


# ==============================================================================
# 7. CALCULATE ABSOLUTE DIFFERENCE FROM CONTROL
#
# Difference = individual concentration - Control mean
#
# Calculated separately for each Element x Extract combination
# Units remain mg/g soil
# ==============================================================================

control_means <- tessier_combined %>%
  filter(Tmt == "Control") %>%
  group_by(Element, Extract) %>%
  summarise(
    Control_mean = if (all(is.na(`Concentration (mg g-1 soil)`))) {
      NA_real_
    } else {
      mean(`Concentration (mg g-1 soil)`, na.rm = TRUE)
    },
    .groups = "drop"
  )


# Join Control means and calculate absolute differences

tessier_combined <- tessier_combined %>%
  left_join(
    control_means,
    by = c("Element", "Extract")
  ) %>%
  mutate(
    Difference_mg_g =
      `Concentration (mg g-1 soil)` - Control_mean
  )


# Check Control means and normalization

control_means

tessier_combined %>%
  filter(Tmt == "Control") %>%
  group_by(Element, Extract) %>%
  summarise(
    Mean_difference = mean(Difference_mg_g, na.rm = TRUE),
    .groups = "drop"
  )


# ==============================================================================
# 8. CONSISTENT EXTRACT COLOURS
# ==============================================================================

extract_colors <- c(
  "NH4OAc"        = "#E69F00",
  "Acetic"        = "#56B4E9",
  "Hydroxylamine" = "#009E73",
  "H2O2"          = "#CC79A7"
)


# ==============================================================================
# 9. FULL DATASET OVERVIEW — ALL APPLICATION RATES
# Absolute difference from Control
# ==============================================================================

ggplot(
  tessier_combined,
  aes(
    x = App_rate,
    y = Difference_mg_g,
    fill = Extract
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    position = position_dodge(0.8)
  ) +
  
  geom_point(
    aes(color = Extract),
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.8
    ),
    size = 1.2,
    alpha = 0.7
  ) +
  
  facet_grid(Element ~ Tmt, scales = "free_y") +
  
  scale_fill_manual(values = extract_colors) +
  
  scale_color_manual(values = extract_colors) +
  
  scale_x_discrete(
    labels = c(
      "0" = "0 (Control)",
      "2" = "2",
      "4" = "4",
      "8" = "8",
      "12" = "12",
      "20" = "20",
      "30" = "30",
      "50" = "50"
    )
  ) +
  
  labs(
    x = expression("Application rate (t ha"^-1*")"),
    y = expression("Difference from Control (mg g"^-1*" soil)"),
    title = "Tessier Element Concentrations — Difference from Control",
    subtitle = "Values are relative to the Control mean for each Element × Extract",
    fill = "Tessier Extract",
    color = "Tessier Extract"
  ) +
  
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )


# ==============================================================================
# 10. 50 t/ha TREATMENT COMPARISONS
# ANOVA + TUKEY HSD
# ==============================================================================

plot50_tessier <- tessier_combined %>%
  filter(
    analysis_group %in% c(
      "Control",
      "Bolsdorfer_50",
      "Eifelgold_50",
      "Huhnerberg_50",
      "Lime_2"
    )
  ) %>%
  mutate(
    analysis_group = factor(
      analysis_group,
      levels = c(
        "Control",
        "Bolsdorfer_50",
        "Eifelgold_50",
        "Huhnerberg_50",
        "Lime_2"
      )
    )
  )


# ==============================================================================
# 10a. ANOVA FOR EACH ELEMENT x EXTRACT
# ==============================================================================

anova_tessier_50 <- plot50_tessier %>%
  group_by(Element, Extract) %>%
  group_modify(~ {
    
    dat <- .x %>%
      filter(!is.na(Difference_mg_g)) %>%
      mutate(analysis_group = droplevels(analysis_group))
    
    n_groups <- n_distinct(dat$analysis_group)
    
    # Skip if fewer than 2 treatment groups
    if (n_groups < 2) {
      return(tibble(
        term = "analysis_group",
        df = NA_real_,
        statistic = NA_real_,
        p.value = NA_real_,
        n_groups = n_groups,
        status = "Skipped: fewer than 2 treatment groups"
      ))
    }
    
    # Fit ANOVA
    mod <- aov(
      Difference_mg_g ~ analysis_group,
      data = dat
    )
    
    # Extract results
    broom::tidy(mod) %>%
      mutate(
        n_groups = n_groups,
        status = "OK"
      )
  }) %>%
  ungroup()


# View all ANOVA results

anova_tessier_50


# View significant treatment effects

anova_tessier_50 %>%
  filter(
    term == "analysis_group",
    p.value < 0.05
  )


# ==============================================================================
# 10b. TUKEY HSD LETTERS FOR EACH ELEMENT x EXTRACT
# ==============================================================================

letters_tessier_50 <- plot50_tessier %>%
  group_by(Element, Extract) %>%
  group_modify(~ {
    
    dat <- .x %>%
      filter(!is.na(Difference_mg_g)) %>%
      mutate(analysis_group = droplevels(analysis_group))
    
    n_groups <- n_distinct(dat$analysis_group)
    
    # Skip if fewer than 2 groups or insufficient residual df
    if (
      n_groups < 2 ||
      nrow(dat) <= n_groups
    ) {
      return(tibble(
        analysis_group = factor(character()),
        .group = character(),
        y = numeric()
      ))
    }
    
    mod <- aov(
      Difference_mg_g ~ analysis_group,
      data = dat
    )
    
    emm <- emmeans(mod, ~ analysis_group)
    
    cld_res <- multcomp::cld(
      emm,
      Letters = letters,
      adjust = "tukey"
    )
    
    # Dynamic label height that also works with negative differences
    data_range <- range(dat$Difference_mg_g, na.rm = TRUE)
    padding <- diff(data_range) * 0.15
    
    if (padding == 0 || !is.finite(padding)) {
      padding <- 0.02
    }
    
    label_y <- max(dat$Difference_mg_g, na.rm = TRUE) + padding
    
    cld_res %>%
      as.data.frame() %>%
      dplyr::select(analysis_group, .group) %>%
      mutate(
        .group = stringr::str_trim(.group),
        y = label_y
      )
  }) %>%
  ungroup()


# View Tukey letters

letters_tessier_50


# ==============================================================================
# 10c. PLOT 50 t/ha COMPARISONS WITH TUKEY LETTERS
# ==============================================================================

plot50_tessier_labeled <- plot50_tessier %>%
  left_join(
    letters_tessier_50,
    by = c("Element", "Extract", "analysis_group")
  )


ggplot(
  plot50_tessier_labeled,
  aes(
    x = analysis_group,
    y = Difference_mg_g,
    fill = Extract
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    width = 0.65
  ) +
  
  geom_point(
    aes(color = Extract),
    position = position_jitter(width = 0.08),
    size = 1.5,
    alpha = 0.8
  ) +
  
  geom_text(
    data = letters_tessier_50,
    aes(
      x = analysis_group,
      y = y,
      label = .group
    ),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = 0
  ) +
  
  facet_grid(Element ~ Extract, scales = "free_y") +
  
  scale_fill_manual(values = extract_colors) +
  
  scale_color_manual(values = extract_colors) +
  
  scale_x_discrete(
    labels = c(
      "Control" = "Control",
      "Bolsdorfer_50" = "Bolsdorfer\n50 t/ha",
      "Eifelgold_50" = "Eifelgold\n50 t/ha",
      "Huhnerberg_50" = "Hühnerberg\n50 t/ha",
      "Lime_2" = "Lime\n2 t/ha"
    )
  ) +
  # Add extra space above each facet for Tukey letters
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.30))
  ) +
  
  # Prevent labels from being clipped at panel boundaries
  coord_cartesian(clip = "off") +
  labs(
    x = "Treatment",
    y = expression("Difference from Control (mg g"^-1*" soil)"),
    title = "Tessier Extract Comparisons — 50 t/ha Treatments",
    subtitle = "Absolute difference from Control; Tukey HSD within each Element × Extract",
    fill = "Extract",
    color = "Extract"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.margin = margin(
      t = 15,
      r = 15,
      b = 15,
      l = 15
    )
  )


# ==============================================================================
# 11. HÜHNERBERG DOSE-RESPONSE
# CONTROL VS BASALT APPLICATION RATES
# NO TUKEY: UNREPLICATED DOSE-RESPONSE TREATMENTS
# ==============================================================================

dose_resp_tessier <- tessier_combined %>%
  filter(Tmt %in% c("Control", "Huhnerberg")) %>%
  mutate(
    Treatment = ifelse(
      Tmt == "Control",
      "Control (0 t/ha)",
      paste0(as.character(App_rate), " t/ha")
    ),
    
    Treatment = factor(
      Treatment,
      levels = c(
        "Control (0 t/ha)",
        "2 t/ha",
        "4 t/ha",
        "8 t/ha",
        "12 t/ha",
        "20 t/ha",
        "30 t/ha",
        "50 t/ha"
      )
    )
  )


# ==============================================================================
# 11a. DOSE-RESPONSE PLOT
# Points and lines; no Tukey testing
# ==============================================================================

ggplot(
  dose_resp_tessier,
  aes(
    x = Treatment,
    y = Difference_mg_g,
    color = Extract,
    group = Extract
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  
  geom_point(
    position = position_dodge(width = 0.5),
    size = 2,
    alpha = 0.85
  ) +
  
  geom_line(
    position = position_dodge(width = 0.5),
    linewidth = 0.6,
    alpha = 0.7
  ) +
  
  facet_grid(Element ~ ., scales = "free_y") +
  
  scale_color_manual(values = extract_colors) +
  
  labs(
    x = "Hühnerberg application rate (t/ha)",
    y = expression("Difference from Control (mg g"^-1*" soil)"),
    title = "Tessier Element Concentrations — Hühnerberg Dose-Response",
    subtitle = "Absolute difference from Control; no Tukey testing for unreplicated dose rates",
    color = "Tessier Extract"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

