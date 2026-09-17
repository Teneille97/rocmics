# ==============================================================================
# 1. LIBRARIES & SETUP
# ==============================================================================
library(tidyverse)
library(here)
library(emmeans)
library(multcomp)
library(multcompView)
library(ggplot2)
library(dplyr)
library(broom)

# Set global ggplot theme
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
# 3. IMPORT & CLEAN OXALATE DATA (mg g^-1)
# ==============================================================================
oxalate_raw <- read.csv(here("csv_files", "oxalate_sesq_feb2026.csv"))

oxalate_long <- cbind(oxalate_raw, metadata) %>%
  separate(Treatment, into = c("Tmt", "App_rate"), sep = "_", fill = "right") %>%
  mutate(
    App_rate = replace_na(App_rate, "0"),
    Tmt = factor(Tmt),
    App_rate = factor(App_rate, levels = c("0", "2", "4", "8", "12", "20", "30", "50")),
    Method = "Oxalate"
  ) %>%
  dplyr::select(sample, analysis_group, Tmt, App_rate, Method,
                Fe = Fe..mg.g.1.soil.,
                Al = Al..mg.g.1.soil.,
                Mn = Mn..mg.g.1.soil.,
                Si = Si..mg.g.1.soil.) %>%
  pivot_longer(
    cols = c(Fe, Al, Mn, Si),
    names_to = "Element",
    values_to = "Concentration"
  )

# ==============================================================================
# 4. IMPORT & CLEAN DCB DATA (Converted to mg g^-1: % * 10)
# ==============================================================================
dcb_raw <- read.csv(here("csv_files", "DCB_Feb2026.csv"))

dcb_long <- cbind(dcb_raw, metadata) %>%
  separate(Treatment, into = c("Tmt", "App_rate"), sep = "_", fill = "right") %>%
  mutate(
    App_rate = replace_na(App_rate, "0"),
    Tmt = factor(Tmt),
    App_rate = factor(App_rate, levels = c("0", "2", "4", "8", "12", "20", "30", "50")),
    Method = "DCB",
    Fe = if("X.Fe" %in% names(.)) X.Fe * 10 else NA_real_,
    Al = X.Al * 10,
    Mn = X.Mn * 10,
    Si = X.Si * 10
  ) %>%
  dplyr::select(sample, analysis_group, Tmt, App_rate, Method, Fe, Al, Mn, Si) %>%
  pivot_longer(
    cols = c(Fe, Al, Mn, Si),
    names_to = "Element",
    values_to = "Concentration"
  )

# ==============================================================================
# 5. UNIFIED DATASET & METRIC CALCULATIONS (DIFFERENCES & RATIOS)
# ==============================================================================
combined_extractions <- bind_rows(oxalate_long, dcb_long)

extraction_metrics <- combined_extractions %>%
  pivot_wider(
    names_from = Method,
    values_from = Concentration
  ) %>%
  mutate(
    Diff_DCB_Ox = DCB - Oxalate,
    Ratio_Ox_DCB = Oxalate / DCB
  )

metrics_long <- extraction_metrics %>%
  pivot_longer(
    cols = c(Diff_DCB_Ox, Ratio_Ox_DCB),
    names_to = "Metric",
    values_to = "Value"
  )

# ==============================================================================
# 6. COMBINED PLOTS: OXALATE VS DCB (SHARED SCALE)
# ==============================================================================
method_colors <- c("Oxalate" = "#E69F00", "DCB" = "#56B4E9")

# --- Plot 1: Full Dataset Overview by Application Rate & Method ---
ggplot(combined_extractions, aes(x = App_rate, y = Concentration, fill = Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.8)) +
  geom_point(aes(color = Method), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 1.2, alpha = 0.7) +
  facet_grid(Element ~ Tmt, scales = "free_y") +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  labs(
    x = expression(paste("Application Rate (t ha"^-1, ")")),
    y = expression(paste("Concentration (mg g"^-1, " soil)")),
    title = "Comparison of Oxalate and DCB Extractions Across Treatments"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# --- Plot 2: 50 t/ha Treatments Comparison (ANOVA & Tukey Letters) ---
plot50_combined <- combined_extractions %>%
  filter(analysis_group %in% c("Control", "Bolsdorfer_50", "Eifelgold_50", "Huhnerberg_50", "Lime_2")) %>%
  mutate(analysis_group = factor(analysis_group, levels = c("Control", "Bolsdorfer_50", "Eifelgold_50", "Huhnerberg_50", "Lime_2")))

letters_df_combined <- plot50_combined %>%
  group_by(Element, Method) %>%
  group_modify(~{
    mod <- aov(Concentration ~ analysis_group, data = .x)
    emm <- emmeans(mod, ~ analysis_group)
    cld_res <- cld(emm, Letters = letters, adjust = "tukey")
    
    cld_res %>%
      as.data.frame() %>%
      dplyr::select(analysis_group, .group) %>%
      mutate(
        .group = str_trim(.group),
        y = max(.x$Concentration, na.rm = TRUE) * 1.08
      )
  }) %>%
  ungroup()

ggplot(plot50_combined, aes(x = analysis_group, y = Concentration, fill = Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.8)) +
  geom_point(aes(color = Method), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 1.5) +
  geom_text(
    data = letters_df_combined,
    aes(x = analysis_group, y = y, label = .group, group = Method),
    position = position_dodge(0.8),
    size = 3.5,
    vjust = 0
  ) +
  facet_grid(Element ~ Method, scales = "free_y") +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  labs(
    x = "",
    y = expression(paste("Concentration (mg g"^-1, " soil)")),
    title = "50 t/ha Group Comparisons with Post-Hoc Tukey Groups"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# --- Plot 3: Dose-Response (Huhnerberg vs Control) ---
dose_resp_combined <- combined_extractions %>%
  filter(Tmt == "Control" | Tmt == "Huhnerberg") %>%
  mutate(
    Treatment = ifelse(Tmt == "Control", "Control", paste0(App_rate, " t ha^-1")),
    Treatment = factor(Treatment, levels = c("Control", "2 t ha^-1", "4 t ha^-1", "8 t ha^-1", "12 t ha^-1", "20 t ha^-1", "30 t ha^-1", "50 t ha^-1"))
  )

ggplot(dose_resp_combined, aes(x = Treatment, y = Concentration, fill = Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.8)) +
  geom_point(aes(color = Method), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1.8) +
  facet_grid(Element ~ ., scales = "free_y") +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  labs(
    x = "",
    y = expression(paste("Concentration (mg g"^-1, " soil)")),
    title = "Dose-Response Comparison (Huhnerberg)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# ==============================================================================
# 7. DIFFERENCE (DCB - Oxalate) AND RATIO (Oxalate / DCB) PLOTS
# ==============================================================================

# --- Plot 4: Absolute Difference (DCB - Oxalate) across 50 t/ha groups ---
metrics_50 <- extraction_metrics %>%
  filter(analysis_group %in% c("Control", "Bolsdorfer_50", "Eifelgold_50", "Huhnerberg_50", "Lime_2")) %>%
  mutate(analysis_group = factor(analysis_group, levels = c("Control", "Bolsdorfer_50", "Eifelgold_50", "Huhnerberg_50", "Lime_2")))

ggplot(metrics_50, aes(x = analysis_group, y = Diff_DCB_Ox, fill = analysis_group)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_grid(Element ~ ., scales = "free_y") +
  labs(
    x = "",
    y = expression(paste("Difference: DCB - Oxalate (mg g"^-1, " soil)")),
    title = "Crystalline/Residual Fraction (DCB - Oxalate)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# --- Plot 5: Extraction Ratio (Oxalate / DCB) across 50 t/ha groups ---
ggplot(metrics_50, aes(x = analysis_group, y = Ratio_Ox_DCB, fill = analysis_group)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_grid(Element ~ ., scales = "free_y") +
  labs(
    x = "",
    y = expression("Active Ratio (Oxalate / DCB)"),
    title = "Amorphous Fraction Ratio (Oxalate / DCB)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# ==============================================================================
# 8. TUKEY STATISTICAL ANALYSIS FOR AMORPHOUS FRACTION RATIO (Oxalate / DCB)
# ==============================================================================

# 1. Filter metrics data for the 50 t/ha comparison group
ratio_50 <- extraction_metrics %>%
  filter(analysis_group %in% c("Control", "Bolsdorfer_50", "Eifelgold_50", "Huhnerberg_50", "Lime_2")) %>%
  mutate(analysis_group = factor(analysis_group, levels = c("Control", "Bolsdorfer_50", "Eifelgold_50", "Huhnerberg_50", "Lime_2")))

# 2. Calculate Tukey HSD compact letter displays (CLD) per Element for Ratio_Ox_DCB
letters_ratio_df <- ratio_50 %>%
  group_by(Element) %>%
  group_modify(~{
    mod <- aov(Ratio_Ox_DCB ~ analysis_group, data = .x)
    emm <- emmeans(mod, ~ analysis_group)
    cld_results <- cld(
      emm,
      Letters = letters,
      adjust = "tukey"
    )
    
    cld_results %>%
      as.data.frame() %>%
      dplyr::select(analysis_group, .group) %>%
      mutate(
        .group = str_trim(.group),
        y = max(.x$Ratio_Ox_DCB, na.rm = TRUE) * 1.10
      )
  }) %>%
  ungroup()

# 3. Plot Amorphous Fraction Ratio with Tukey Letters
ggplot(ratio_50, aes(x = analysis_group, y = Ratio_Ox_DCB)) +
  geom_boxplot(outlier.shape = NA, fill = "grey92", color = "grey30") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8, color = "#2B5C8F") +
  geom_text(
    data = letters_ratio_df,
    aes(x = analysis_group, y = y, label = .group),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold"
  ) +
  facet_grid(Element ~ ., scales = "free_y") +
  labs(
    x = "",
    y = expression(paste("Active Ratio (Oxalate / DCB)")),
    title = "Amorphous Fraction Ratio (Oxalate / DCB) with Tukey HSD Groups",
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 13)
  )

# --- Check Normality ---
normality_results <- ratio_50 %>%
  filter(!is.na(Ratio_Ox_DCB)) %>%
  group_by(Element) %>%
  group_modify(~{
    mod <- aov(Ratio_Ox_DCB ~ analysis_group, data = .x)
    resids <- na.omit(residuals(mod))
    
    if (length(resids) >= 3 && var(resids) > 0) {
      st <- shapiro.test(resids)
      data.frame(
        W_statistic = unname(st$statistic),
        p_value     = st$p.value,
        Is_Normal   = ifelse(st$p.value > 0.05, "Yes (p > 0.05)", "No (p <= 0.05)")
      )
    } else {
      data.frame(
        W_statistic = NA_real_,
        p_value     = NA_real_,
        Is_Normal   = "Insufficient Data"
      )
    }
  }) %>%
  ungroup()

print(normality_results)

# ==============================================================================
# DOSE-RESPONSE RATIO PLOT (HÜHNERBERG)
# ==============================================================================
dose_resp_ratio <- extraction_metrics %>%
  filter(Tmt == "Control" | Tmt == "Huhnerberg") %>%
  mutate(
    Treatment = ifelse(Tmt == "Control", "Control", paste0(App_rate, " t ha^-1")),
    Treatment = factor(
      Treatment,
      levels = c("Control", "2 t ha^-1", "4 t ha^-1", "8 t ha^-1", "12 t ha^-1", "20 t ha^-1", "30 t ha^-1", "50 t ha^-1")
    )
  )

ggplot(dose_resp_ratio, aes(x = Treatment, y = Ratio_Ox_DCB, group = 1)) +
  geom_point(size = 3, alpha = 0.9) +
  facet_grid(Element ~ ., scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(
    x = "",
    y = expression(paste("Active Ratio (Oxalate / DCB)")),
    title = "Dose-Response Comparison (Hühnerberg)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 13)
  )

# ==============================================================================
# SILICATE VS NON-SILICATE & TYPE COMPARISONS ON OX/DCB RATIOS
# ==============================================================================

# Prepare factor levels & classes using extraction_metrics directly
stats_50_ratio <- extraction_metrics %>%
  filter(
    analysis_group %in% c("Control", "Bolsdorfer_50", "Eifelgold_50", "Huhnerberg_50", "Lime_2")
  ) %>%
  mutate(
    analysis_group = factor(
      analysis_group,
      levels = c("Control", "Bolsdorfer_50", "Eifelgold_50", "Huhnerberg_50", "Lime_2")
    ),
    amendment_class = case_when(
      analysis_group == "Control" ~ "Control",
      analysis_group == "Lime_2" ~ "Lime",
      analysis_group %in% c("Bolsdorfer_50", "Eifelgold_50", "Huhnerberg_50") ~ "Silicate"
    ),
    silicate_type = case_when(
      analysis_group == "Bolsdorfer_50" ~ "Bolsdorfer",
      analysis_group == "Eifelgold_50" ~ "Eifelgold",
      analysis_group == "Huhnerberg_50" ~ "Huhnerberg",
      TRUE ~ NA_character_
    )
  )

# 1. Overall silicate effect: Control vs Silicate (evaluated on Ratio_Ox_DCB per Element)
silicate_effect_results <- stats_50_ratio %>%
  filter(amendment_class %in% c("Control", "Silicate")) %>%
  group_by(Element) %>%
  group_modify(~ {
    mod <- aov(Ratio_Ox_DCB ~ amendment_class, data = .x)
    anova_tab <- summary(mod)[[1]]
    
    data.frame(
      silicate_F = anova_tab["amendment_class", "F value"],
      silicate_p = anova_tab["amendment_class", "Pr(>F)"]
    )
  }) %>%
  ungroup()

print(silicate_effect_results) #no sig diff

# 2. Silicate type effect: Bolsdorfer vs Eifelgold vs Huhnerberg (evaluated on Ratio_Ox_DCB per Element)
silicate_type_results <- stats_50_ratio %>%
  filter(amendment_class == "Silicate") %>%
  group_by(Element) %>%
  group_modify(~ {
    # ANOVA across silicate types for Ratio_Ox_DCB
    mod <- aov(Ratio_Ox_DCB ~ silicate_type, data = .x)
    
    # Estimated marginal means and Tukey comparisons
    emm <- emmeans(mod, ~ silicate_type)
    cld_res <- cld(emm, Letters = letters, adjust = "tukey")
    
    anova_tab <- summary(mod)[[1]]
    type_p <- anova_tab["silicate_type", "Pr(>F)"]
    
    cld_res %>%
      as.data.frame() %>%
      transmute(
        silicate_type = as.character(silicate_type),
        silicate_type_p = type_p,
        silicate_letter = str_trim(.group)
      )
  }) %>%
  ungroup()

print(silicate_type_results)

############### quick check CEC cor dcb ox
CEC_Feb_26 <- read.csv(here("csv_files", CEC_Feb_26.csv)) 
  
