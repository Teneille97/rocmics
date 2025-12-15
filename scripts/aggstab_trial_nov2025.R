#load packages
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(car)
library(lme4)
library(ggthemr)
#set theme
ggthemr('fresh')
library(rstatix)
library(emmeans)
library(multcompView)

#import data
aggstab_data <- read.csv( here("csv_files", "aggstab_trial_nov2025.csv"), header = TRUE)

#data cleaning
aggstab_data[, 4:8] <- lapply(aggstab_data[, 4:8], as.numeric)

aggstab_data <- aggstab_data %>%
  rename(recovery = recovery....,
         coarse = X...2mm.fraction,
         macro = X...0.25mm.fraction,
         micro = X...0.063mm.fraction,
         siltclay = X...0.063mm.fraction.1)

#statistical test
agg_long_raw <- aggstab_data %>%
  select(extract, coarse, macro, micro, siltclay) %>%
  pivot_longer(
    cols = c(coarse, macro, micro, siltclay),
    names_to = "fraction",
    values_to = "value"
  ) %>%
  mutate(
    fraction = factor(
      fraction,
      levels = c("coarse", "macro", "micro", "siltclay")
    ),
    extract = factor(extract)
  )

normality_tests <- agg_long_raw %>%
  group_by(fraction) %>%
  do({
    model <- lm(value ~ extract, data = .)
    tibble(
      shapiro_p = shapiro.test(residuals(model))$p.value
    )
  })

normality_tests #all normal

variance_tests <- agg_long_raw %>%
  group_by(fraction) %>%
  do({
    test <- car::leveneTest(value ~ extract, data = .)
    tibble(
      levene_p = test$`Pr(>F)`[1]
    )
  })

variance_tests #all homoscedastic

#choose comparative test based on variance and normality of residual distribution
get_letters <- function(df) {
  
  model <- lm(value ~ extract, data = df)
  shapiro_p <- shapiro.test(residuals(model))$p.value
  levene_p  <- car::leveneTest(value ~ extract, data = df)$`Pr(>F)`[1]
  
  if (shapiro_p > 0.05 & levene_p > 0.05) {
    # Parametric
    aov_mod <- aov(value ~ extract, data = df)
    emm <- emmeans(aov_mod, ~ extract)
    cld <- multcomp::cld(
      emm,
      Letters = letters,
      adjust = "tukey"
    )
    out <- as.data.frame(cld)[, c("extract", ".group")]
  } else {
    # Non-parametric
    dunn <- dunn_test(df, value ~ extract, p.adjust.method = "bonferroni")
    letters <- multcompView::multcompLetters(
      setNames(dunn$p.adj, paste(dunn$group1, dunn$group2, sep = "-"))
    )$Letters
    out <- tibble(
      extract = names(letters),
      .group = letters
    )
  }
  
  out
}

#obtain cld letters from comparison tests
letters_df <- agg_long_raw %>%
  group_by(fraction) %>%
  group_modify(~ get_letters(.x)) %>%
  ungroup()

#add cld letters to df
agg_long_plot <- agg_long %>%
  left_join(letters_df, by = c("fraction", "extract"))


### - plot
aggstab_data_summary <- aggstab_data %>%
  group_by(extract) %>%
  summarise(
    mean_recovery = mean(recovery, na.rm = TRUE),
    se_recovery = sd(recovery, na.rm = TRUE) / sqrt(n()),
    mean_coarse = mean(coarse, na.rm = TRUE),
    se_coarse = sd(coarse, na.rm = TRUE) / sqrt(n()),
    mean_macro = mean(macro, na.rm = TRUE),
    se_macro = sd(macro, na.rm = TRUE) / sqrt(n()),
    mean_micro = mean(micro, na.rm = TRUE),
    se_micro = sd(micro, na.rm = TRUE) / sqrt(n()),
    mean_siltclay = mean(siltclay, na.rm = TRUE),
    se_siltclay = sd(siltclay, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

agg_long <- aggstab_data_summary %>%
  select(
    extract,
    mean_coarse, se_coarse,
    mean_macro,  se_macro,
    mean_micro,  se_micro,
    mean_siltclay, se_siltclay
  ) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "fraction",
    values_to = "mean"
  ) %>%
  mutate(
    fraction = str_remove(fraction, "mean_"),
    se = case_when(
      fraction == "coarse"   ~ se_coarse,
      fraction == "macro"    ~ se_macro,
      fraction == "micro"    ~ se_micro,
      fraction == "siltclay" ~ se_siltclay
    ),
    fraction = factor(
      fraction,
      levels = c("coarse", "macro", "micro", "siltclay")
    )
  )


ggplot(agg_long_plot, aes(x = fraction, y = mean, fill = fraction)) +
  geom_col(width = 0.7, colour = "grey20") +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = .group, y = mean + se + 2),
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(~ extract) +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(
    x = NULL,
    y = "Aggregate size class (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

#view statistical p-value results
get_pvalues <- function(df) {
  
  model <- lm(value ~ extract, data = df)
  shapiro_p <- shapiro.test(residuals(model))$p.value
  levene_p  <- car::leveneTest(value ~ extract, data = df)$`Pr(>F)`[1]
  
  if (shapiro_p > 0.05 & levene_p > 0.05) {
    
    # ANOVA
    aov_mod <- aov(value ~ extract, data = df)
    
    list(
      omnibus = tidy(aov_mod),
      pairwise = as.data.frame(TukeyHSD(aov_mod)$extract) %>%
        mutate(comparison = rownames(.))
    )
    
  } else {
    
    # Kruskal–Wallis
    list(
      omnibus = kruskal_test(df, value ~ extract),
      pairwise = dunn_test(
        df,
        value ~ extract,
        p.adjust.method = "bonferroni"
      )
    )
  }
}

pval_results <- agg_long_raw %>%
  group_by(fraction) %>%
  group_map(~ get_pvalues(.x))

#view results

pval_results[[which(levels(agg_long_raw$fraction) == "coarse")]]$omnibus
pval_results[[which(levels(agg_long_raw$fraction) == "coarse")]]$pairwise

pval_results[[which(levels(agg_long_raw$fraction) == "macro")]]$omnibus
pval_results[[which(levels(agg_long_raw$fraction) == "macro")]]$pairwise

pval_results[[which(levels(agg_long_raw$fraction) == "micro")]]$omnibus
pval_results[[which(levels(agg_long_raw$fraction) == "micro")]]$pairwise

pval_results[[which(levels(agg_long_raw$fraction) == "siltclay")]]$omnibus
pval_results[[which(levels(agg_long_raw$fraction) == "siltclay")]]$pairwise