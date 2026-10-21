# ==============================================================================
# LIBRARIES & SETUP
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
# IMPORT & CLEAN DATA
# ==============================================================================

zeta_potential <- read.csv(
  here("csv_files", "zeta_potential.csv")
)

zeta_potential$Treatment<-as.factor(zeta_potential$Treatment)

library(ggplot2)
library(dplyr)
library(multcompView)

# ============================================================
# 1. Check normality of model residuals
# ============================================================

model <- aov(Zeta ~ Treatment, data = zeta_potential)
shapiro_result <- shapiro.test(residuals(model))

print(shapiro_result)

# ============================================================
# 2. Select statistical test and calculate CLD letters
# ============================================================

if (shapiro_result$p.value >= 0.05) {
  
  # Normal residuals: ANOVA + Tukey HSD
  test_result <- anova(model)
  print(test_result)
  
  tukey_result <- TukeyHSD(model, "Treatment")
  print(tukey_result)
  
  p_values <- tukey_result$Treatment[, "p adj"]
  names(p_values) <- rownames(tukey_result$Treatment)
  
  cld <- multcompLetters(p_values, threshold = 0.05)$Letters
  
  test_used <- "ANOVA + Tukey HSD"
  
} else {
  
  # Non-normal residuals: Kruskal-Wallis + pairwise Wilcoxon
  test_result <- kruskal.test(Zeta ~ Treatment, data = zeta_potential)
  print(test_result)
  
  treatments <- levels(droplevels(zeta_potential$Treatment))
  pairs <- combn(treatments, 2, simplify = FALSE)
  
  p_values <- sapply(pairs, function(pair) {
    x <- zeta_potential$Zeta[
      zeta_potential$Treatment == pair[1]
    ]
    y <- zeta_potential$Zeta[
      zeta_potential$Treatment == pair[2]
    ]
    
    wilcox.test(x, y, exact = FALSE)$p.value
  })
  
  names(p_values) <- sapply(pairs, function(pair) {
    paste(pair[1], pair[2], sep = "-")
  })
  
  # Holm correction for multiple comparisons
  p_values <- p.adjust(p_values, method = "holm")
  print(p_values)
  
  cld <- multcompLetters(p_values, threshold = 0.05)$Letters
  
  test_used <- "Kruskal-Wallis + pairwise Wilcoxon (Holm)"
}

# ============================================================
# 3. Prepare letter positions
# ============================================================

letter_data <- data.frame(
  Treatment = names(cld),
  Letters = cld,
  stringsAsFactors = FALSE
)

letter_data$Treatment <- factor(
  letter_data$Treatment,
  levels = levels(zeta_potential$Treatment)
)

letter_data <- letter_data %>%
  mutate(
    y_position = max(zeta_potential$Zeta, na.rm = TRUE) +
      0.08 * diff(range(zeta_potential$Zeta, na.rm = TRUE))
  )

# ============================================================
# 4. Plot zeta potential with CLD letters
# ============================================================

p_zeta <- ggplot(
  zeta_potential,
  aes(x = Treatment, y = Zeta, fill = Treatment)
) +
  geom_boxplot(
    alpha = 0.5,
    outlier.shape = NA,
    width = 0.65
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.5,
    alpha = 0.8
  ) +
  geom_text(
    data = letter_data,
    aes(x = Treatment, y = y_position, label = Letters),
    inherit.aes = FALSE,
    colour = "black",
    size = 5
  ) +
  labs(
    x = "Treatment",
    y = "Zeta potential (mV)",
    title = "Zeta potential by treatment",
    subtitle = test_used
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )

print(p_zeta)
