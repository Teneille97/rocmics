#load packages
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(car)
library(lme4)
library(lmerTest)
library(ggthemr)
#set theme
ggthemr('fresh')
library(rstatix)
library(emmeans)
library(multcomp)
library(multcompView)
library(here)
library(ARTool)
library(tidyr)

# Read treatment names once
treatment_names <- read.csv(here("csv_files", "treatment_names.csv"), header = TRUE)
tmt_names_vec <- rep(treatment_names$Treatment, each = 2)

# -------------------------------------------------
# FUNCTION
# -------------------------------------------------

clean_aggstab_extract <- function(filename, extract_label) {
  
  # Import data
  df <- read.csv(here("csv_files", filename), header = TRUE)
  
  # Add treatment names
  df$treatment <- tmt_names_vec
  
  # Convert to factors
  df[c("extract", "rep", "treatment")] <- 
    lapply(df[c("extract", "rep", "treatment")], as.factor)
  
  # Select columns
  df <- df %>%
    dplyr::select(sample, extract, rep, recovery....,
                  X...2mm.fraction,
                  X...0.25mm.fraction,
                  X...0.063mm.fraction,
                  X...0.063mm.fraction.1,
                  treatment)
  
  # Rename columns
  df <- df %>%
    rename(
      recovery = recovery....,
      `>2mm.fraction`   = X...2mm.fraction,
      `>0.25mm.fraction` = X...0.25mm.fraction,
      `>0.063mm.fraction` = X...0.063mm.fraction,
      `<0.063mm.fraction` = X...0.063mm.fraction.1
    )
  
  # Calculate MWD (mean weight diameter) **before pivoting**
  df <- df %>%
    mutate(
      MWD = `>2mm.fraction` * 3 +
        `>0.25mm.fraction` * 1.125 +
        `>0.063mm.fraction` * 0.1565 +
        `<0.063mm.fraction` * 0.0315
    )
  
  # Pivot to long format
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = c(`>2mm.fraction`,
               `>0.25mm.fraction`,
               `>0.063mm.fraction`,
               `<0.063mm.fraction`),
      names_to = "fraction_type",
      values_to = "fraction"
    ) %>%
    mutate(
      fraction_type = case_when(
        fraction_type == ">2mm.fraction"   ~ ">2",
        fraction_type == ">0.25mm.fraction" ~ ">0.25",
        fraction_type == ">0.063mm.fraction" ~ ">0.063",
        fraction_type == "<0.063mm.fraction" ~ "<0.063"
      ),
      fraction_type = factor(
        fraction_type,
        levels = c(">2", ">0.25", ">0.063", "<0.063")
      ),
      extract_type = factor(
        extract_label,
        levels = c("h2o", "oxalate", "h2o2", "naoh")
      )
    )
  
  return(as.data.frame(df_long))
}


#apply to my dataframes
aggstab_h2o      <- clean_aggstab_extract("aggstab_2026_h2o.csv", "h2o")
aggstab_oxalate  <- clean_aggstab_extract("aggstab_2026_oxalate.csv", "oxalate")
aggstab_h2o2     <- clean_aggstab_extract("aggstab_2026_h2o2.csv", "h2o2")
aggstab_naoh     <- clean_aggstab_extract("aggstab_2026_naoh.csv", "naoh")

#merge into one dataframe
aggstab_all <- bind_rows(
  aggstab_h2o,
  aggstab_oxalate,
  aggstab_h2o2,
  aggstab_naoh
)

#aggstab_all <- aggstab_all %>%
  #select(-extract) #remove old column; temporarily not working after restart

#check final df
str(aggstab_all)

#statistical comparison tests
run_fraction_logit_lmm <- function(data) {
  
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(dplyr)
  
  # Ensure factors
  data <- data %>%
    mutate(
      fraction_type = factor(
        fraction_type,
        levels = c(">2", ">0.25", ">0.063", "<0.063")
      ),
      treatment = factor(treatment),
      extract_type = factor(extract_type),
      sample = factor(sample)
    )
  
  # Convert % to proportion
  data <- data %>%
    mutate(
      # Ensure values are strictly between 0 and 1
      fraction_prop = (fraction / 100),
      fraction_prop = pmax(0.0001, pmin(0.9999, fraction_prop)), 
      fraction_logit = qlogis(fraction_prop)
    )
  
  # Fit mixed model
  lmm <- lmer(fraction_logit ~ extract_type * fraction_type * treatment + (1 | sample), data = data)
  
  cat("\nSingular fit?: ", isSingular(lmm), "\n")
  
  # Type III ANOVA (more appropriate for factorial design)
  anova_res <- anova(lmm, type = 3)
  
  # Estimated marginal means (back-transformed)
  emm <- emmeans(
    lmm,
    ~ extract_type | fraction_type * treatment
  )
  
  emm_df <- as.data.frame(
    summary(emm, type = "response")
  )
  
  # Summary for plotting (raw scale)
  summary_df <- data %>%
    group_by(treatment, extract_type, fraction_type) %>%
    summarise(
      mean = mean(fraction, na.rm = TRUE),
      se = sd(fraction, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  return(list(
    model = lmm,
    anova = anova_res,
    emmeans = emm_df,
    summary = summary_df,
    data = data
  ))
}


#run the statistical analyses
results <- run_fraction_logit_lmm(aggstab_all)

results$anova

# Generate CLD letters from the emmeans object
# Note: we use the 'emmeans' object, not the summary dataframe
# 1. Generate CLD for Extracts (comparing extracts within each fraction/treatment)
cld_res <- multcomp::cld(
  emmeans(results$model, ~ extract_type | fraction_type * treatment), 
  Letters = letters, 
  adjust = "tukey"
)

# 2. Clean up and convert to data frame
cld_df <- as.data.frame(cld_res, level = 0.9) %>%
  mutate(.group = trimws(.group))

# 3. Merge back to your summary_df
# Note: We use dplyr::select to avoid the 'unused arguments' error from earlier
results$summary <- results$summary %>%
  dplyr::select(-any_of(".group")) %>% # Remove old letters if they exist
  left_join(cld_df %>% 
              dplyr::select(treatment, extract_type, fraction_type, .group), 
            by = c("treatment", "extract_type", "fraction_type"))

#quick plot of raw values
ggplot(results$data,
       aes(x = fraction_type,
           y = fraction,
           colour = fraction_type)) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  facet_grid(treatment ~ extract_type) +
  theme(legend.position = "none")

#proper results plots function
plot_fraction_results <- function(summary_df) {
  
  ggplot(summary_df, 
         aes(x = fraction_type, 
             y = mean, 
             fill = fraction_type)) +
    
    geom_col(position = position_dodge(width = 0.7), 
             width = 0.6, 
             colour = "black") +
    
    geom_errorbar(
      aes(ymin = mean - se, 
          ymax = mean + se),
      width = 0.2,
      position = position_dodge(width = 0.7)
    ) +
    
    # ADDED: This layer places the letters above the error bars
    geom_text(
      aes(label = .group, 
          y = mean + se), 
      vjust = -0.5,           # Adjust this to move letters up/down
      size = 3.5, 
      fontface = "bold",
      position = position_dodge(width = 0.7)
    ) +
    
    facet_grid(treatment ~ extract_type) +
    
    scale_fill_viridis_d(option = "D", end = 0.9) +
    
    # Expand y-axis slightly so letters don't get cut off at the top
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    
    labs(
      x = NULL,
      y = "Aggregate size class (%)"
    ) +
    
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Run the updated plot
plot_fraction_results(results$summary)


#new plot for extracts side-by-side
# Create the plot
agg_plot <- ggplot(results$summary, aes(x = extract_type, y = mean, fill = extract_type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(0.7)) +
  # Facet_grid creates the matrix: Treatments in rows, Fractions in columns
  facet_grid(treatment ~ fraction_type) + 
  labs(
    title = "Aggregate Stability: Extract Comparison within Fractions and Treatments",
    x = "Extraction Method",
    y = "Mean Mass Fraction (%)",
    fill = "Extract"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom"
  )

# Display the plot
print(agg_plot)

#side-by-side view
ggplot(cld_res, aes(x = extract_type, y = emmean, fill = extract_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  # Using lower.CL and upper.CL for the error bars
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  # Add the letters above the upper confidence limit
  geom_text(aes(label = .group, y = upper.CL), vjust = -0.5, size = 3.5, fontface = "bold") +
  facet_grid(treatment ~ fraction_type) +
  labs(
    title = "Soil Aggregate Stability",
    subtitle = "Letters indicate significant differences (p < 0.05) within groups",
    x = "Extraction Method",
    y = "Back-transformed Proportion",
    fill = "Extract"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

#still want to check what this test looks like if only compare pooled basalts with control


#MWD test
# Function to compare extracts within each treatment
analyze_MWD_by_treatment <- function(df) {
  
  results_list <- list()
  plot_list <- list()
  
  for(tmt in unique(df$treatment)) {
    
    # Subset for this treatment
    MWD_df <- df %>% filter(treatment == tmt)
    
    # Ensure extract_type is factor
    MWD_df <- MWD_df %>% mutate(extract_type = factor(extract_type))
    
    # Check assumptions
    mwd_mod <- lm(MWD ~ extract_type, data = MWD_df)
    shapiro_p <- shapiro.test(residuals(mwd_mod))$p.value
    levene_p  <- car::leveneTest(MWD ~ extract_type, data = MWD_df)$`Pr(>F)`[1]
    
    if(shapiro_p > 0.05 & levene_p > 0.05) {
      
      # Parametric: ANOVA + Tukey
      aov_mwd <- aov(MWD ~ extract_type, data = MWD_df)
      emm_mwd <- emmeans(aov_mwd, ~ extract_type)
      
      letters_mwd <- multcomp::cld(
        emm_mwd,
        Letters = letters,
        adjust = "tukey"
      ) %>%
        as.data.frame() %>%
        select(extract_type, .group)
      
    } else {
      
      # Non-parametric: Kruskal-Wallis + Dunn
      dunn_mwd <- dunn_test(
        MWD_df,
        MWD ~ extract_type,
        p.adjust.method = "bonferroni"
      )
      
      letters_vec <- multcompView::multcompLetters(
        setNames(
          dunn_mwd$p.adj,
          paste(dunn_mwd$group1, dunn_mwd$group2, sep = "-")
        )
      )$Letters
      
      letters_mwd <- tibble(
        extract_type = names(letters_vec),
        .group = letters_vec
      )
    }
    
    # Summarise mean ± SE
    MWD_summary <- MWD_df %>%
      group_by(extract_type) %>%
      summarise(
        mean_MWD = mean(MWD, na.rm = TRUE),
        se_MWD = sd(MWD, na.rm = TRUE)/sqrt(n()),
        .groups = "drop"
      ) %>%
      left_join(letters_mwd, by = "extract_type") %>%
      mutate(treatment = tmt)
    
    results_list[[tmt]] <- list(
      summary = MWD_summary,
      shapiro_p = shapiro_p,
      levene_p  = levene_p
    )
    
    # Plot
    p <- ggplot(MWD_summary, aes(x = extract_type, y = mean_MWD, fill = extract_type)) +
      geom_col(width = 0.7, colour = "grey20") +
      geom_errorbar(aes(ymin = mean_MWD - se_MWD, ymax = mean_MWD + se_MWD),
                    width = 0.2, linewidth = 0.6) +
      geom_text(aes(label = .group, y = mean_MWD + se_MWD + 1),
                size = 4, fontface = "bold") +
      scale_fill_viridis_d(option = "D", end = 0.9) +
      labs(x = NULL, y = "Mean weight diameter (µm)", title = paste("Treatment:", tmt)) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot_list[[tmt]] <- p
  }
  
  return(list(results = results_list, plots = plot_list))
}
