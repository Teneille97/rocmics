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

# data cleaning function

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
      MWD = 
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
# -------------------------------------------------
# REMOVE problematic sample from h2o dataset
# -------------------------------------------------

aggstab_h2o <- aggstab_h2o %>%
  filter(!(sample == 4 &
             rep == "b" &
             treatment == "Bolsdorfer_50"))
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

# MWD comparison of extracts with figures of merit
analyze_MWD_by_treatment <- function(df) {
  
  results_list <- list()
  plot_list <- list()
  
  for(tmt in unique(df$treatment)) {
    
    # Subset
    MWD_df <- df %>% 
      filter(treatment == tmt) %>%
      mutate(extract_type = factor(extract_type))
    
    # Assumption checks
    mwd_mod <- lm(MWD ~ extract_type, data = MWD_df)
    shapiro_p <- shapiro.test(residuals(mwd_mod))$p.value
    levene_p  <- car::leveneTest(MWD ~ extract_type, data = MWD_df)$`Pr(>F)`[1]
    
    # -----------------------------
    # PARAMETRIC
    # -----------------------------
    if(shapiro_p > 0.05 & levene_p > 0.05) {
      
      test_type <- "ANOVA"
      
      aov_mwd <- aov(MWD ~ extract_type, data = MWD_df)
      aov_tab <- summary(aov_mwd)[[1]]
      
      F_value  <- aov_tab$`F value`[1]
      p_value  <- aov_tab$`Pr(>F)`[1]
      df1      <- aov_tab$Df[1]
      df2      <- aov_tab$Df[2]
      
      emm_mwd <- emmeans(aov_mwd, ~ extract_type)
      
      letters_mwd <- multcomp::cld(
        emm_mwd,
        Letters = letters,
        adjust = "tukey"
      ) %>%
        as.data.frame() %>%
        dplyr::select(extract_type, .group)
      
    } else {
      
      # -----------------------------
      # NON-PARAMETRIC
      # -----------------------------
      
      test_type <- "Kruskal-Wallis"
      
      kw_test <- kruskal.test(MWD ~ extract_type, data = MWD_df)
      
      F_value <- kw_test$statistic
      p_value <- kw_test$p.value
      df1     <- kw_test$parameter
      df2     <- NA
      
      dunn_mwd <- rstatix::dunn_test(
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
    
    # -----------------------------
    # Summary table
    # -----------------------------
    
    MWD_summary <- MWD_df %>%
      group_by(extract_type) %>%
      summarise(
        mean_MWD = mean(MWD, na.rm = TRUE),
        se_MWD   = sd(MWD, na.rm = TRUE)/sqrt(n()),
        .groups = "drop"
      ) %>%
      left_join(letters_mwd, by = "extract_type") %>%
      mutate(treatment = tmt)
    
    # Store full results
    results_list[[tmt]] <- list(
      summary    = MWD_summary,
      test_type  = test_type,
      statistic  = F_value,
      df_between = df1,
      df_within  = df2,
      p_value    = p_value,
      shapiro_p  = shapiro_p,
      levene_p   = levene_p
    )
    
    # -----------------------------
    # Plot
    # -----------------------------
    
    p_label <- if(test_type == "ANOVA") {
      paste0("ANOVA: F(", df1, ", ", df2, ") = ",
             round(F_value, 2),
             ", p = ",
             signif(p_value, 3))
    } else {
      paste0("Kruskal-Wallis: χ²(", df1, ") = ",
             round(F_value, 2),
             ", p = ",
             signif(p_value, 3))
    }
    
    p <- ggplot(MWD_summary, aes(x = extract_type, y = mean_MWD, fill = extract_type)) +
      geom_col(width = 0.7, colour = "grey20") +
      geom_errorbar(aes(ymin = mean_MWD - se_MWD,
                        ymax = mean_MWD + se_MWD),
                    width = 0.2, linewidth = 0.6) +
      geom_text(aes(label = .group,
                    y = mean_MWD + se_MWD + 1),
                size = 4, fontface = "bold") +
      annotate("text",
               x = Inf, y = Inf,
               label = p_label,
               hjust = 1.1, vjust = 1.5,
               size = 3.5) +
      scale_fill_viridis_d(option = "D", end = 0.9) +
      labs(x = NULL,
           y = "Mean weight diameter (µm)",
           title = paste("Treatment:", tmt)) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot_list[[tmt]] <- p
  }
  
  return(list(results = results_list, plots = plot_list))
}

# Run function
MWD_analysis <- analyze_MWD_by_treatment(aggstab_all)

# Combine all treatment summaries into one dataframe
MWD_all_summary <- bind_rows(
  lapply(MWD_analysis$results, function(x) x$summary)
)

# Make treatment a factor (preserve original order)
MWD_all_summary$treatment <- factor(MWD_all_summary$treatment, 
                                    levels = unique(MWD_all_summary$treatment))

# Facetted plot
ggplot(MWD_all_summary, aes(x = extract_type, y = mean_MWD, fill = extract_type)) +
  geom_col(width = 0.7, colour = "grey20") +
  geom_errorbar(aes(ymin = mean_MWD - se_MWD, ymax = mean_MWD + se_MWD),
                width = 0.2, linewidth = 0.6) +
  geom_text(aes(label = .group, y = mean_MWD + se_MWD + 1),
            size = 4, fontface = "bold") +
  facet_wrap(~ treatment, scales = "free_y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(x = NULL, y = "Mean weight diameter (µm)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

#access individual test results like so:
MWD_analysis$results[["Bolsdorfer_50"]]$summary

##access individual plot facets like so: 
MWD_analysis$plots

#function to compare mwd of tmts with figures of merit
analyze_MWD_by_extract <- function(df) {
  
  results_list <- list()
  plot_list <- list()
  
  for(xtract in unique(df$extract_type)) {
    
    # Subset
    MWD_df <- df %>% 
      filter(extract_type == xtract) %>%
      mutate(treatment = factor(treatment))
    
    # Assumption checks
    mwd_mod <- lm(MWD ~ treatment, data = MWD_df)
    shapiro_p <- shapiro.test(residuals(mwd_mod))$p.value
    levene_p  <- car::leveneTest(MWD ~ treatment, data = MWD_df)$`Pr(>F)`[1]
    
    # -----------------------------
    # PARAMETRIC
    # -----------------------------
    if(shapiro_p > 0.05 & levene_p > 0.05) {
      
      test_type <- "ANOVA"
      
      aov_mwd <- aov(MWD ~ treatment, data = MWD_df)
      aov_tab <- summary(aov_mwd)[[1]]
      
      F_value  <- aov_tab$`F value`[1]
      p_value  <- aov_tab$`Pr(>F)`[1]
      df1      <- aov_tab$Df[1]
      df2      <- aov_tab$Df[2]
      
      emm_mwd <- emmeans(aov_mwd, ~ treatment)
      
      letters_mwd <- multcomp::cld(
        emm_mwd,
        Letters = letters,
        adjust = "tukey"
      ) %>%
        as.data.frame() %>%
        dplyr::select(treatment, .group)
      
    } else {
      
      # -----------------------------
      # NON-PARAMETRIC
      # -----------------------------
      
      test_type <- "Kruskal-Wallis"
      
      kw_test <- kruskal.test(MWD ~ treatment, data = MWD_df)
      
      F_value <- kw_test$statistic
      p_value <- kw_test$p.value
      df1     <- kw_test$parameter
      df2     <- NA
      
      dunn_mwd <- rstatix::dunn_test(
        MWD_df,
        MWD ~ treatment,
        p.adjust.method = "bonferroni"
      )
      
      letters_vec <- multcompView::multcompLetters(
        setNames(
          dunn_mwd$p.adj,
          paste(dunn_mwd$group1, dunn_mwd$group2, sep = "-")
        )
      )$Letters
      
      letters_mwd <- tibble(
        treatment = names(letters_vec),
        .group = letters_vec
      )
    }
    
    # -----------------------------
    # Summary table
    # -----------------------------
    
    MWD_summary <- MWD_df %>%
      group_by(treatment) %>%
      summarise(
        mean_MWD = mean(MWD, na.rm = TRUE),
        se_MWD   = sd(MWD, na.rm = TRUE)/sqrt(n()),
        .groups = "drop"
      ) %>%
      left_join(letters_mwd, by = "treatment") %>%
      mutate(extract_type = xtract)
    
    # Store full results
    results_list[[xtract]] <- list(
      summary    = MWD_summary,
      test_type  = test_type,
      statistic  = F_value,
      df_between = df1,
      df_within  = df2,
      p_value    = p_value,
      shapiro_p  = shapiro_p,
      levene_p   = levene_p
    )
    
    # -----------------------------
    # Plot
    # -----------------------------
    
    p_label <- if(test_type == "ANOVA") {
      paste0("ANOVA: F(", df1, ", ", df2, ") = ",
             round(F_value, 2),
             ", p = ",
             signif(p_value, 3))
    } else {
      paste0("Kruskal-Wallis: χ²(", df1, ") = ",
             round(F_value, 2),
             ", p = ",
             signif(p_value, 3))
    }
    
    p <- ggplot(MWD_summary, aes(x = treatment, y = mean_MWD, fill = treatment)) +
      geom_col(width = 0.7, colour = "grey20") +
      geom_errorbar(aes(ymin = mean_MWD - se_MWD,
                        ymax = mean_MWD + se_MWD),
                    width = 0.2, linewidth = 0.6) +
      geom_text(aes(label = .group,
                    y = mean_MWD + se_MWD + 1),
                size = 4, fontface = "bold") +
      annotate("text",
               x = Inf, y = Inf,
               label = p_label,
               hjust = 1.1, vjust = 1.5,
               size = 3.5) +
      scale_fill_viridis_d(option = "D", end = 0.9) +
      labs(x = NULL,
           y = "Mean weight diameter (µm)",
           title = paste("Extract:", xtract)) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot_list[[xtract]] <- p
  }
  
  return(list(results2 = results_list, plots2 = plot_list))
}

# Run function
MWD_analysis2 <- analyze_MWD_by_extract(aggstab_all)

# Combine all treatment summaries into one dataframe
MWD_all_summary2 <- bind_rows(
  lapply(MWD_analysis2$results2, function(x) x$summary)
)

# Make treatment a factor (preserve original order)
MWD_all_summary2$treatment <- factor(MWD_all_summary2$treatment, 
                                    levels = unique(MWD_all_summary2$treatment))

# Facetted plot
ggplot(MWD_all_summary2, aes(x = treatment, y = mean_MWD, fill = treatment)) +
  geom_col(width = 0.7, colour = "grey20") +
  geom_errorbar(aes(ymin = mean_MWD - se_MWD, ymax = mean_MWD + se_MWD),
                width = 0.2, linewidth = 0.6) +
  geom_text(aes(label = .group, y = mean_MWD + se_MWD + 1),
            size = 4, fontface = "bold") +
  facet_wrap(~ extract_type, scales = "free_y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(x = NULL, y = "Mean weight diameter (µm)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )


#access individual test results like so:
MWD_analysis2$results2[["h2o"]]$summary

#access individual plot facets like so: 
MWD_analysis2$plots2

#load soil pH, EC, rhizon and xrf files
Alldata_Soil_phEC_summary <- read.csv(here("outputs", "Alldata_Soil_phEC_summary.csv"), header = TRUE)
Alldata_Soil_phEC_Oct25 <- Alldata_Soil_phEC_summary %>% filter(Date == "2025-10-28")
Alldata_Rhizon_summary <- read.csv(here("outputs", "Alldata_Rhizon_summary.csv"), header = TRUE)



#format aggstab files to have same format as ph ec
aggstab_h2o_formatpH <- aggstab_h2o %>%
separate(treatment, into = c("Tmt", "App_rate"), sep = "_") %>%
  mutate(
    App_rate = replace_na(App_rate, "0"),
    Tmt = as.factor(Tmt)) %>%
  group_by(Tmt, App_rate, extract_type) %>%
  summarise(across(
    c(MWD),
    list(mean = ~ mean(.x, na.rm = TRUE),
         se   = ~ sd(.x, na.rm = TRUE) / sqrt(n())),
    .names = "{.fn}_{.col}"),
    .groups = "drop")


cor(Alldata_Soil_phEC_Oct25$mean_EC, aggstab_h2o_formatpH$mean_MWD) #no cor
cor(Alldata_Soil_phEC_Oct25$mean_pH, aggstab_h2o_formatpH$mean_MWD) #no cor

#check cor composition
