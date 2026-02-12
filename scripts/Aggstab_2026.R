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
library(here)

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
  
  # Convert to factors  (FIXED small typo from your code)
  df[c("extract", "rep", "treatment")] <- 
    lapply(df[c("extract", "rep", "treatment")], as.factor)
  
  # Select columns
  df <- df %>%
    select(sample, extract, rep, recovery....,
           X...2mm.fraction,
           X...0.25mm.fraction,
           X...0.063mm.fraction,
           X...0.063mm.fraction.1,
           treatment)
  
  # Rename columns
  df <- df %>%
    rename(
      recovery = recovery....,
      `>2mm.fraction` = X...2mm.fraction,
      `>0.25mm.fraction` = X...0.25mm.fraction,
      `>0.063mm.fraction` = X...0.063mm.fraction,
      `<0.063mm.fraction` = X...0.063mm.fraction.1
    )
  
  # Pivot to long
  df_long <- df %>%
    pivot_longer(
      cols = c(`>2mm.fraction`,
               `>0.25mm.fraction`,
               `>0.063mm.fraction`,
               `<0.063mm.fraction`),
      names_to = "fraction_type",
      values_to = "fraction"
    ) %>%
    mutate(
      fraction_type = case_when(
        fraction_type == ">2mm.fraction" ~ ">2",
        fraction_type == ">0.25mm.fraction" ~ ">0.25",
        fraction_type == ">0.063mm.fraction" ~ ">0.063",
        fraction_type == "<0.063mm.fraction" ~ "<0.063"
      ),
      extract_type = factor(extract_label,
                            levels = c("h2o", "oxalate", "h2o2", "naoh"))
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

#check final df
str(aggstab_all)

