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

#import data
soil_BD <- read.csv(here("csv_files", "soil_BD.csv"), header = TRUE)

# Helper vector for filtering treatments
valid_treatments <- tribble(
  ~Tmt,          ~App_rate,
  "Control",     "0",
  "Lime",        "02",
  "Bolsdorfer",  "50",
  "Eifelgold",   "50",
  "Huhnerberg",  "50"
)

# Filter helper function to filter by valid treatments
filter_valid_treatments <- function(df, date_col) {
  df %>%
    semi_join(valid_treatments, by = c("Tmt", "App_rate")) %>%
    arrange(!!sym(date_col))
}

# Data cleaning and transformation function 
clean_data <- function(df) {
  df %>%
    mutate(across(c(Treatment, Dose_resp_exp), as.factor)) %>%
    separate(Treatment, into = c("Tmt", "App_rate"), sep = "_") %>%
    mutate(
      App_rate = replace_na(App_rate, "0"),
      Tmt = as.factor(Tmt),
      Date = as.Date(Date, format = "%Y/%m/%d"),
    )
}

soil_BD <- clean_data(soil_BD)

soil_BD_summary <- soil_BD %>%
  group_by(Date, Tmt, App_rate) %>%
  summarise(
    mean_BD = mean(BD..g.cm.3., na.rm = TRUE),
    se_BD = sd(BD..g.cm.3., na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write.csv(soil_BD_summary,file=here::here("outputs","soil_BD_summary.csv"), row.names=FALSE)


# General plotting function for soil and rhizon data
plot_time_series <- function(df, date_col, yvar_mean, yvar_se, y_label, plot_title, ylim = NULL) {
  df_filtered <- filter_valid_treatments(df, date_col)
  
  ggplot(df_filtered, aes_string(x = date_col, y = yvar_mean, colour = "Tmt")) +
    geom_point(size = 4) +
    geom_errorbar(aes_string(ymin = paste0(yvar_mean, " - ", yvar_se),
                             ymax = paste0(yvar_mean, " + ", yvar_se)),
                  width = 0.1) +
    scale_x_date(date_labels = "%b %Y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = y_label, title = plot_title, colour = "Treatment") +
    {if (!is.null(ylim)) ylim(ylim)}
}


#BD plot time series
soil_BD_plot <- plot_time_series(soil_BD_summary, "Date", "mean_BD", "se_BD", "Soil BD g cm^-3", "Soil bulk density over time")
