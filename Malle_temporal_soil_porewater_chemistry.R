#load packages
library(ggplot2)
library(gganimate)
library(viridis)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(car)
library(lme4)
library(zoo)
library(ggthemr)
theme_update(
  legend.title = element_text(size = 10),
  legend.text  = element_text(size = 8)
)
library(gridExtra)

#set theme
ggthemr('flat dark')

library(tidyverse)

# Import files in one step
Alldata_Rhizon <- read.csv("csv_files/Alldata_Rhizon.csv", header = TRUE)[1:252, ]
Alldata_Soil_phEC <- read.csv("csv_files/Alldata_Soil_phEC.csv", header = TRUE)
Alldata_PRS <- read.csv("csv_files/Alldata_PRS.csv", header = TRUE)

# Helper vector for filtering treatments
valid_treatments <- tribble(
  ~Tmt,          ~App_rate,
  "Control",     "0",
  "Lime",        "02",
  "Bolsdorfer",  "50",
  "Eifelgold",   "50",
  "Huhnerberg",  "50"
)

# Data cleaning and transformation function for Alldata_Rhizon
clean_rhizon <- function(df) {
  df %>%
    mutate(across(c(Treatment, Dose_response_exp.), as.factor)) %>%
    separate(Treatment, into = c("Tmt", "App_rate"), sep = "_") %>%
    mutate(
      App_rate = replace_na(App_rate, "0"),
      Tmt = as.factor(Tmt),
      Sampled_on = as.Date(Sampled_on, format = "%Y/%m/%d")
    ) %>%
    # Convert columns from 'Alkalinity_mg_l_CaCO3' onward to numeric after trimming whitespace
    { 
      start_col <- which(names(.) == "Alkalinity_mg_l_CaCO3")
      .[start_col:ncol(.)] <- lapply(.[start_col:ncol(.)], function(x) as.numeric(trimws(x)))
      .
    }
}

Alldata_Rhizon <- clean_rhizon(Alldata_Rhizon)

# Summarise Alldata_Rhizon
Alldata_Rhizon_summary <- Alldata_Rhizon %>%
  group_by(Sampled_on, Tmt, App_rate) %>%
  summarise(across(
    c(pH, EC_µS_cm, DOC, DIC, Alkalinity_meq_l, NO2.N_mgN_l, NH4.N_mgN_l,
      NO3.N_mgN_l, PO4.P_mgP_l, Cl_mgCl_l, Ca_mg_l, Fe_mg_l, K_mg_l, 
      Mg_mg_l, Na_mg_l, Ni_mg_l, P_mg_l, Al_mg_l, Si_mg_l),
    list(mean = ~ mean(.x, na.rm = TRUE),
         se   = ~ sd(.x, na.rm = TRUE) / sqrt(n())),
    .names = "{.fn}_{.col}"),
    .groups = "drop")

# Clean and transform Alldata_Soil_phEC similarly
Alldata_Soil_phEC <- Alldata_Soil_phEC %>%
  mutate(across(c(Pot, Plot, Treatment, Dose_response_exp.), as.factor)) %>%
  separate(Treatment, into = c("Tmt", "App_rate"), sep = "_") %>%
  mutate(
    App_rate = replace_na(App_rate, "0"),
    Date = as.Date(Date, format = "%Y/%m/%d"),
    Tmt = as.factor(Tmt)
  )

Alldata_Soil_phEC_summary <- Alldata_Soil_phEC %>%
  group_by(Date, Tmt, App_rate) %>%
  summarise(
    mean_pH = mean(pH, na.rm = TRUE),
    se_pH = sd(pH, na.rm = TRUE) / sqrt(n()),
    mean_EC = mean(EC, na.rm = TRUE),
    se_EC = sd(EC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Filter helper function to filter by valid treatments
filter_valid_treatments <- function(df, date_col) {
  df %>%
    semi_join(valid_treatments, by = c("Tmt", "App_rate")) %>%
    arrange(!!sym(date_col))
}

# General plotting function for soil and rhizon data
plot_time_series <- function(df, date_col, yvar_mean, yvar_se, y_label, plot_title, ylim = NULL) {
  df_filtered <- filter_valid_treatments(df, date_col)
  
  ggplot(df_filtered, aes_string(x = date_col, y = yvar_mean, colour = "Tmt")) +
    geom_point(size = 4) +
    geom_errorbar(aes_string(ymin = paste0(yvar_mean, " - ", yvar_se),
                             ymax = paste0(yvar_mean, " + ", yvar_se)),
                  width = 0.1) +
    scale_x_date(breaks = "3 months", date_labels = "%b %Y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = y_label, title = plot_title, colour = "Treatment") +
    {if (!is.null(ylim)) ylim(ylim)}
}

# Create plots with the function
soilpH_plot <- plot_time_series(Alldata_Soil_phEC_summary, "Date", "mean_pH", "se_pH", "Soil pH", "Soil pH over time")
soilEC_plot <- plot_time_series(Alldata_Soil_phEC_summary, "Date", "mean_EC", "se_EC", "Soil EC (uS/cm)", "Soil EC over time")

rhizon_pH_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_pH", "se_pH", "Porewater pH", "Porewater pH over time")
rhizon_EC_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_EC_µS_cm", "se_EC_µS_cm", "Porewater EC (uS/cm)", "Porewater EC over time (no outliers 2nd sampling)", ylim = c(0, 750))
rhizon_DIC_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_DIC", "se_DIC", "Porewater DIC (ppm)", "Porewater DIC over time")
rhizon_DOC_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_DOC", "se_DOC", "Porewater DOC (ppm)", "Porewater DOC over time")
rhizon_alkalinity_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_Alkalinity_meq_l", "se_Alkalinity_meq_l", "Porewater alkalinity (meq/L)", "Porewater alkalinity over time")
rhizon_Ca_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_Ca_mg_l", "se_Ca_mg_l", "Porewater Ca (mg/L)", "Porewater Ca over time")
rhizon_Mg_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_Mg_mg_l", "se_Mg_mg_l", "Porewater Mg (mg/L)", "Porewater Mg over time")
rhizon_K_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_K_mg_l", "se_K_mg_l", "Porewater K (mg/L)", "Porewater K over time")
rhizon_Na_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_Na_mg_l", "se_Na_mg_l", "Porewater Na (mg/L)", "Porewater Na over time")
rhizon_Fe_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_Fe_mg_l", "se_Fe_mg_l", "Porewater Fe (mg/L)", "Porewater Fe over time")
rhizon_Al_plot <- plot_time_series(Alldata_Rhizon_summary, "Sampled_on", "mean_Al_mg_l", "se_Al_mg_l", "Porewater Al (mg/L)", "Porewater Al over time")


#data cleaning PRS
Alldata_PRS<-Alldata_PRS[,1:22]

Alldata_PRS <- Alldata_PRS %>%
  separate(
    col = Treatment,       
    into = c("Tmt", "App_rate"), 
    sep = "_"              
  ) %>%
  mutate(
    App_rate = replace_na(App_rate, "0") 
  )

Alldata_PRS$Retrieval.Date<-as.Date(Alldata_PRS$Retrieval.Date,format = "%Y/%m/%d") 
Alldata_PRS$Tmt<-as.factor(Alldata_PRS$Tmt)

# Label generator

make_lab <- function(element) {
  as.expression(
    bquote(.(element) ~ "supply rate (" * mu * "g / 10 cm"^2 * "/ week)")
  )
}

elements <- c("Ca","Mg","K","Fe","Mn","Cu","Zn","B","S","Cd","Pb")
labels   <- setNames(lapply(elements, make_lab), elements)

# Common filter 

PRS_filtered <- Alldata_PRS %>%
  dplyr::filter(
    (Tmt == "Control" & App_rate == "0") |
      (Tmt == "Lime"    & App_rate == "02") |
      (Tmt %in% c("Bolsdorfer","Eifelgold","Huhnerberg") & App_rate == "50")
  )

# Plot function

make_prs_plot <- function(df, element) {
  ggplot(df, aes(x = Retrieval.Date, y = .data[[element]], colour = Tmt)) +
    geom_point(size = 4) +
    scale_x_date(breaks = "3 months", date_labels = "%b %Y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(
      y = labels[[element]],
      x = "Date",
      title = paste(element, "supply rate over time"),
      colour = "Treatment"
    )
}

# Generate all plots in a list

PRS_plots <- lapply(elements, \(el) make_prs_plot(PRS_filtered, el))
names(PRS_plots) <- elements


plots <- list(soilpH_plot, soilEC_plot, rhizon_pH_plot, rhizon_EC_plot, rhizon_DIC_plot, rhizon_DOC_plot, rhizon_alkalinity_plot, rhizon_Ca_plot, rhizon_Mg_plot, rhizon_Na_plot, rhizon_K_plot, rhizon_Al_plot, rhizon_Fe_plot,
              PRS_plots$Ca, PRS_plots$Mg, PRS_plots$K, PRS_plots$Fe, PRS_plots$Mn, PRS_plots$Cu, PRS_plots$Zn, PRS_plots$B, PRS_plots$S, PRS_plots$Pb)  

#print in console
for (p in plots) print(p)

#render with quarto
multi_page <- marrangeGrob(
  grobs = plots,
  nrow = 1,
  ncol = 1,
  top = NULL,      # removes “page x” title inside gridExtra
  layout_matrix = matrix(1) # ensures full-page plot
)
multi_page
