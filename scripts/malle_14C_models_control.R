# SOM decomposition models for Malle soil

## load libraries
library(here)
library(SoilR)
library(FME)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(emmeans)

## import data
radiocarbon_samples <- read.csv(here("outputs", "radiocarbon_samples.csv"))
soil_BD_summary<-read.csv(here("outputs","soil_BD_summary.csv"))

## data cleaning
radiocarbon_samples <- radiocarbon_samples %>%
  mutate(across(c("SOM_fraction", "Tmt", "App_rate"), as.factor))  %>% 
  filter(Tmt != "Huhnerberg" | App_rate == 50) #drop dose-response

soil_BD_summary_Feb2026<-soil_BD_summary[soil_BD_summary$Date=="2026-02-28",] %>% #keep only Feb 2026 BD 
  mutate(across(c("Tmt", "App_rate"), as.factor)) %>% 
  filter(Tmt != "Huhnerberg" | App_rate == 50) #drop dose-response

radiocarbon_samples <- radiocarbon_samples %>%
  left_join(
    soil_BD_summary_Feb2026 %>%
      select(Tmt, App_rate, mean_BD) %>%
      rename(bulk_density = mean_BD),
    by = c("Tmt", "App_rate")
  )

names(radiocarbon_samples)[names(radiocarbon_samples) == 'Sampling_year'] <- 'Year'

# CHOOSE TREATMENT
radiocarbon_samples  %>% 
  filter(Tmt == "Control") 

## inputs
mean_C_inputs <- 100 # g m-2

## Atmospheric radiocarbon ---------------------------------------------------

Atm14C <- Hua2021$NHZone1[,1:2]
fAtm14C <- read.csv(here("csv_files", "NHZ1forecast.csv"))
Atm14C <- rbind(
  Atm14C,
  data.frame(
    Year = fAtm14C$time,
    mean.Delta14C = Delta14C_from_AbsoluteFractionModern(fAtm14C$F14C)
  )
)

# add C and N stocks
radiocarbon_samples$C_stocks_gm2 <- radiocarbon_samples$Ctotal * radiocarbon_samples$bulk_density * 0.2 * 10000 * (radiocarbon_samples$fraction_prop/100)
radiocarbon_samples$N_stocks_gm2 <- radiocarbon_samples$Ntotal * radiocarbon_samples$bulk_density * 0.2 * 10000 * (radiocarbon_samples$fraction_prop/100)

summary_df <- radiocarbon_samples %>%
  group_by(Year, SOM_fraction) %>%
  summarise(
    C_stocks_mean = mean(C_stocks_gm2, na.rm = TRUE),
    C_stocks_sd   = sd(C_stocks_gm2, na.rm = TRUE),
    N_stocks_mean = mean(N_stocks_gm2, na.rm = TRUE),
    N_stocks_sd   = sd(N_stocks_gm2, na.rm = TRUE),
    d14C_mean     = mean(X.14C....., na.rm = TRUE),
    d14C_sd       = sd(X.14C....., na.rm = TRUE),
    CN_mean       = mean(CN, na.rm = TRUE),
    CN_sd       = sd(CN, na.rm = TRUE),
    .groups = "drop"
  )

bulk <- subset(summary_df, SOM_fraction == "Bulk")
fPOM <- subset(summary_df, SOM_fraction == "fPOM")
oPOM<- subset(summary_df, SOM_fraction == "oPOM")
MAOM<- subset(summary_df, SOM_fraction == "MAOM")

# obs data for cost func
Cobs_bulk <- data.frame(Year = bulk$Year, Ct = bulk$C_stocks_mean, Ct_sd = bulk$C_stocks_sd)
C14obs_bulk <- data.frame(Year = bulk$Year, C14t = bulk$d14C_mean, C14t_sd = bulk$d14C_sd)
Cobs_fPOM <- data.frame(Year = fPOM$Year, Ct_fPOM = fPOM$C_stocks_mean, Ct_fPOM_sd = fPOM$C_stocks_sd)
C14obs_fPOM <-data.frame(Year = fPOM$Year, C14t_fPOM = fPOM$d14C_mean, C14t_fPOM_sd = fPOM$d14C_sd)
Cobs_oPOM<- data.frame(Year = oPOM$Year, Ct_oPOM = oPOM$C_stocks_mean, Ct_oPOM_sd = oPOM$C_stocks_sd)
C14obs_oPOM <-data.frame(Year = oPOM$Year, C14t_oPOM = oPOM$d14C_mean, C14t_oPOM_sd = oPOM$d14C_sd)
Cobs_MAOM <- data.frame(Year = MAOM$Year, Ct_MAOM = MAOM$C_stocks_mean, Ct_MAOM_sd = MAOM$C_stocks_sd)
C14obs_MAOM <- data.frame(Year = MAOM$Year, C14t_MAOM = MAOM$d14C_mean, C14t_MAOM_sd = MAOM$d14C_sd)

# obs N stocks (not for models, only final age distributions)
Nobs_bulk <- data.frame(Year = bulk$Year, Nt = bulk$N_stocks_mean, Nt_sd = bulk$N_stocks_sd)
Nobs_fPOM <- data.frame(Year = fPOM$Year, Nt_fPOM = fPOM$N_stocks_mean, Nt_fPOM_sd = fPOM$N_stocks_sd)
Nobs_oPOM<- data.frame(Year = oPOM$Year, Nt_oPOM = oPOM$N_stocks_mean, Nt_oPOM_sd = oPOM$N_stocks_sd)
Nobs_MAOM <- data.frame(Year = MAOM$Year, Nt_MAOM = MAOM$N_stocks_mean, Nt_MAOM_sd = MAOM$N_stocks_sd)

# initial values
yr <- seq(2023, 2026, by = 1/12)
C0_bulk <- mean(Cobs_bulk[Cobs_bulk$Year==2023,]$Ct)
C0_fPOM <-  mean(Cobs_fPOM[Cobs_fPOM$Year==2023,]$Ct_fPOM)
C0_oPOM <-  mean(Cobs_oPOM[Cobs_oPOM$Year==2023,]$Ct_oPOM)
C0_MAOM <-  mean(Cobs_MAOM[Cobs_MAOM$Year==2023,]$Ct_MAOM)
F0_bulk <- mean(C14obs_bulk[C14obs_bulk$Year==2023,]$C14t) 
F0_fPOM  <- mean(C14obs_fPOM[C14obs_fPOM$Year==2023,]$C14t_fPOM) 
F0_oPOM  <- mean(C14obs_oPOM[C14obs_oPOM$Year==2023,]$C14t_oPOM) 
F0_MAOM  <- mean(C14obs_MAOM[C14obs_MAOM$Year==2023,]$C14t_MAOM)

# func to run mod
run_mod <- function(pars){ #kf, ki, ks, alpha 21, alpha 32
  c14_atm <- BoundFc(Atm14C, format = "Delta14C")
  c14_initial <- ConstFc(
    values = c(F0_fPOM, F0_oPOM, F0_MAOM), 
    format = "Delta14C"
  )
  A3 <- diag(-c(pars[1:3]))
  A3[2,1] <- pars[1]*pars[4]
  A3[3,2] <- pars[2]*pars[5]
  mod<-GeneralModel_14(
    t = yr,
    A = A3,
    ivList = c(C0_fPOM, C0_oPOM, C0_MAOM),
    initialValF = c14_initial,
    inputFluxes = c(0.8*mean_C_inputs,0.2*mean_C_inputs,0),
    inputFc = c14_atm
  )
  
  Ct_pools <- getC(mod)
  C14_pools <- getF14(mod)
  C14t <- getF14C(mod)
  
  mod_results<-data.frame(
    Year = yr,
    Ct = rowSums(Ct_pools),
    C14t = C14t,
    Ct_fPOM = Ct_pools[,1],
    Ct_oPOM = Ct_pools[,2],
    Ct_MAOM = Ct_pools[,3],
    C14t_fPOM = C14_pools[,1],
    C14t_oPOM = C14_pools[,2],
    C14t_MAOM = C14_pools[,3]
  )
  
  return(mod_results)
}

inipars <- c(0.1, 0.05, 0.001, 0.05, 0.05) #kf, ki, ks, alpha 21, alpha 32

# cost func
mc <- function(pars){
  out = run_mod(pars)
  Cost1 <- modCost(out, Cobs_bulk, x = "Year", err = "Ct_sd")
  Cost2 <- modCost(out, C14obs_bulk, x = "Year", cost = Cost1, err = "C14t_sd")
  Cost3 <- modCost(out, Cobs_fPOM, x = "Year", cost = Cost2, err = "Ct_fPOM_sd")
  Cost4 <- modCost(out, C14obs_fPOM, x = "Year", cost = Cost3, err = "C14t_fPOM_sd")
  Cost5 <- modCost(out, Cobs_oPOM, x = "Year", cost = Cost4, err = "Ct_oPOM_sd")
  Cost6 <- modCost(out, C14obs_oPOM, x = "Year", cost = Cost5, err = "C14t_oPOM_sd")
  Cost7 <- modCost(out, Cobs_MAOM, x = "Year", cost = Cost6, err = "Ct_MAOM_sd")
  modCost(out, C14obs_MAOM, x = "Year", cost = Cost7, err = "C14t_MAOM_sd")
}                                                                                                                                                                                                                                  

#fit mod
#mFit <- modFit(
#  f = mc,
#  p = inipars, 
#  method = "Nelder-Mead",
#  upper = c(1.5, 0.5 ,0.01, 0.5, 0.5),
#  lower = c(0.05, 0.005, 0.0001, 0, 0) 
#)

#save(mFit, file = file.path("mod_runs", "mFit_3p_Control.Rdata"))
load(here::here("mod_runs/mFit_3p_Control.Rdata"))


# create mod with best pars
bestpars <- mFit$par
out_best <- run_mod(bestpars)


# plots
Delta14Clabel <- expression(Delta^14*C)
stocks_label <- expression(C ~ (g ~ m^{-2}))

col_bulk <- "black"
col_fPOM <- "#1b9e77"
col_oPOM <- "#0000FF"
col_MAOM <- "#FF0000"

plot_C <- ggplot() +
  
  ## --- observed uncertainty ribbons FIRST ---
  geom_ribbon(
    data = Cobs_bulk,
    aes(
      x = Year,
      ymin = Ct - Ct_sd,
      ymax = Ct + Ct_sd,
      fill = "Bulk"
    ),
    alpha = 0.15
  ) +
  
  geom_ribbon(
    data = Cobs_fPOM,
    aes(
      x = Year,
      ymin = Ct_fPOM - Ct_fPOM_sd,
      ymax = Ct_fPOM + Ct_fPOM_sd,
      fill = "fPOM"
    ),
    alpha = 0.15
  ) +
  
  geom_ribbon(
    data = Cobs_MAOM,
    aes(
      x = Year,
      ymin = Ct_MAOM - Ct_MAOM_sd,
      ymax = Ct_MAOM + Ct_MAOM_sd,
      fill = "MAOM"
    ),
    alpha = 0.15
  ) +
  
  geom_ribbon(
    data = Cobs_oPOM,
    aes(
      x = Year,
      ymin = Ct_oPOM - Ct_oPOM_sd,
      ymax = Ct_oPOM + Ct_oPOM_sd,
      fill = "oPOM"
    ),
    alpha = 0.15
  ) +
  
  
  ## --- modelled lines ---
  geom_line(
    data = out_best,
    aes(x = Year, y = Ct, colour = "Bulk", linetype = "Bulk"),
    linewidth = 1
  ) +
  
  geom_line(
    data = out_best,
    aes(x = Year, y = Ct_fPOM, colour = "fPOM", linetype = "fPOM"),
    linewidth = 0.8
  ) +
  
  geom_line(
    data = out_best,
    aes(x = Year, y = Ct_MAOM, colour = "MAOM", linetype = "MAOM"),
    linewidth = 0.8
  ) +
  
  geom_line(
    data = out_best,
    aes(x = Year, y = Ct_oPOM, colour = "oPOM", linetype = "oPOM"),
    linewidth = 0.8
  ) +
  
  
  ## --- observed points ---
  geom_point(
    data = Cobs_bulk,
    aes(x = Year, y = Ct, colour = "Bulk"),
    size = 2
  ) +
  
  geom_point(
    data = Cobs_fPOM,
    aes(x = Year, y = Ct_fPOM, colour = "fPOM"),
    size = 2
  ) +
  
  geom_point(
    data = Cobs_MAOM,
    aes(x = Year, y = Ct_MAOM, colour = "MAOM"),
    size = 2
  ) +
  
  geom_point(
    data = Cobs_oPOM,
    aes(x = Year, y = Ct_oPOM, colour = "oPOM"),
    size = 2
  ) +
  
  
  ## --- colour scale ---
  scale_colour_manual(
    name = "Pool",
    values = c(
      "Bulk" = col_bulk,
      "fPOM" = col_fPOM,
      "MAOM" = col_MAOM,
      "oPOM" = col_oPOM
    )
  ) +
  
  
  ## --- fill scale ---
  scale_fill_manual(
    name = "Pool",
    values = c(
      "Bulk" = col_bulk,
      "fPOM" = col_fPOM,
      "MAOM" = col_MAOM,
      "oPOM" = col_oPOM
    ),
    guide = "none"
  ) +
  
  
  ## --- line types ---
  scale_linetype_manual(
    name = "Pool",
    values = c(
      "Bulk" = "solid",
      "fPOM" = "dotted",
      "MAOM" = "dotted",
      "oPOM" = "dotted"
    ),
    guide = "none"
  ) +
  
  
  ## --- labels ---
  ggtitle("C stocks") +
  ylab(stocks_label) +
  xlab("Year") +
  
  
  ## --- theme ---
  theme_minimal()

plot_14C <- ggplot() +
  
  ## --- observed uncertainty ribbons FIRST ---
  geom_ribbon(
    data = C14obs_bulk,
    aes(
      x = Year,
      ymin = C14t - C14t_sd,
      ymax = C14t + C14t_sd,
      fill = "Bulk"
    ),
    alpha = 0.15
  ) +
  
  geom_ribbon(
    data = C14obs_fPOM,
    aes(
      x = Year,
      ymin = C14t_fPOM - C14t_fPOM_sd,
      ymax = C14t_fPOM + C14t_fPOM_sd,
      fill = "fPOM"
    ),
    alpha = 0.15
  ) +
  
  geom_ribbon(
    data = C14obs_MAOM,
    aes(
      x = Year,
      ymin = C14t_MAOM - C14t_MAOM_sd,
      ymax = C14t_MAOM + C14t_MAOM_sd,
      fill = "MAOM"
    ),
    alpha = 0.15
  ) +
  
  geom_ribbon(
    data = C14obs_oPOM,
    aes(
      x = Year,
      ymin = C14t_oPOM - C14t_oPOM_sd,
      ymax = C14t_oPOM + C14t_oPOM_sd,
      fill = "oPOM"
    ),
    alpha = 0.15
  ) +
  
  
  ## --- modelled lines ---
  geom_line(
    data = out_best,
    aes(
      x = Year,
      y = C14t,
      colour = "Bulk",
      linetype = "Bulk"
    ),
    linewidth = 1
  ) +
  
  geom_line(
    data = out_best,
    aes(
      x = Year,
      y = C14t_fPOM,
      colour = "fPOM",
      linetype = "fPOM"
    ),
    linewidth = 0.8
  ) +
  
  geom_line(
    data = out_best,
    aes(
      x = Year,
      y = C14t_MAOM,
      colour = "MAOM",
      linetype = "MAOM"
    ),
    linewidth = 0.8
  ) +
  
  geom_line(
    data = out_best,
    aes(
      x = Year,
      y = C14t_oPOM,
      colour = "oPOM",
      linetype = "oPOM"
    ),
    linewidth = 0.8
  ) +
  
  
  ## --- observed points ---
  geom_point(
    data = C14obs_bulk,
    aes(x = Year, y = C14t, colour = "Bulk"),
    size = 2
  ) +
  
  geom_point(
    data = C14obs_fPOM,
    aes(x = Year, y = C14t_fPOM, colour = "fPOM"),
    size = 2
  ) +
  
  geom_point(
    data = C14obs_MAOM,
    aes(x = Year, y = C14t_MAOM, colour = "MAOM"),
    size = 2
  ) +
  
  geom_point(
    data = C14obs_oPOM,
    aes(x = Year, y = C14t_oPOM, colour = "oPOM"),
    size = 2
  ) +
  
  
  ## --- colour scale ---
  scale_colour_manual(
    name = "Pool",
    values = c(
      "Bulk" = col_bulk,
      "fPOM" = col_fPOM,
      "MAOM" = col_MAOM,
      "oPOM" = col_oPOM
    )
  ) +
  
  
  ## --- fill scale for uncertainty ribbons ---
  scale_fill_manual(
    values = c(
      "Bulk" = col_bulk,
      "fPOM" = col_fPOM,
      "MAOM" = col_MAOM,
      "oPOM" = col_oPOM
    ),
    guide = "none"
  ) +
  
  
  ## --- line types ---
  scale_linetype_manual(
    name = "Pool",
    values = c(
      "Bulk" = "solid",
      "fPOM" = "dotted",
      "MAOM" = "dotted",
      "oPOM" = "dotted"
    )
  ) +
  
  
  ## --- combine colour + linetype into one legend ---
  guides(
    colour = guide_legend(
      order = 1,
      override.aes = list(
        linetype = c("solid", "dotted", "dotted", "dotted"),
        linewidth = c(1, 0.8, 0.8, 0.8)
      )
    ),
    linetype = "none"
  ) +
  
  
  ## --- labels ---
  ggtitle(expression(Delta^14*C)) +
  ylab(Delta14Clabel) +
  xlab("Year") +
  
  
  ## --- theme ---
  theme_minimal()

# save results and plots
output_dir <- "plots/3p/control"

file_C <- file.path(output_dir, "C.png")
file_14C <- file.path(output_dir, "14C.png")

ggsave(filename = file_C, plot = plot_C,
       width = 7, height = 5, dpi = 300, bg = "white")
ggsave(filename = file_14C, plot = plot_14C,
       width = 7, height = 5, dpi = 300, bg = "white")

saveRDS(plot_C,
        file = file.path(output_dir, "C.rds"))
saveRDS(plot_14C,
        file = file.path(output_dir, "14C.rds"))

########################### edit from here
########## MCMC ############

var0 <- mFit$var_ms_unweighted
cov0 <- summary(mFit)$cov.scaled #  cov matrix can be used for jump 

# run 
#MCMC <- modMCMC(f=mc, p = bestpars, niter = 50000, jump = cov0*0.001, var0 = var0, wvar0 = 1, updatecov = 1000, burninlength =  1000, 
#                upper = c(1.5, 0.5 ,0.01, 0.5, 0.5),
#                lower = c(0.05, 0.005, 0.0001, 0, 0)) 
#save(MCMC, file = file.path("mod_runs", "MCMC_3p_Control.Rdata"))
load(here::here("mod_runs/MCMC_3p_Control.Rdata"))

# view distribution
class(MCMC) <- "modMCMC"
parsMCMC<-summary(MCMC)

# plot convergence
convergence_plot<-plot(MCMC) 

output_dir <- "plots/3p/control"
file_convergence_plot <- file.path(output_dir, "convergence_plot.png")

png(file_convergence_plot, width = 10, height = 8, units = "in", res = 300)
plot(MCMC)
dev.off()

# plot densitites
mcmc_long <- as.data.frame(MCMC$pars) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

density_plot <- ggplot(mcmc_long, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  facet_wrap(~ parameter, scales = "free") +
  theme_minimal() +
  labs(title = "MCMC Parameter Densities", x = "Value", y = "Density")

file_density_plot <- file.path(output_dir, "density_plot.png")
ggsave(filename = file_density_plot, plot = density_plot,
       width = 10, height = 8, dpi = 300, bg = "white")
saveRDS(density_plot,
        file = file.path(output_dir, "density_plot.rds"))


# check inter-dependence of parameters using 500 MCMC runs
pairs_plot<-pairs(MCMC,nsample=500) 

file_pairs_plot <- file.path(output_dir, "pairs_plot.png")

png(file_pairs_plot, width = 10, height = 8, units = "in", res = 300)
pairs(MCMC, nsample = 500)
dev.off()
saveRDS(pairs_plot,
        file = file.path(output_dir, "pairs_plot.rds"))

# performance
par(mfrow = c(1, 1))  
hist(MCMC$SS, breaks = 50)
percentage_accepted<-(100*MCMC$naccepted)/50000
MCMC_bestpars<-MCMC$bestpar
cost_modfit <- mc(bestpars)$model 
cost_mcmc   <- mc(MCMC_bestpars)$model  

# plots showing 95% credible interval of model outputs given the posterior parameter distribution
# i.e. posterior predictive uncertainty conditional on sampled MCMC chains

set.seed(1)
pars_sub <- MCMC$pars[sample(1:nrow(MCMC$pars), 500), ]

#runs_3p_control <- lapply(1:nrow(pars_sub), function(i) run_mod(pars_sub[i, ]))
#save(runs_3p_control, file = file.path("mod_runs", "runs_3p_control.Rdata"))
load(here::here("mod_runs/runs_3p_control.Rdata"))

# convert to array-like structure
extract_var <- function(var){
  sapply(runs_3p_control, function(x) x[[var]])
}

vars <- c("Ct","Ct_fPOM","Ct_MAOM","Ct_oPOM", 
          "C14t","C14t_fPOM","C14t_MAOM","C14t_oPOM")

unc_list <- lapply(vars, function(v){
  mat <- extract_var(v)
  
  data.frame(
    Year = runs_3p_control[[1]]$Year,
    var = v,
    Mean = rowMeans(mat, na.rm = TRUE),
    Low  = apply(mat, 1, quantile, 0.025, na.rm = TRUE),
    High = apply(mat, 1, quantile, 0.975, na.rm = TRUE)
  )
})

unc_df <- do.call(rbind, unc_list)

# var names for plotting
unc_df$Pool <- dplyr::case_when(
  unc_df$var == "Ct" ~ "Bulk",
  unc_df$var == "Ct_fPOM" ~ "fPOM",
  unc_df$var == "Ct_MAOM" ~ "MAOM",
  unc_df$var == "Ct_oPOM" ~ "oPOM",
  unc_df$var == "C14t" ~ "Bulk",
  unc_df$var == "C14t_fPOM" ~ "fPOM",
  unc_df$var == "C14t_MAOM" ~ "MAOM",
  unc_df$var == "C14t_oPOM" ~ "oPOM"
)
unc_C    <- subset(unc_df, grepl("^Ct", var))
unc_C14  <- subset(unc_df, grepl("^C14t", var))


############  C stocks plot with MCMC uncertainty ################
out_best_MCMC <- run_mod(MCMC$bestpar)

# plot using predictions made by mcmc

plot_C_final <- ggplot() +
  
  ## --- MCMC ribbons (model uncertainty) ---
  
  geom_ribbon(data = unc_C,
              aes(x = Year, ymin = Low, ymax = High, fill = Pool),
              alpha = 0.2) +
  
  ## --- model lines (All solid lines) ---
  
  geom_line(data = out_best_MCMC, aes(Year, Ct, colour = "Bulk"), linewidth = 1) +
  geom_line(data = out_best_MCMC, aes(Year, Ct_fPOM, colour = "fPOM"), linewidth = 0.8) +
  geom_line(data = out_best_MCMC, aes(Year, Ct_MAOM, colour = "MAOM"), linewidth = 0.8) +
  geom_line(data = out_best_MCMC, aes(Year, Ct_oPOM, colour = "oPOM"), linewidth = 0.8) +

  ## --- observations: points + error bars ---
  
  geom_point(data = Cobs_bulk, aes(Year, Ct, colour = "Bulk")) +
  geom_errorbar(data = Cobs_bulk,
                aes(Year, ymin = Ct - Ct_sd, ymax = Ct + Ct_sd, colour = "Bulk"),
                width = 0.5) +
  
  geom_point(data = Cobs_fPOM, aes(Year, Ct_fPOM, colour = "fPOM")) +
  geom_errorbar(data = Cobs_fPOM,
                aes(Year, ymin = Ct_fPOM - Ct_fPOM_sd, ymax = Ct_fPOM + Ct_fPOM_sd, colour = "fPOM"),
                width = 0.5) +
  
  geom_point(data = Cobs_MAOM, aes(Year, Ct_MAOM, colour = "MAOM")) +
  geom_errorbar(data = Cobs_MAOM,
                aes(Year, ymin = Ct_MAOM - Ct_MAOM_sd, ymax = Ct_MAOM + Ct_MAOM_sd, colour = "MAOM"),
                width = 0.5) +
  
  geom_point(data = Cobs_oPOM, aes(Year, Ct_oPOM, colour = "oPOM")) +
  geom_errorbar(data = Cobs_oPOM,
                aes(Year, ymin = Ct_oPOM - Ct_oPOM_sd, ymax = Ct_oPOM + Ct_oPOM_sd, colour = "oPOM"),
                width = 0.5) +
  
  ## --- scales ---
  
  scale_colour_manual(
    name = "Pool",
    values = c("Bulk" = col_bulk,
               "fPOM" = col_fPOM,
               "MAOM" = col_MAOM,
               "oPOM" = col_oPOM),
    breaks = c("Bulk", "fPOM", "MAOM", "oPOM"),
    labels = c("Bulk", "fPOM", "MAOM", "oPOM")
  ) +
  
  scale_fill_manual(
    name = "Pool",
    values = c("Bulk" = col_bulk,
               "fPOM" = col_fPOM,
               "MAOM" = col_MAOM,
               "oPOM" = col_oPOM),
    breaks = c("Bulk", "fPOM", "MAOM", "oPOM"),
    labels = c("Bulk", "fPOM", "MAOM", "oPOM")
  ) +
  guides(
    colour = guide_legend(order = 1),
    fill   = "none"
  ) +
  
  ylab(stocks_label) +
  xlab("Year") +
  
  ## ---publication-style theme ---
  
  theme_minimal() +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    axis.title = element_text(
      colour = "black",
      size = 11
    ),
    plot.title = element_text(
      colour = "black",
      size = 12,
      face = "bold"
    )
  )


############# 14C plot final ######################

plot_14C_final <- ggplot() +
  
  ## --- MCMC ribbons (model uncertainty) ---
  
  geom_ribbon(data = unc_C14,
              aes(x = Year, ymin = Low, ymax = High, fill = Pool),
              alpha = 0.2) +
  
  ## --- model lines ---
  
  geom_line(data = out_best_MCMC,
            aes(Year, C14t, colour = "Bulk"),
            linewidth = 1) +
  
  geom_line(data = out_best_MCMC,
            aes(Year, C14t_fPOM, colour = "fPOM"),
            linewidth = 0.8) +
  
  geom_line(data = out_best_MCMC,
            aes(Year, C14t_MAOM, colour = "MAOM"),
            linewidth = 0.8) +
  
  geom_line(data = out_best_MCMC,
            aes(Year, C14t_oPOM, colour = "oPOM"),
            linewidth = 0.8) +
  
  ## --- observations: points + error bars ---
  
  geom_point(data = C14obs_bulk,
             aes(Year, C14t, colour = "Bulk")) +
  
  geom_errorbar(data = C14obs_bulk,
                aes(Year,
                    ymin = C14t - C14t_sd,
                    ymax = C14t + C14t_sd,
                    colour = "Bulk"),
                width = 0.5) +
  
  geom_point(data = C14obs_fPOM,
             aes(Year, C14t_fPOM, colour = "fPOM")) +
  
  geom_errorbar(data = C14obs_fPOM,
                aes(Year,
                    ymin = C14t_fPOM - C14t_fPOM_sd,
                    ymax = C14t_fPOM + C14t_fPOM_sd,
                    colour = "fPOM"),
                width = 0.5) +
  
  geom_point(data = C14obs_MAOM,
             aes(Year, C14t_MAOM, colour = "MAOM")) +
  
  geom_errorbar(data = C14obs_MAOM,
                aes(Year,
                    ymin = C14t_MAOM - C14t_MAOM_sd,
                    ymax = C14t_MAOM + C14t_MAOM_sd,
                    colour = "MAOM"),
                width = 0.5) +
  
  geom_point(data = C14obs_oPOM,
             aes(Year, C14t_oPOM, colour = "oPOM")) +
  
  geom_errorbar(data = C14obs_oPOM,
                aes(Year,
                    ymin = C14t_oPOM - C14t_oPOM_sd,
                    ymax = C14t_oPOM + C14t_oPOM_sd,
                    colour = "oPOM"),
                width = 0.5) +
  
  ## --- scales ---
  
  scale_colour_manual(
    name = "Pool",
    values = c("Bulk" = col_bulk,
               "fPOM" = col_fPOM,
               "MAOM" = col_MAOM,
               "oPOM" = col_oPOM),
    breaks = c("Bulk", "fPOM", "MAOM", "oPOM"),
    labels = c("Bulk", "fPOM", "MAOM", "oPOM")
  ) +
  
  scale_fill_manual(
    name = "Pool",
    values = c("Bulk" = col_bulk,
               "fPOM" = col_fPOM,
               "MAOM" = col_MAOM,
               "oPOM" = col_oPOM),
    breaks = c("Bulk", "fPOM", "MAOM", "oPOM"),
    labels = c("Bulk", "fPOM", "MAOM", "oPOM")
  ) +
  
  ylab(Delta14Clabel) +
  xlab("Year") +
  
  ## --- publication-style theme ---
  
  theme_minimal() +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    axis.title = element_text(
      colour = "black",
      size = 11
    ),
    plot.title = element_text(
      colour = "black",
      size = 12,
      face = "bold"
    )
  )

# save individual plots
output_dir <- "plots/3p/control"
file_C_final <- file.path(output_dir, "C_final.png")
file_14C_final <- file.path(output_dir, "14C_final.png")

ggsave(filename = file_C_final, plot = plot_C_final, width = 7, height = 5, dpi = 300, bg = "white")
ggsave(filename = file_14C_final, plot = plot_14C_final, width = 7, height = 5, dpi = 300, bg = "white")

saveRDS(plot_C_final, file = file.path(output_dir, "C_final.rds"))
saveRDS(plot_14C_final, file = file.path(output_dir, "14C_final.rds"))

# combined MCMC plots

plot_C_notitle <- plot_C_final + ggtitle(NULL)
plot_14C_notitle <- plot_14C_final + ggtitle(NULL)

combined_plot <- (plot_C_notitle / plot_14C_notitle) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

combined_plot <- combined_plot +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = file.path(output_dir, "C_and_14C_combined.png"),
  plot = combined_plot,
  width = 7,
  height = 10,
  dpi = 300,
  bg = "white"
)

saveRDS(
  combined_plot,
  file = file.path(output_dir, "C_and_14C_combined.rds")
)



############# Average SD and prediction error over years ################

# Function to calculate mean uncertainty
calc_uncertainty <- function(obs_sd, pred_df, pool_name, variable){
  
  pred <- pred_df %>%
    filter(Pool == pool_name)
  
  data.frame(
    Variable = variable,
    Pool = pool_name,
    Observed_SD = ifelse(length(obs_sd) == 0, NA, mean(obs_sd, na.rm = TRUE)),
    Prediction_error_95CI = mean(pred$High - pred$Low, na.rm = TRUE),
    Prediction_error_95CI_half = mean((pred$High - pred$Low)/2, na.rm = TRUE)
  )
}

########## C stocks ##########

C_unc_summary <- bind_rows(
  
  calc_uncertainty(
    Cobs_bulk$Ct_sd,
    unc_C,
    "Bulk",
    "C stocks"
  ),
  
  calc_uncertainty(
    Cobs_fPOM$Ct_fPOM_sd,
    unc_C,
    "FPOM",
    "C stocks"
  ),
  
  calc_uncertainty(
    Cobs_oPOM$Ct_oPOM_sd,
    unc_C,
    "OPOM",
    "C stocks"
  ),
  
  calc_uncertainty(
    Cobs_MAOM$Ct_MAOM_sd,
    unc_C,
    "MAOM",
    "C stocks"
  )
)

########## 14C ##########

C14_unc_summary <- bind_rows(
  
  calc_uncertainty(
    C14obs_bulk$C14t_sd,
    unc_C14,
    "Bulk",
    "Delta14C"
  ),
  
  calc_uncertainty(
    NULL,
    unc_C14,
    "FPOM",
    "Delta14C"
  ),
  
  calc_uncertainty(
    NULL,
    unc_C14,
    "OPOM",
    "Delta14C"
  ),
  
  calc_uncertainty(
    C14obs_MAOM$C14t_MAOM_sd,
    unc_C14,
    "MAOM",
    "Delta14C"
  )
)
########## Combine results ##########

uncertainty_summary <- bind_rows(
  C_unc_summary,
  C14_unc_summary
)

# round for reporting
uncertainty_summary <- uncertainty_summary %>%
  mutate(
    across(
      c(Observed_SD,
        Prediction_error_95CI,
        Prediction_error_95CI_half),
      ~round(.x, 2)
    )
  )

print(uncertainty_summary)

########################################################################

################# Ages and transit times ###############

# sample mcmc pars
set.seed(1)

burnin <- 20000  # or more, see below
pars_full <- as.data.frame(MCMC$pars)
pars_post <- pars_full[-(1:burnin), ]
pars_sub <- pars_post[sample(1:nrow(pars_post), 700), ]

# sample ~700 parameter sets
n_samp <- 700
pars_sub <- pars_full[sample(1:nrow(pars_full), n_samp), ]

# func to build 3 pool system matrix
build_A_u <- function(pars){
  # kf, ki, ks, alpha_fi, alpha_is
  kf <- pars[1]
  ki <- pars[2]
  ks <- pars[3]
  alpha_fi <- pars[4]
  alpha_is <- pars[5]
  A <- diag(-c(kf, ki, ks)) # 3 pool system ---
  A[2,1] <- kf * alpha_fi
  A[3,2] <- ki * alpha_is
  u <- matrix(c(0.8*mean_C_inputs, 0.2*mean_C_inputs, 0), ncol = 1) # input vector 
  return(list(A = A, u = u))
}


# compute densities for one par set
get_age_tt_dens <- function(pars, ages){
  AU <- build_A_u(pars)
  SA <- systemAge(A = AU$A, u = AU$u, a = ages)
  TT <- transitTime(A = AU$A, u = AU$u, a = ages)
  data.frame(
    age = ages,
    system_age = SA$systemAgeDensity,
    fPOM = SA$poolAgeDensity[,1],
    oPOM = SA$poolAgeDensity[,2],
    MAOM = SA$poolAgeDensity[,3],
    transit_time = TT$transitTimeDensity
  )
}


# run for pars 
MCMC_bestpars<-MCMC$bestpar
ages <- seq(0, 500, by = 1)
dens_best <- get_age_tt_dens(MCMC_bestpars, ages)

# run for sampled mcmc pars 
#dens_list_3p_control <- lapply(1:nrow(pars_sub), function(i){
#  get_age_tt_dens(as.numeric(pars_sub[i, ]), ages)
#})

#save(dens_list_3p_control, file = file.path("mod_runs", "dens_list_3p_control.Rdata"))
load(here::here("mod_runs/dens_list_3p_control.Rdata"))

# uncertainty envelopes
extract_var <- function(var){
  sapply(dens_list_3p_control, function(x) x[[var]])
}


vars <- c("system_age","fPOM","oPOM","MAOM","transit_time")

unc_dens <- lapply(vars, function(v){
  
  mat <- extract_var(v)
  
  data.frame(
    age = ages,
    variable = v,
    mean = rowMeans(mat, na.rm = TRUE),
    median = apply(mat, 1, median, na.rm = TRUE),
    low  = apply(mat, 1, quantile, 0.025, na.rm = TRUE),
    high = apply(mat, 1, quantile, 0.975, na.rm = TRUE)
  )
})

unc_dens_df <- do.call(rbind, unc_dens)

# best pars for plotting
dens_best_long <- dens_best %>%
  pivot_longer(-age, names_to = "variable", values_to = "value")

# summary stats func
sum_fun <- function(age, dens){
  dens <- dens / sum(dens)   # <-- THIS is the key
  mean_age <- sum(age * dens)
  cdf <- cumsum(dens)
  
  data.frame(
    mean   = mean_age,
    median = age[which.min(abs(cdf - 0.5))]
  )
}

# summary stats from mean of MCMC densities
summary_list <- lapply(unique(unc_dens_df$variable), function(v){
  df <- unc_dens_df %>% filter(variable == v)
  stats <- sum_fun(df$age, df$mean)
  cbind(variable = v, stats)
})

summary_table <- do.call(rbind, summary_list)

summary_table_age_TT_3p_Control <- summary_table %>%
  mutate(across(-variable, ~round(., 1)))

save(summary_table_age_TT_3p_Control, file = file.path("mod_runs", "summary_table_age_TT_3p_Control.Rdata"))


# Prepare text labels for the mean & median values per facet
summary_labels <- summary_table_fmt %>%
  mutate(
    label_text = paste0("Mean: ", mean, "\nMedian: ", median)
  )

# Prepare panel tags data frame (A, B, C, D, E)
tag_df <- data.frame(
  variable = unique(unc_dens_df$variable),
  tag = LETTERS[1:length(unique(unc_dens_df$variable))]
)


# ==========================================
# 1. Plot of age and transit time
# ==========================================
plot_log_age_tt <- ggplot() +
  # --- ribbon ---
  geom_ribbon(data = unc_dens_df,
              aes(x = age, ymin = low, ymax = high),
              fill = "grey70",
              alpha = 0.5) +
  # --- mean density curve ---
  geom_line(data = unc_dens_df,
            aes(x = age, y = mean),
            colour = "black",
            linewidth = 0.8) +
  facet_wrap(~variable, scales = "free", ncol = 2) +
  
  # --- MEDIAN ---
  geom_vline(data = summary_table,
             aes(xintercept = median, colour = "Median"),
             linetype = "dotted",
             linewidth = 1) +
  # --- MEAN ---
  geom_vline(data = summary_table,
             aes(xintercept = mean, colour = "Mean"),
             linetype = "dashed",
             linewidth = 1) +
  
  # --- PANEL LABELS (A, B, C...) ---
  geom_text(data = tag_df,
            aes(x = -Inf, y = Inf, label = tag),
            hjust = -0.5, vjust = 1.3, size = 5, fontface = "bold",
            inherit.aes = FALSE) +
  
  # --- DIRECT TEXT ANNOTATION ---
  geom_text(data = summary_labels,
            aes(x = Inf, y = Inf, label = label_text),
            hjust = 1.1, vjust = 1.2, size = 4, fontface = "italic") +
  
  # --- legend control ---
  scale_colour_manual(
    name = "Statistic",
    values = c("Median" = "red", "Mean" = "blue")
  ) +
  scale_x_log10(
    limits = c(1, 500),
    breaks = c(1, 5, 10, 100, 500)
  ) +
  labs(x = "Age (years), log-scale", y = "Density") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# save plot
output_dir <- "plots/3p/control"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

file_plot_log_age_tt <- file.path(output_dir, "plot_log_age_tt.png")

png(file_plot_log_age_tt,
    width = 10,
    height = 8,
    units = "in",
    res = 300)

print(plot_log_age_tt)

dev.off()

################ --- Calculate fluxes and stocks for all pools ######################

flux_list <- lapply(1:length(runs_3p_control), function(i){
  
  pars <- as.numeric(MCMC$pars[i, ])
  run  <- runs_3p_control[[i]]
  
  # take pool sizes at final year (quasi steady state)
  C_fPOM  <- tail(run$Ct_fPOM, 1)
  C_oPOM <- tail(run$Ct_oPOM, 1)
  C_MAOM  <- tail(run$Ct_MAOM, 1)
  
  kf <- pars[1]
  ki <- pars[2]
  ks <- pars[3]
  alpha_fi <- pars[4]
  alpha_is <- pars[5]

  data.frame(
    
    # Gross decomposition
    fPOM_decomp  = kf * C_fPOM,
    oPOM_decomp = ki * C_oPOM,
    MAOM_decomp  = ks * C_MAOM,
    
    # Transfers
    fPOM_to_oPOM = alpha_fi * kf * C_fPOM,
    oPOM_to_MAOM = alpha_is * ki * C_oPOM,
    
    # Respiration (remainder)
    fPOM_resp  = (1 - alpha_fi) * kf * C_fPOM,
    oPOM_resp = (1 - alpha_is) * ki * C_oPOM,
    MAOM_resp  = ks * C_MAOM
  )
})

flux_df <- do.call(rbind, flux_list)

flux_summary <- flux_df %>%
  summarise(across(everything(),
                   list(
                     median = ~median(.),
                     lci = ~quantile(., 0.05),
                     uci = ~quantile(., 0.95)
                   )
  )) %>%
  pivot_longer(everything(),
               names_to = "Metric",
               values_to = "Value")

#visualize

############################################
# Prepare flux summary for carbon fate plot
############################################

############################################
# Prepare carbon fate data for plotting
############################################

flux_plot_df <- flux_df %>%
  
  # Convert the individual flux columns into rows
  pivot_longer(
    cols = everything(),
    names_to = "Flux",
    values_to = "Value"
  ) %>%
  
  # Calculate posterior summaries for each flux
  group_by(Flux) %>%
  summarise(
    median = median(Value, na.rm = TRUE),
    lci = quantile(Value, 0.05, na.rm = TRUE),
    uci = quantile(Value, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # Assign carbon fate
  mutate(
    Fate = case_when(
      grepl("_resp$", Flux) ~ "Respiration",
      grepl("_to_", Flux) ~ "Stabilization",
      TRUE ~ NA_character_
    ),
    
    # Assign source pool
    Pool = case_when(
      grepl("^fPOM", Flux) ~ "fPOM",
      grepl("^oPOM", Flux) ~ "oPOM",
      grepl("^MAOM", Flux) ~ "MAOM",
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(!is.na(Fate), !is.na(Pool)) %>%
  
  mutate(
    Fate = factor(
      Fate,
      levels = c("Respiration", "Stabilization")
    ),
    Pool = factor(
      Pool,
      levels = c("fPOM", "oPOM", "MAOM")
    )
  )
############################################
# Assign flux fate and source SOM pool
############################################

flux_plot_df <- flux_plot_df %>%
  mutate(
    
    # Carbon fate
    Fate = case_when(
      
      grepl("_resp$", Flux) ~ "Respiration",
      
      grepl("_to_", Flux) ~ "Stabilization",
      
      TRUE ~ NA_character_
    ),
    
    
    # Source pool
    Pool = case_when(
      
      grepl("^fPOM", Flux) ~ "fPOM",
      
      grepl("^oPOM", Flux) ~ "oPOM",
      
      grepl("^MAOM", Flux) ~ "MAOM",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    !is.na(Fate),
    !is.na(Pool)
  )



############################################
# Set factor order
############################################
flux_plot_df <- flux_plot_df %>%
  mutate(
    Fate = case_when(
      grepl("_resp$", Flux) ~ "Respiration",
      grepl("_to_", Flux) ~ "Stabilization",
      TRUE ~ NA_character_
    ),
    
    Pool = case_when(
      grepl("^fPOM", Flux) ~ "fPOM",
      grepl("^oPOM", Flux) ~ "oPOM",
      grepl("^MAOM", Flux) ~ "MAOM",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Fate), !is.na(Pool))

############################################
# Plot: carbon fate by source SOM pool
############################################

fluxplot<-ggplot(
  flux_plot_df,
  aes(
    x = Fate,
    y = median,
    fill = Pool
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7
  ) +
  
  geom_errorbar(
    aes(
      ymin = lci,
      ymax = uci
    ),
    position = position_dodge(
      width = 0.8
    ),
    width = 0.2,
    linewidth = 0.8
  ) +
  
  scale_fill_manual(
    values = c(
      "fPOM" = col_fPOM,
      "oPOM" = col_oPOM,
      "MAOM" = col_MAOM
    )
  ) +
  
  labs(
    x = NULL,
    y = expression(
      "Carbon flux (g C m"^-2*" yr"^-1*")"
    ),
    fill = "Source SOM pool"
  ) +
  
  theme_minimal(
    base_size = 14
  ) +
  
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )


fluxplot_path <- file.path(output_dir, "flux_plot.png")

png(fluxplot_path, width = 10, height = 8, units = "in", res = 300)
print(fluxplot)
dev.off()


############## efficiency
## ============================================================
## Calculate posterior efficiencies
## ============================================================

efficiency_list <- lapply(1:length(runs_3p_control), function(i){
  
  pars <- as.numeric(MCMC$pars[i, ])
  run  <- runs_3p_control[[i]]
  
  # final pool sizes (quasi steady state)
  C_fPOM  <- tail(run$Ct_fPOM, 1)
  C_oPOM <- tail(run$Ct_oPOM, 1)
  C_MAOM  <- tail(run$Ct_MAOM, 1)
  
  # decomposition rates
  kf <- pars[1]
  ki <- pars[2]
  ks <- pars[3]
  
  # transfer coefficients
  alpha_fi <- pars[4]   # fPOM -> oPOMmediate
  alpha_is <- pars[5]   # oPOMmediate -> MAOM
  
  
  ## -------------------------------
  ## FPOM POOL
  ## -------------------------------
  
  fPOM_decomp <- kf * C_fPOM
  
  fPOM_stabilization <- alpha_fi * kf * C_fPOM
  

  fPOM_resp <- fPOM_decomp -
    fPOM_stabilization 

  
  ## -------------------------------
  ## OPOMMEDIATE POOL
  ## -------------------------------
  
  oPOM_decomp <- ki * C_oPOM
  
  oPOM_stabilization <- alpha_is * ki * C_oPOM
  

  oPOM_resp <- oPOM_decomp -
    oPOM_stabilization 
  
  
  ## -------------------------------
  ## MAOM POOL
  ## -------------------------------
  
  MAOM_decomp <- ks * C_MAOM
  
  MAOM_resp <- MAOM_decomp
  
  
  data.frame(
    
    fPOM_resp_eff =
      fPOM_resp / fPOM_decomp,
    
    fPOM_stab_eff =
      fPOM_stabilization / fPOM_decomp,
    
    
    oPOM_resp_eff =
      oPOM_resp / oPOM_decomp,
    
    oPOM_stab_eff =
      oPOM_stabilization / oPOM_decomp,
    
    
    MAOM_resp_eff = 1
    
  )
})


efficiency_df <- do.call(rbind, efficiency_list)


## ============================================================
## Summarize posterior distributions
## ============================================================

efficiency_summary <- efficiency_df %>%
  summarise(
    across(
      everything(),
      list(
        median = ~median(.),
        lci = ~quantile(., 0.05),
        uci = ~quantile(., 0.95)
      )
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = "Metric",
    values_to = "Value"
  )


efficiency_summary


# possibly not necessary
# ============================================================
# Weighted C and N age density distributions
# THREE-POOL SYSTEM: fPOM, oPOM, MAOM
# ============================================================


# ============================================================
# 1. Final stocks for best parameter set
# ============================================================

best_run <- out_best_MCMC
final_yr_idx <- nrow(best_run)


# Best-fit C pools
best_C_stocks <- c(
  fPOM = best_run$Ct_fPOM[final_yr_idx],
  oPOM = best_run$Ct_oPOM[final_yr_idx],
  MAOM = best_run$Ct_MAOM[final_yr_idx]
)


# observed final C and N stocks

final_Cobs <- c(
  fPOM = tail(Cobs_fPOM$Ct_fPOM, 1),
  oPOM = tail(Cobs_oPOM$Ct_oPOM, 1),
  MAOM = tail(Cobs_MAOM$Ct_MAOM, 1)
)

final_Nobs <- c(
  fPOM = tail(Nobs_fPOM$Nt_fPOM, 1),
  oPOM = tail(Nobs_oPOM$Nt_oPOM, 1),
  MAOM = tail(Nobs_MAOM$Nt_MAOM, 1)
)


# ============================================================
# Pool-specific C:N ratios
# ============================================================

CN_bulk_final <- tail(bulk$CN_mean, 1)

CN_fPOM_final <- final_Cobs["fPOM"] / final_Nobs["fPOM"]
CN_oPOM_final <- final_Cobs["oPOM"] / final_Nobs["oPOM"]
CN_MAOM_final <- final_Cobs["MAOM"] / final_Nobs["MAOM"]


# ============================================================
# 2. BEST-FIT density-weighted system
# ============================================================

dens_best_weighted <- dens_best %>%
  
  mutate(
    
    # density x Carbon stocks
    C_fPOM_stock = fPOM * best_C_stocks["fPOM"],
    C_oPOM_stock = oPOM * best_C_stocks["oPOM"],
    C_MAOM_stock = MAOM * best_C_stocks["MAOM"],
    
    C_system_stock =
      C_fPOM_stock +
      C_oPOM_stock +
      C_MAOM_stock
    
  ) %>%
  
  mutate(
    
    # density x Nitrogen stocks
    N_fPOM_stock =
      fPOM * best_C_stocks["fPOM"] * 1 / CN_fPOM_final,
    
    N_oPOM_stock =
      oPOM * best_C_stocks["oPOM"] * 1 / CN_oPOM_final,
    
    N_MAOM_stock =
      MAOM * best_C_stocks["MAOM"] * 1 / CN_MAOM_final,
    
    N_system_stock =
      N_fPOM_stock +
      N_oPOM_stock +
      N_MAOM_stock
    
  )


# ============================================================
# 3. MCMC ENSEMBLE
# ============================================================

stock_vars <- c(
  "C_fPOM_stock",
  "C_oPOM_stock",
  "C_MAOM_stock",
  "C_system_stock",
  
  "N_fPOM_stock",
  "N_oPOM_stock",
  "N_MAOM_stock",
  "N_system_stock"
)


mcmc_stock_list <- lapply(seq_along(runs_3p_control), function(i){
  
  r <- runs_3p_control[[i]]
  d <- dens_list_3p_control[[i]]
  f_idx <- nrow(r)
  
  
  # ----------------------------------------------------------
  # Carbon pools
  # ----------------------------------------------------------
  
  Cf <- d$fPOM * r$Ct_fPOM[f_idx]
  Co <- d$oPOM * r$Ct_oPOM[f_idx]
  Cm <- d$MAOM * r$Ct_MAOM[f_idx]
  
  C_sys <- Cf + Co + Cm
  
  
  # ----------------------------------------------------------
  # Nitrogen pools
  # using final observed C:N ratios
  # ----------------------------------------------------------
  
  Nf <- d$fPOM * r$Ct_fPOM[f_idx] * 1 / CN_fPOM_final
  No <- d$oPOM * r$Ct_oPOM[f_idx] * 1 / CN_oPOM_final
  Nm <- d$MAOM * r$Ct_MAOM[f_idx] * 1 / CN_MAOM_final
  
  N_sys <- Nf + No + Nm
  
  
  data.frame(
    
    age = d$age,
    
    # Carbon
    C_fPOM_stock = Cf,
    C_oPOM_stock = Co,
    C_MAOM_stock = Cm,
    C_system_stock = C_sys,
    
    # Nitrogen
    N_fPOM_stock = Nf,
    N_oPOM_stock = No,
    N_MAOM_stock = Nm,
    N_system_stock = N_sys
    
  )
})


# ============================================================
# 4. Uncertainty envelopes
# ============================================================

ages <- dens_list_3p_control[[1]]$age


unc_stocks_list <- lapply(stock_vars, function(v){
  
  mat <- sapply(
    mcmc_stock_list,
    function(x) x[[v]]
  )
  
  data.frame(
    age = ages,
    variable = v,
    low = apply(
      mat,
      1,
      quantile,
      0.025,
      na.rm = TRUE
    ),
    high = apply(
      mat,
      1,
      quantile,
      0.975,
      na.rm = TRUE
    )
  )
})


unc_stocks_df <- do.call(
  rbind,
  unc_stocks_list
)


unc_stocks_df <- unc_stocks_df %>%
  mutate(
    type = ifelse(
      grepl("^C_", variable),
      "C",
      "N"
    )
  )


# ============================================================
# 5. Merge best-fit + uncertainty
# ============================================================

plot_df_final <- dens_best_weighted %>%
  
  select(
    age,
    contains("_stock")
  ) %>%
  
  pivot_longer(
    -age,
    names_to = "variable",
    values_to = "best_val"
  ) %>%
  
  mutate(
    type = ifelse(
      grepl("^C_", variable),
      "C",
      "N"
    )
  ) %>%
  
  left_join(
    unc_stocks_df,
    by = c(
      "age",
      "variable",
      "type"
    )
  )


plot_df_final$variable <- recode(
  plot_df_final$variable,
  
  "C_system_stock" = "C system (reconstructed)",
  
  "N_system_stock" = "N system (reconstructed)",
  
  "C_fPOM_stock" = "C fPOM",
  "C_oPOM_stock" = "C oPOM",
  "C_MAOM_stock" = "C MAOM",
  
  "N_fPOM_stock" = "N fPOM",
  "N_oPOM_stock" = "N oPOM",
  "N_MAOM_stock" = "N MAOM"
)


# ============================================================
# 6. Summary statistics
# ============================================================

fPOM_ref <- summary_table_fmt$median[
  summary_table_fmt$variable == "fPOM"
]

MAOM_ref <- summary_table_fmt$median[
  summary_table_fmt$variable == "MAOM"
]


final_stats <- dens_best_weighted %>%
  
  select(
    age,
    contains("_stock")
  ) %>%
  
  pivot_longer(
    -age,
    names_to = "variable",
    values_to = "value"
  ) %>%
  
  group_by(variable) %>%
  
  summarise(
    
    total_stock = sum(value),
    
    mean_val =
      sum(age * value) /
      sum(value),
    
    median_val =
      age[
        which.min(
          abs(
            cumsum(value) /
              sum(value) -
              0.5
          )
        )
      ],
    
    frac_younger_fPOM =
      sum(
        value[age <= fPOM_ref]
      ) /
      sum(value),
    
    frac_older_MAOM =
      sum(
        value[age >= MAOM_ref]
      ) /
      sum(value)
    
  )


final_stats$variable <- recode(
  final_stats$variable,
  
  "C_system_stock" = "C system (reconstructed)",
  "N_system_stock" = "N system (reconstructed)",
  
  "C_fPOM_stock" = "C fPOM",
  "C_oPOM_stock" = "C oPOM",
  "C_MAOM_stock" = "C MAOM",
  
  "N_fPOM_stock" = "N fPOM",
  "N_oPOM_stock" = "N oPOM",
  "N_MAOM_stock" = "N MAOM"
)


print(final_stats)


save(
  final_stats,
  file = file.path(
    "mod_runs",
    "final_stats_weighted_3p_control.Rdata"
  )
)


# ============================================================
# WEIGHTED AGE DISTRIBUTION PLOTS
# ============================================================


# ============================================================
# Unified colour system
# ============================================================

col_fPOM <- "#1b9e77"
col_MAOM <- "#BF40BF"
col_oPOM <- "#0000FF"


pool_colors <- c(
  
  "C fPOM" = col_fPOM,
  "C oPOM" = col_oPOM,
  "C MAOM" = col_MAOM,
  
  "N fPOM" = col_fPOM,
  "N oPOM" = col_oPOM,
  "N MAOM" = col_MAOM,
  
  "C system (reconstructed)" = "black",
  "N system (reconstructed)" = "black",
  
  "Mean" = "blue",
  "Median" = "red"
)


# ============================================================
# FUNCTION
# ============================================================

create_stock_age_plot <- function(
    data_subset,
    y_label,
    text_x_coords = NULL,
    var_labels = NULL
) {
  
  stats_subset <- final_stats %>%
    filter(
      variable %in%
        unique(data_subset$variable)
    )
  
  
  # Assign custom x position per panel text
  if (!is.null(text_x_coords)) {
    
    stats_subset$text_x <-
      text_x_coords[
        as.character(
          stats_subset$variable
        )
      ]
    
  } else if (
    !"text_x" %in%
    colnames(stats_subset)
  ) {
    
    stats_subset$text_x <- 10
    
  }
  
  
  ggplot(
    data_subset,
    aes(x = age)
  ) +
    
    # ribbons
    geom_ribbon(
      aes(
        ymin = low,
        ymax = high,
        fill = variable
      ),
      alpha = 0.4,
      show.legend = FALSE
    ) +
    
    # best-fit lines
    geom_line(
      aes(
        y = best_val,
        colour = variable
      ),
      linewidth = 0.8
    ) +
    
    # mean
    geom_vline(
      data = stats_subset,
      aes(
        xintercept = mean_val,
        colour = "Mean"
      ),
      linetype = "dotted",
      linewidth = 1.0
    ) +
    
    # median
    geom_vline(
      data = stats_subset,
      aes(
        xintercept = median_val,
        colour = "Median"
      ),
      linetype = "dotted",
      linewidth = 1.0
    ) +
    
    # mean / median text
    geom_text(
      data = stats_subset,
      aes(
        x = text_x,
        y = Inf,
        label = paste0(
          "Mean = ",
          round(mean_val, 0),
          " y\n",
          "Median = ",
          round(median_val, 0),
          " y"
        )
      ),
      hjust = 0,
      vjust = 2.5,
      size = 4.2,
      colour = "black",
      inherit.aes = FALSE
    ) +
    
    facet_wrap(
      ~variable,
      scales = "free_y",
      ncol = 2,
      labeller =
        if (!is.null(var_labels))
          labeller(
            variable = var_labels
          )
      else
        "label_value"
    ) +
    
    scale_x_log10(
      limits = c(1, 500),
      breaks = c(
        1,
        5,
        10,
        100,
        500
      )
    ) +
    
    scale_fill_manual(
      values = pool_colors,
      guide = "none"
    ) +
    
    scale_colour_manual(
      values = pool_colors,
      breaks = c(
        "Mean",
        "Median"
      ),
      name = "Statistics"
    ) +
    
    labs(
      x = "Age (years, log-scale)",
      y = y_label
    ) +
    
    theme_minimal(
      base_size = 14
    ) +
    
    theme(
      
      legend.position = "bottom",
      
      plot.title =
        element_blank(),
      
      panel.border =
        element_rect(
          colour = "black",
          fill = NA,
          linewidth = 0.5
        ),
      
      axis.ticks =
        element_line(
          colour = "black",
          linewidth = 0.6
        ),
      
      axis.ticks.length =
        unit(3, "pt"),
      
      axis.text =
        element_text(
          colour = "black"
        ),
      
      axis.title =
        element_text(
          colour = "black"
        )
    )
}


# ============================================================
# Custom text positions
# ============================================================

text_x_C <- c(
  
  "C fPOM" = 80,
  
  "C oPOM" = 1,
  
  "C MAOM" = 1,
  
  "C system (reconstructed)" = 1
)


text_x_N <- c(
  
  "N fPOM" = 80,
  
  "N oPOM" = 1,
  
  "N MAOM" = 1,
  
  "N system (reconstructed)" = 1
)


# ============================================================
# Clean titles
# ============================================================

var_labels_C <- c(
  
  "C fPOM" =
    "(A) Carbon fPOM",
  
  "C oPOM" =
    "(B) Carbon oPOM",
  
  "C MAOM" =
    "(C) Carbon MAOM",
  
  "C system (reconstructed)" =
    "(D) Carbon system"
)


var_labels_N <- c(
  
  "N fPOM" =
    "(E) Nitrogen fPOM",
  
  "N oPOM" =
    "(F) Nitrogen oPOM",
  
  "N MAOM" =
    "(G) Nitrogen MAOM",
  
  "N system (reconstructed)" =
    "(H) Nitrogen system"
)


# ============================================================
# Create Carbon Plot
# ============================================================

plot_C_weighted_final <- create_stock_age_plot(
  
  data_subset =
    subset(
      plot_df_final,
      type == "C"
    ),
  
  y_label =
    expression(
      Carbon~stocks~(g~m^{-2})
    ),
  
  text_x_coords =
    text_x_C,
  
  var_labels =
    var_labels_C
)


# ============================================================
# Create Nitrogen Plot
# ============================================================

plot_N_weighted_final <- create_stock_age_plot(
  
  data_subset =
    subset(
      plot_df_final,
      type == "N"
    ),
  
  y_label =
    expression(
      Nitrogen~stocks~(g~m^{-2})
    ),
  
  text_x_coords =
    text_x_N,
  
  var_labels =
    var_labels_N
)


# ============================================================
# Combined Plot
# ============================================================

combined_CN_plot <-
  (plot_C_weighted_final /
     plot_N_weighted_final) +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )


# ============================================================
# SAVE
# ============================================================

output_dir <- "plots/3p/control"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


ggsave(
  file.path(
    output_dir,
    "plot_C_weighted_final.png"
  ),
  plot_C_weighted_final,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)


ggsave(
  file.path(
    output_dir,
    "plot_N_weighted_final.png"
  ),
  plot_N_weighted_final,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)


ggsave(
  file.path(
    output_dir,
    "CN_age_distribution_combined.png"
  ),
  combined_CN_plot,
  width = 9,
  height = 11,
  dpi = 300,
  bg = "white"
)


saveRDS(
  plot_C_weighted_final,
  file.path(
    output_dir,
    "plot_C_weighted_final.rds"
  )
)


saveRDS(
  plot_N_weighted_final,
  file.path(
    output_dir,
    "plot_N_weighted_final.rds"
  )
)


saveRDS(
  combined_CN_plot,
  file.path(
    output_dir,
    "CN_age_distribution_combined.rds"
  )
)


# ============================================================
# CREATE C:N AGE DISTRIBUTIONS
# ============================================================

cn_pool_list <- lapply(
  seq_along(mcmc_stock_list),
  function(i){
    
    x <- mcmc_stock_list[[i]]
    
    data.frame(
      
      age = x$age,
      
      CN_fPOM =
        x$C_fPOM_stock /
        x$N_fPOM_stock,
      
      CN_oPOM =
        x$C_oPOM_stock /
        x$N_oPOM_stock,
      
      CN_MAOM =
        x$C_MAOM_stock /
        x$N_MAOM_stock,
      
      CN_system =
        x$C_system_stock /
        x$N_system_stock
    )
  }
)


# ============================================================
# SYSTEM C:N
# ============================================================

ages <- cn_pool_list[[1]]$age

cn_sys_mat <-
  sapply(
    cn_pool_list,
    function(x)
      x$CN_system
  )


cn_sys_df <- data.frame(
  
  age = ages,
  
  mean =
    rowMeans(
      cn_sys_mat,
      na.rm = TRUE
    ),
  
  low =
    apply(
      cn_sys_mat,
      1,
      quantile,
      0.025,
      na.rm = TRUE
    ),
  
  high =
    apply(
      cn_sys_mat,
      1,
      quantile,
      0.975,
      na.rm = TRUE
    ),
  
  median =
    apply(
      cn_sys_mat,
      1,
      median,
      na.rm = TRUE
    )
)


# Extract system stats
cn_stats_subset <- final_stats %>%
  filter(
    variable ==
      "C system (reconstructed)"
  )


# Calculate CN statistics for annotation
cn_stats_labels <- cn_stats_subset %>%
  mutate(
    
    label_text = paste0(
      
      "Mean = ",
      round(mean_val, 1),
      " y\n",
      
      "\nMedian = ",
      round(median_val, 1),
      " y"
    )
  )


# ============================================================
# LOG AGE VERSION
# ============================================================

plot_CN_system <- ggplot(
  cn_sys_df,
  aes(x = age)
) +
  
  # uncertainty ribbon
  geom_ribbon(
    aes(
      ymin = low,
      ymax = high
    ),
    fill = "black",
    alpha = 0.2
  ) +
  
  # mean CN trajectory
  geom_line(
    aes(y = mean),
    colour = "black",
    linewidth = 0.9
  ) +
  
  # mean
  geom_vline(
    data = cn_stats_subset,
    aes(
      xintercept = mean_val,
      colour = "Mean"
    ),
    linetype = "dotted",
    linewidth = 1.2
  ) +
  
  # median
  geom_vline(
    data = cn_stats_subset,
    aes(
      xintercept = median_val,
      colour = "Median"
    ),
    linetype = "dotted",
    linewidth = 1.2
  ) +
  
  # annotation
  geom_text(
    data = cn_stats_labels,
    aes(
      x = mean_val,
      y = Inf,
      label = label_text
    ),
    hjust = -0.05,
    vjust = 1.3,
    size = 5,
    fontface = "italic",
    inherit.aes = FALSE
  ) +
  
  scale_x_log10(
    limits = c(1, 500),
    breaks = c(
      1,
      5,
      10,
      100,
      500
    )
  ) +
  
  scale_colour_manual(
    values = c(
      "Mean" = "blue",
      "Median" = "red"
    ),
    breaks = c(
      "Mean",
      "Median"
    ),
    name = "Statistics"
  ) +
  
  labs(
    x = "Age (years, log-scale)",
    y = "C:N ratio (system)"
  ) +
  
  theme_minimal(
    base_size = 14
  ) +
  
  theme(
    
    panel.border =
      element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),
    
    axis.line =
      element_line(
        colour = "black",
        linewidth = 0.4
      ),
    
    axis.ticks =
      element_line(
        colour = "black",
        linewidth = 0.4
      ),
    
    axis.text =
      element_text(
        colour = "black",
        size = 12
      ),
    
    axis.title =
      element_text(
        colour = "black",
        size = 14
      ),
    
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "bottom"
  )


# display
plot_CN_system



# ============================================================
# SAVE C:N PLOTS
# ============================================================

ggsave(
  file.path(
    output_dir,
    "plot_CN_system.png"
  ),
  plot_CN_system,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)


saveRDS(
  plot_CN_system,
  file.path(
    output_dir,
    "plot_CN_system.rds"
  )
)


##########################
#compare CN

# ------------------------------------------------------
# Final year observations
# ------------------------------------------------------

final_year <- 2026

final_dat <- all_data %>%
  filter(Year == final_year)

# create C and N stocks
final_pools <- final_dat %>%
  select(LTE,
         Temperature...C.,
         C_stocks_gm2,
         N_stocks_gm2) %>%
  pivot_wider(
    names_from = Temperature...C.,
    values_from = c(C_stocks_gm2, N_stocks_gm2)
  ) %>%
  mutate(
    
    C_fPOM  = C_stocks_gm2_Soil - C_stocks_gm2_325,
    C_oPOM = C_stocks_gm2_325  - C_stocks_gm2_400,
    C_MAOM  = C_stocks_gm2_400,
    
    N_fPOM  = N_stocks_gm2_Soil - N_stocks_gm2_325,
    N_oPOM = N_stocks_gm2_325  - N_stocks_gm2_400,
    N_MAOM  = N_stocks_gm2_400,
    
    CN_fPOM  = C_fPOM / N_fPOM,
    CN_oPOM = C_oPOM / N_oPOM,
    CN_MAOM  = C_MAOM / N_MAOM
  )

# ------------------------------------------------------
# Long format
# ------------------------------------------------------

CN_long <- final_pools %>%
  select(LTE, CN_fPOM, CN_oPOM, CN_MAOM) %>%
  pivot_longer(
    -LTE,
    names_to = "Pool",
    values_to = "CN"
  )

CN_long$Pool <- factor(
  CN_long$Pool,
  levels = c("CN_fPOM","CN_oPOM","CN_MAOM"),
  labels = c("FPOM","OPOM","MAOM")
)

# ------------------------------------------------------
# Repeated measures ANOVA
# ------------------------------------------------------

mod <- aov(
  CN ~ Pool + Error(LTE/Pool),
  data = CN_long
)

summary(mod)

# ------------------------------------------------------
# Pairwise comparisons
# ------------------------------------------------------

lm_mod <- lm(CN ~ LTE + Pool, data = CN_long)

pairs(
  emmeans(lm_mod, "Pool"),
  adjust = "tukey"
)

friedman.test(
  CN ~ Pool | LTE,
  data = CN_complete
)

pairwise.wilcox.test( #doesnt make sense to report because too few samples
  CN_long$CN,
  CN_long$Pool,
  paired = TRUE,
  p.adjust.method = "holm"
)

#check 
CN_long %>%
  arrange(LTE, Pool)

# report bulk CN vs system age relationship
cn_age_lm <- lm(
  mean ~ log10(age),
  data = cn_sys_df %>% filter(age > 0)
)

summary(cn_age_lm)

