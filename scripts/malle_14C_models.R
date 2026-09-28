# SOM decomposition models for Malle soil

## load libraries
library(here)
library(SoilR)
library(FME)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

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


## inputs
mean_C_inputs <- 638.2 # g m-2

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
radiocarbon_samples$C_stocks_gm2 <- radiocarbon_samples$Ctotal * radiocarbon_samples$bulk_density * 0.2 * 10000
radiocarbon_samples$N_stocks_gm2 <- radiocarbon_samples$Ntotal * radiocarbon_samples$bulk_density * 0.2 * 10000

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
    inputFluxes = c(mean_C_inputs,0,0),
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
mFit <- modFit(
  f = mc,
  p = inipars, 
  method = "Nelder-Mead",
  upper = c(0.5, 0.2 ,0.005, 0.5, 0.5),
  lower = c(0.05, 0.01, 0.0001, 0, 0) 
)

save(mFit, file = file.path("mod_runs", "mFit_3p.Rdata"))
load(here::here("mod_runs/mFit_3p.Rdata"))