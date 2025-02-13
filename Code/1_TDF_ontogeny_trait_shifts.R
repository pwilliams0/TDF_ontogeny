#======================================================#
# COMPARE TRAITS OF SEEDLINGS AND ADULTS, TRAIT SHIFTS #
#======================================================#
library(tidyverse)
library(nimble)
library(Matrix)
library(MASS)

# ========== LOAD FULL TRAIT DATASET ==========

# Load combined seedling and adult traits
# Includes species with incomplete trait data
# Log-transform skewed traits
traits_data <- read.csv("Data/traits_data.csv") %>%
    mutate(log_LA = log(LA_cm2),
           log_SLA = log(SLA_cm2g),
           log_dh = log(dh_um),
           log_VD = log(VD_nomm2),
           log_DApit = log(DApit_um),
           id = paste(code,site,stage,X))

# ========== ANALYZE SHIFT IN TRAIT VALUES ==========

# ----- Ontogenetic trait shifts ---------- 

# Calculate shift in trait values from seedling to adult
# Arguments:
#   Trait = name of trait (log_SLA, LDMC_gg, etc.)
trait_shift <- function(Trait){

  # Create data list
  data <- traits_data %>%
    dplyr::filter(get({{Trait}}) > 0) %>%
    mutate(Trait = scale(get(Trait))) %>%
    dplyr::select(stage, code, site, Trait,
                  leaf_phenology, life_form) %>%
    na.omit() %>%
    mutate(pop_id = paste(site, code, sep="_")) %>%
    group_by(pop_id) %>%
    arrange(stage) %>%
    mutate(change_trait = Trait[1] - Trait[2]) %>%
    slice(1) %>% ungroup %>%
    dplyr::select(-stage) %>%
    na.omit()
  data_list <- list(change_trait = data$change_trait,
                    n = nrow(data))
  str(data_list)
  
  # Write model in NIMBLE
  model <- nimbleCode({
    
    # Priors
    mu ~ dnorm(0,sd=100) # Mean change in trait value
    sigma ~ dlnorm(0, sd=100) # SD
    
    # Likelihood
    for (i in 1:n){
      change_trait[i] ~ dnorm(mu, sd = sigma)
    }
    
  })
  
  # Set initial values
  inits <- list(mu = rnorm(1), sigma = rlnorm(1))
  
  # Run model
  trait_shift_model <- nimbleMCMC(code = model, 
                                  constants = data_list,
                                  inits = inits,
                                  niter = 4000, nburnin = 1000,
                                  nchains = 3,
                                  samplesAsCodaMCMC = TRUE)
  
  print(paste("Sample size:",data_list$n))
  return(trait_shift_model)
  
}

# LA
output <- trait_shift("log_LA")
summary(output)
saveRDS(output, "Data/shift_LA.RDS")
# LDMC
output <- trait_shift("LDMC_gg")
summary(output)
saveRDS(output, "Data/shift_LDMC.RDS")
# SLA
output <- trait_shift("log_SLA")
summary(output)
saveRDS(output, "Data/shift_SLA.RDS")
# WSG
output <- trait_shift("WSG_gcm3")
summary(output)
saveRDS(output, "Data/shift_WSG.RDS")
# TFW
output <- trait_shift("TFW_um")
summary(output)
saveRDS(output, "Data/shift_TFW.RDS")
# dh
output <- trait_shift("log_dh")
summary(output)
saveRDS(output, "Data/shift_dh.RDS")
# VD
output <- trait_shift("log_VD")
summary(output)
saveRDS(output, "Data/shift_VD.RDS")
# DApit
output <- trait_shift("log_DApit")
summary(output)
saveRDS(output, "Data/shift_DApit.RDS")

# ----- Trait shifts by site, phenology, & growth form -----

# Calculate shift in trait values from seedling to adult
#   Compare shifts by site, phenology, & growth form
# Arguments:
#   Trait = name of trait (log_SLA, LDMC_gg, etc.)
compare_shift <- function(Trait){
  
  # Create data list
  data <- traits_data %>%
    dplyr::filter(get({{Trait}}) > 0) %>%
    mutate(Trait = scale(get(Trait))) %>%
    dplyr::select(stage, code, site, Trait,
                  leaf_phenology, life_form) %>%
    na.omit() %>%
    mutate(pop_id = paste(site, code, sep="_")) %>%
    group_by(pop_id) %>%
    arrange(stage) %>%
    mutate(change_trait = Trait[1] - Trait[2]) %>%
    slice(1) %>% ungroup %>%
    dplyr::select(-stage) %>%
    na.omit()
  change_trait <- data$change_trait
  decid <- as.numeric(data$leaf_phenology == "deciduous")
  liana <- as.numeric(data$life_form == "liana")
  site1 <- as.numeric(data$site == "Jabiru")
  site2 <- as.numeric(data$site == "Colorados")
  site3 <- as.numeric(data$site == "Cotove")
  site4 <- as.numeric(data$site == "Tayrona")
  data_list <- list(change_trait = change_trait,
                    site1 = site1, site2 = site2,
                    site3 = site3, site4 = site4,
                    decid = decid, liana = liana,
                    n = nrow(data))
  str(data_list)
  
  # Write model in NIMBLE
  model <- nimbleCode({
    
    # Priors
    a1_Jabi ~ dnorm(0,sd=100) # Intercept for Jabiru
    a2_Colo ~ dnorm(0,sd=100) # Intercept for Colorados
    a3_Coto ~ dnorm(0,sd=100) # Intercept for Cotove
    a4_Tayr ~ dnorm(0,sd=100) # Intercept for Tayrona
    b1_decid ~ dnorm(0,sd=100) # Effect of deciduous
    b2_liana ~ dnorm(0,sd=100) # Effect of liana
    sigma ~ dlnorm(0, sd=100) # SD
    
    # Likelihood
    for (i in 1:n){
      change_trait[i] ~ dnorm(mu[i], sd = sigma)
      mu[i] <- a1_Jabi * site1[i] + a2_Colo * site2[i] +
        a3_Coto * site3[i] + a4_Tayr * site4[i] +
        b1_decid * decid[i] + b2_liana * liana[i]
    }
    
    # Derived quantities
    diff_Jabi_Colo <- a1_Jabi - a2_Colo
    diff_Jabi_Coto <- a1_Jabi - a3_Coto
    diff_Jabi_Tayr <- a1_Jabi - a4_Tayr
    diff_Colo_Coto <- a2_Colo - a3_Coto
    diff_Colo_Tayr <- a2_Colo - a4_Tayr
    diff_Coto_Tayr <- a3_Coto - a4_Tayr
    
  })
  
  # Set initial values
  inits <- list(a1_Jabi = rnorm(1), a2_Colo = rnorm(1),
                a3_Coto = rnorm(1), a4_Tayr = rnorm(1),
                b1_decid = rnorm(1), b2_liana = rnorm(1),
                sigma = rlnorm(1))
  
  # Parameters to estimate
  params <- c("a1_Jabi", "a2_Colo", "a3_Coto", "a4_Tayr",
              "b1_decid", "b2_liana", "sigma",
              "diff_Jabi_Colo", "diff_Jabi_Coto",
              "diff_Jabi_Tayr", "diff_Colo_Coto",
              "diff_Colo_Tayr", "diff_Coto_Tayr")
  
  # Run model
  compare_shift_model <- nimbleMCMC(code = model, 
                                    constants = data_list,
                                    inits = inits,
                                    monitors = params,
                                    niter = 4000, nburnin = 1000,
                                    nchains = 3,
                                    samplesAsCodaMCMC = TRUE)
  
  print(paste("Sample size:",data_list$n))
  print(paste("Jabiru (site 1):",sum(site1)))
  print(paste("Colorados (site 2):",sum(site2)))
  print(paste("Cotove (site 3):",sum(site3)))
  print(paste("Tayrona (site 4):",sum(site4)))
  return(compare_shift_model)
  
}

# LA
output <- compare_shift("log_LA")
summary(output)
saveRDS(output, "Data/compare_shift_LA.RDS")
# LDMC
output <- compare_shift("LDMC_gg")
summary(output)
saveRDS(output, "Data/compare_shift_LDMC.RDS")
# SLA
output <- compare_shift("log_SLA")
summary(output)
saveRDS(output, "Data/compare_shift_SLA.RDS")
# WSG
output <- compare_shift("WSG_gcm3")
summary(output)
saveRDS(output, "Data/compare_shift_WSG.RDS")
# TFW
output <- compare_shift("TFW_um")
summary(output)
saveRDS(output, "Data/compare_shift_TFW.RDS")
# dh
output <- compare_shift("log_dh")
summary(output)
saveRDS(output, "Data/compare_shift_dh.RDS")
# VD
output <- compare_shift("log_VD")
summary(output)
saveRDS(output, "Data/compare_shift_VD.RDS")
# DApit
output <- compare_shift("log_DApit")
summary(output)
saveRDS(output, "Data/compare_shift_DApit.RDS")

# ----- Trait shifts by site (no phenology or growth form) -----

# Calculate shift in trait values from seedling to adult
#   Compare shifts by site, phenology, & growth form
# Arguments:
#   Trait = name of trait (log_SLA, LDMC_gg, etc.)
compare_shift <- function(Trait){
  
  # Create data list
  data <- traits_data %>%
    dplyr::filter(get({{Trait}}) > 0) %>%
    mutate(Trait = scale(get(Trait))) %>%
    dplyr::select(stage, code, site, Trait) %>%
    na.omit() %>%
    mutate(pop_id = paste(site, code, sep="_")) %>%
    group_by(pop_id) %>%
    arrange(stage) %>%
    mutate(change_trait = Trait[1] - Trait[2]) %>%
    slice(1) %>% ungroup %>%
    dplyr::select(-stage) %>%
    na.omit()
  change_trait <- data$change_trait
  site1 <- as.numeric(data$site == "Jabiru")
  site2 <- as.numeric(data$site == "Colorados")
  site3 <- as.numeric(data$site == "Cotove")
  site4 <- as.numeric(data$site == "Tayrona")
  data_list <- list(change_trait = change_trait,
                    site1 = site1, site2 = site2,
                    site3 = site3, site4 = site4,
                    n = nrow(data))
  str(data_list)
  
  # Write model in NIMBLE
  model <- nimbleCode({
    
    # Priors
    a1_Jabi ~ dnorm(0,sd=100) # Intercept for Jabiru
    a2_Colo ~ dnorm(0,sd=100) # Intercept for Colorados
    a3_Coto ~ dnorm(0,sd=100) # Intercept for Cotove
    a4_Tayr ~ dnorm(0,sd=100) # Intercept for Tayrona
    sigma ~ dlnorm(0, sd=100) # SD
    
    # Likelihood
    for (i in 1:n){
      change_trait[i] ~ dnorm(mu[i], sd = sigma)
      mu[i] <- a1_Jabi * site1[i] + a2_Colo * site2[i] +
        a3_Coto * site3[i] + a4_Tayr * site4[i]
    }
    
    # Derived quantities
    diff_Jabi_Colo <- a1_Jabi - a2_Colo
    diff_Jabi_Coto <- a1_Jabi - a3_Coto
    diff_Jabi_Tayr <- a1_Jabi - a4_Tayr
    diff_Colo_Coto <- a2_Colo - a3_Coto
    diff_Colo_Tayr <- a2_Colo - a4_Tayr
    diff_Coto_Tayr <- a3_Coto - a4_Tayr
    
  })
  
  # Set initial values
  inits <- list(a1_Jabi = rnorm(1), a2_Colo = rnorm(1),
                a3_Coto = rnorm(1), a4_Tayr = rnorm(1),
                sigma = rlnorm(1))
  
  # Parameters to estimate
  params <- c("a1_Jabi", "a2_Colo", "a3_Coto", "a4_Tayr",
              "sigma",
              "diff_Jabi_Colo", "diff_Jabi_Coto",
              "diff_Jabi_Tayr", "diff_Colo_Coto",
              "diff_Colo_Tayr", "diff_Coto_Tayr")
  
  # Run model
  compare_shift_model <- nimbleMCMC(code = model, 
                                    constants = data_list,
                                    inits = inits,
                                    monitors = params,
                                    niter = 4000, nburnin = 1000,
                                    nchains = 3,
                                    samplesAsCodaMCMC = TRUE)
  
  print(paste("Sample size:",data_list$n))
  print(paste("Jabiru (site 1):",sum(site1)))
  print(paste("Colorados (site 2):",sum(site2)))
  print(paste("Cotove (site 3):",sum(site3)))
  print(paste("Tayrona (site 4):",sum(site4)))
  return(compare_shift_model)
  
}

# LA
output <- compare_shift("log_LA")
summary(output)
saveRDS(output, "Data/compare_shift2_LA.RDS")
# LDMC
output <- compare_shift("LDMC_gg")
summary(output)
saveRDS(output, "Data/compare_shift2_LDMC.RDS")
# SLA
output <- compare_shift("log_SLA")
summary(output)
saveRDS(output, "Data/compare_shift2_SLA.RDS")
# WSG
output <- compare_shift("WSG_gcm3")
summary(output)
saveRDS(output, "Data/compare_shift2_WSG.RDS")
# TFW
output <- compare_shift("TFW_um")
summary(output)
saveRDS(output, "Data/compare_shift2_TFW.RDS")
# dh
output <- compare_shift("log_dh")
summary(output)
saveRDS(output, "Data/compare_shift2_dh.RDS")
# VD
output <- compare_shift("log_VD")
summary(output)
saveRDS(output, "Data/compare_shift2_VD.RDS")
# DApit
output <- compare_shift("log_DApit")
summary(output)
saveRDS(output, "Data/compare_shift2_DApit.RDS")
