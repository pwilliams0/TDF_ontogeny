#============================================================#
# COMPARE TRAITS OF SEEDLINGS AND ADULTS, TRAIT CORRELATIONS #
#============================================================#
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

# ========== ANALYZE TRAIT CORRELATIONS ==========

# Calculate correlation between two traits
# Arguments:
#   Stage = "adult" or "seedling"
#   Trait_1 = name of first trait (log_SLA, LDMC_gg, etc.)
#   Trait_2 = name of second trait (log_SLA, LDMC_gg, etc.)
trait_cor <- function(Stage, Trait_1, Trait_2){
  
  # Create data list
  data <- traits_data %>%
    dplyr::filter(stage == Stage,
                  get({{Trait_1}}) > 0,
                  get({{Trait_2}}) > 0) %>%
    dplyr::select(all_of(c(Trait_1,
                           Trait_2))) %>%
    as.matrix
  data_list <- list(data = data,
                    n = nrow(data))
  str(data_list)
  
  # Write model in NIMBLE
  model <- nimbleCode({
    
    # Priors
    mu[1] ~ dnorm(0, 0.001) # Trait_1 mean
    mu[2] ~ dnorm(0, 0.001) # Trait_2 mean
    sigma[1] ~ dunif(0, 1000) # Trait_1 SD
    sigma[2] ~ dunif(0, 1000) # Trait_2 SD
    rho ~ dunif(-1, 1) # Correlation between Trait_1 & Trait_2
    
    # Constructing the covariance matrix
    cov[1,1] <- sigma[1] * sigma[1]
    cov[1,2] <- sigma[1] * sigma[2] * rho
    cov[2,1] <- sigma[1] * sigma[2] * rho
    cov[2,2] <- sigma[2] * sigma[2]
    
    # Likelihood
    for(i in 1:n){
      data[i,1:2] ~ dmnorm(mu[1:2], cov=cov[1:2,1:2])
    }
    
  })
  
  # Set initial values, use classical estimates for mu & sigma
  inits_list <- list(mu = c(mean(data[, 1]), mean(data[, 2])), 
                     sigma = c(sd(data[, 1]), sd(data[, 1])),
                     runif(1,min=-1,max=1))
  
  # Run model
  cor_model <- nimbleMCMC(code = model, 
                          constants = data_list,
                          inits = inits_list,
                          niter = 4000, nburnin = 1000,
                          nchains = 3,
                          samplesAsCodaMCMC = TRUE)
  
  print(paste("Sample size:",data_list$n))
  return(cor_model)
}

# ----- SLA vs. LDMC ----------
# Adult
output <- trait_cor("adult", "log_SLA", "LDMC_gg")
summary(output)
saveRDS(output, "Data/cor_SLA_LDMC_adult.RDS")
# Seedling
output <- trait_cor("seedling", "log_SLA", "LDMC_gg")
summary(output)
saveRDS(output, "Data/cor_SLA_LDMC_sdlg.RDS")

# ----- DApit vs. dh ----------
# Adult
output <- trait_cor("adult", "log_DApit", "log_dh")
summary(output)
saveRDS(output, "Data/cor_DApit_dh_adult.RDS")
# Seedling
output <- trait_cor("seedling", "log_DApit", "log_dh")
summary(output)
saveRDS(output, "Data/cor_DApit_dh_sdlg.RDS")

# ----- DApit vs. VD ----------
# Adult
output <- trait_cor("adult", "log_DApit", "log_VD")
summary(output)
saveRDS(output, "Data/cor_DApit_VD_adult.RDS")
# Seedling
output <- trait_cor("seedling", "log_DApit", "log_VD")
summary(output)
saveRDS(output, "Data/cor_DApit_VD_sdlg.RDS")

# ----- dh vs. VD ----------
# Adult
output <- trait_cor("adult", "log_dh", "log_VD")
summary(output)
saveRDS(output, "Data/cor_dh_VD_adult.RDS")
# Seedling
output <- trait_cor("seedling", "log_dh", "log_VD")
summary(output)
saveRDS(output, "Data/cor_dh_VD_sdlg.RDS")

# ----- Compare intercepts & slopes ----------

# Compare intercepts & slopes for adults & seedlings
# Arguments:
#   Trait_x = name of predictor trait (log_SLA, LDMC_gg, etc.)
#   Trait_y = name of response trait (log_SLA, LDMC_gg, etc.)
trait_linear <- function(Trait_x, Trait_y){
  
  # Create data list
  data <- traits_data %>%
    dplyr::mutate(sdlg = ifelse(stage=="adult", 0, 1)) %>%
    dplyr::filter(get({{Trait_x}}) > 0,
                  get({{Trait_y}}) > 0)
  data_list <- list(sdlg = data$sdlg,
                    Trait_x = as.numeric(scale(data[,Trait_x])),
                    Trait_y = as.numeric(scale(data[,Trait_y])),
                    n = nrow(data))
  str(data_list)
  
  # Write model in NIMBLE
  model <- nimbleCode({
    
    # Priors
    alpha1 ~ dnorm(0,sd=100) # Sdlg intercept
    alpha2 ~ dnorm(0,sd=100) # Adult intercept
    beta1 ~ dnorm(0,sd=100) # Sdlg slope
    beta2 ~ dnorm(0,sd=100) # Adult slope
    sigma ~ dunif(0, 100) # SD
    
    # Likelihood
    for (i in 1:n){
      Trait_y[i] ~ dnorm(mu[i], sd = sigma)
      mu[i] <- alpha1 * sdlg[i] + beta1 * Trait_x[i] * sdlg[i] +
        alpha2 * (1-sdlg[i]) + beta2 * Trait_x[i] * (1-sdlg[i])
    }
    
    # Derived quantities
    alpha_diff <- alpha1 - alpha2 # Difference in intercepts
    beta_diff <- beta1 - beta2 # Difference in slopes
    
  })
  
  # Set initial values
  inits <- list(alpha1 = rnorm(1), alpha2 = rnorm(1),
                beta1 = rnorm(1), beta2 = rnorm(1),
                sigma = rlnorm(1))
  
  # Parameters to estimate
  params <- c("alpha1", "alpha2", "beta1", "beta2",
              "sigma", "alpha_diff", "beta_diff")
  
  # Run model
  linear_model <- nimbleMCMC(code = model, 
                            constants = data_list,
                            inits = inits,
                            monitors = params,
                            niter = 4000, nburnin = 1000,
                            nchains = 3,
                            samplesAsCodaMCMC = TRUE)
  
  print(paste("Sample size:",data_list$n))
  return(linear_model)
}

# dh vs. VD
output <- trait_linear(Trait_x="log_VD", Trait_y="log_dh")
summary(output)
saveRDS(output, "Data/linear_dh_VD.RDS")

# SLA vs. LDMC
output <- trait_linear(Trait_x="LDMC_gg", Trait_y="log_SLA")
summary(output)
saveRDS(output, "Data/linear_SLA_LDMC.RDS")
