#==========================================================#
# COMPARE TRAITS OF SEEDLINGS AND ADULTS, PCA WTIH VARIMAX #
#==========================================================#
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

# ========== RUN PCA WTIH VARIMAX ==========

# NOTE: Edit lines 30, 223, & 242 for PCA of either seedlings, adults, or combined
#       Edit lines 34, 44, & 48 to exclude DApit from combined seedling/adult PCA 

# ----- Calculate covariance matrix among traits ----------

data <- traits_data %>%
  filter(stage == "adult") %>%
  column_to_rownames("id") %>%
  dplyr::select(log_LA, LDMC_gg, log_SLA, WSG_gcm3,
                TFW_um, log_dh, log_VD,
                log_DApit
                ) %>%
  na.omit() %>%
  mutate(LA = as.numeric(scale(log_LA)),
         LDMC = as.numeric(scale(LDMC_gg)),
         SLA = as.numeric(scale(log_SLA)),
         WSG = as.numeric(scale(WSG_gcm3)),
         TFW = as.numeric(scale(TFW_um)),
         dh = as.numeric(scale(log_dh)),
         VD = as.numeric(scale(log_VD)),
         DApit = as.numeric(scale(log_DApit))
         ) %>%
  dplyr::select(LA, LDMC, SLA, WSG,
                TFW, dh, VD,
                DApit
                ) 

# Create data list
V <- ncol(data) # number of traits
data_list <- list(data = as.matrix(data),
                  n = nrow(data),
                  V = V,
                  Rmat = diag(V),
                  df =  V + 1,
                  mu_prior = rep(0,V))
str(data_list)

# Write model in NIMBLE
model <- nimbleCode({
  
  # Priors
  mu[1:V] ~ dmnorm(mean=mu_prior[1:V], prec=prec[1:V,1:V]) #  Multinormal means
  prec[1:V,1:V] ~ dwish(R = Rmat[1:V,1:V], df = df) # Cov. matrix using precision
  
  # Construct covariance matrix from precision
  cov[1:V,1:V] <- inverse(prec[1:V,1:V])
  
  # Likelihood
  for(i in 1:n){
    data[i,1:V] ~ dmnorm(mean=mu[1:V], prec=prec[1:V,1:V])
  }
  
})

# Set initial values
inits_list <- list(
  mu = rnorm(V),
  prec = rWishart(1, df = data_list$df,
                  Sigma = inverse(data_list$Rmat))[,,1])

# Run model
pca_model <- nimbleMCMC(code = model, 
                        constants = data_list,
                        monitors = c("cov", "mu"),
                        inits = inits_list,
                        niter = 4000, nburnin = 1000,
                        nchains = 3,
                        samplesAsCodaMCMC = TRUE)

print(paste("Number of populations:",data_list$n))
print(paste("Number of traits:",data_list$V))

summary(pca_model)

# ----- Calculate loading matrices with varimax rotation ----------

# Create list of covariance matrices from all samples
samp <- rbind(pca_model$chain1, pca_model$chain2,
              pca_model$chain3)
V <- ncol(data)
samp.cov <- samp[,1:(V^2)]
samp.cov <- split(samp.cov, seq(nrow(samp.cov)))
names(samp.cov) <- NULL

# Function used by lapply to calculate loadings
#    Convert covariance matrix to a correlation matrix
#    Standardize loadings
load.extract <- function(cov, V, data){
  covm <- cov2cor(matrix(cov, V, V))
  eigens <- eigen(covm)
  loadings <- eigen(covm)$vectors %*% sqrt(diag(eigens$values,
                                                nrow = length(eigens$values)))
  row.names(loadings) <- names(data)
  colnames(loadings) <- paste("Comp.", 1:V, sep="")
  return(loadings)
}

# Calculate loadings of each covariance matrix in list
loadings.original <- lapply(X=samp.cov, FUN=load.extract, V=V, data=data)

# Calculate loadings with varimax rotation 
n_axes = 2
do_varimax <- function(loadings, n_axes){
  rotated <- unclass(varimax(loadings[,1:n_axes])$loadings)
  return(rotated)
}
loadings.varimax <- lapply(X=loadings.original,
                           FUN=do_varimax, n_axes=n_axes)

# Identify trait with highest mean loading on each axis
axis_traits <- c()
n_samp <- length(loadings.varimax)
tmp <- matrix(ncol=V, nrow=n_samp)
colnames(tmp) <- colnames(data)
for(j in 1:n_axes){
  for(i in 1:n_samp){
    tmp[i,] <- abs(loadings.varimax[[i]][,j])
  }
  axis_traits[j] <- names(sort(colMeans(tmp),
                               decreasing=TRUE))[1]
}

# Fix the signs and loadings of axes with reference to axis_traits
loadings.fixed <- loadings.varimax
for(i in 1:n_samp){
  mat <- matrix(ncol=n_axes,nrow=V)
  rownames(mat) <- colnames(data)
  assigned_axes <- c()
  for(j in 1:n_axes){
    # For trait with highest loading on axis j overall,
    #   identify axis with highest loading in this iteration
    # Make this axis the new axis j
    df <- abs(loadings.varimax[[i]])
    pc_val <- df[axis_traits[j],]
    pc_val[assigned_axes] <- 0
    col <- as.integer(which(pc_val == max(pc_val))[1])
    mat[,j] <- loadings.varimax[[i]][,col]
    # Start with assigning the first axis
    #   then move through remaining axes
    assigned_axes <- c(assigned_axes,col)
    # Fix signs so that trait with highest loading on axis j is always positive
    loadings.fixed[[i]][,j] <- ifelse(mat[axis_traits[j], j] > 0,
                                      list(mat[,j]),
                                      list(-mat[,j]))[[1]]
    }
}

# ----- Calculate variance explained by each axis ----------

prop_var <- matrix(ncol=n_axes,
                   nrow=length(loadings.fixed))
for(i in 1:length(loadings.fixed)){
  tmp <- loadings.fixed[[i]]
  prop_var[i,1:n_axes] <- colSums(tmp^2)/nrow(tmp)
  colnames(prop_var) <- colnames(tmp)
}
prop_var <- as.data.frame(prop_var)
summary(prop_var)
quantile(prop_var$Comp.1, probs = c(.025,.975))
quantile(prop_var$Comp.2, probs = c(.025,.975))

# ----- Calculate mean loadings with median & 95% CI ----------

# Axis 1
loadings_1 <- matrix(nrow=n_samp, ncol=V)
for(i in 1:n_samp){
  loadings_1[i,1:V] <- loadings.fixed[[i]][,1]
  colnames(loadings_1) <- rownames(loadings.fixed[[i]])
}
loadings_1 <- data.frame(loadings_1)
# Axis 2
loadings_2 <- matrix(nrow=n_samp, ncol=V)
for(i in 1:n_samp){
  loadings_2[i,1:V] <- loadings.fixed[[i]][,2]
  colnames(loadings_2) <- rownames(loadings.fixed[[i]])
}
loadings_2 <- data.frame(loadings_2)
# Summary
loadings_summary <- data.frame("Trait"=colnames(loadings_1),
                               "Comp1.Mean"=NA,
                               "Comp1.Q02.5"=NA,
                               "Comp1.Q25"=NA,
                               "Comp1.Median"=NA,
                               "Comp1.Q75"=NA,
                               "Comp1.Q97.5"=NA,
                               "Comp2.Mean"=NA,
                               "Comp2.Q02.5"=NA,
                               "Comp2.Q25"=NA,
                               "Comp2.Median"=NA,
                               "Comp2.Q75"=NA,
                               "Comp2.Q97.5"=NA)
for(i in 1:V){
  loadings_summary[i,2] <- mean(loadings_1[,i])
  loadings_summary[i,3:7] <- quantile(loadings_1[,i],
                                      probs=c(.025,.25,.5,.75,.975))
  loadings_summary[i,8] <- mean(loadings_2[,i])
  loadings_summary[i,9:13] <- quantile(loadings_2[,i],
                                      probs=c(.025,.25,.5,.75,.975))
}
write.csv(loadings_summary, "Data/pca_loadings_adult.csv")

# ----- Calculate scores ----------

# Calculate median loadings
loadings.median <- data.frame(
  Comp.1 = sapply(loadings_1, FUN=median),
  Comp.2 = sapply(loadings_2, FUN=median))
# Calculate variance explained using median loadings
colSums(loadings.median^2)/nrow(loadings.median)
# Calculate scores
scores <- as.matrix(data) %*% as.matrix(loadings.median)
# Join with full data
scores_df <- as.data.frame(scores) %>%
  rownames_to_column("id") %>%
  separate(id, c("code","site","stage","X")) %>%
  dplyr::select(-X) %>%
  left_join(traits_data, by=c("code","site","stage")) %>%
  arrange(stage)
write.csv(scores_df, "Data/pca_scores_adult.csv")
