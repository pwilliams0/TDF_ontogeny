#================#
# CREATE FIGURES #
#================#
library(tidyverse)
library(MCMCvis)

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
         id = paste(code,site,stage,row_id))

# ========== Fig 1: ONTOGENETIC TRAIT SHIFTS ==========

plot_shift <- function(Trait, log, col, ylab){
  data <- traits_data %>%
    dplyr::filter(get({{Trait}}) > 0) %>%
    dplyr::select(stage, code, site, {Trait},
                  leaf_phenology, life_form) %>%
    na.omit() %>%
    mutate(pop_id = paste(site, code, sep="_")) %>%
    group_by(pop_id) %>%
    arrange(stage) %>%
    mutate(change_trait = get({{Trait}})[1] - get({{Trait}})[2]) %>%
    ungroup %>%
    mutate(stage = factor(stage, levels=c("seedling","adult"))) %>%
    na.omit()
  
  means <- data %>%
    group_by(stage) %>%
    mutate(trait = get({{Trait}}),
           means = ifelse(log==TRUE,
                          exp(mean(log(trait))),
                          mean(trait)),
           group = 1) %>%
    slice(1)
  
  ggplot() +
    geom_path(data=data, aes(x=stage, y=get(Trait), group=pop_id),
              color="gray60", size=.2) +
    geom_path(data=means, aes(x=stage, y=means, group=group),
              color = col, size = 1) +
    geom_point(data=data, aes(x = stage, y = get(Trait)),
               size = .5, alpha=0.5) +
    {if(log)scale_y_continuous(trans = "log10")} +
    theme_classic() +
    ylab(ylab) +
    theme(axis.text = element_text(colour = "black", size=10),
          axis.title = element_text(colour = "black", size=12),
          axis.title.x = element_blank(),
          line = element_line(size = 0.25))
}

svg("Results/shift_LA.svg", width = 2, height = 2)
plot_shift("LA_cm2", log=TRUE, "#b2182b",
           expression(paste("LA (",cm^2,")")))
dev.off()

svg("Results/shift_LDMC.svg", width = 2, height = 2)
plot_shift("LDMC_gg", log=FALSE, "#b2182b",
           expression(paste("LDMC (g ",g^-1,")")))
dev.off()

svg("Results/shift_SLA.svg", width = 2, height = 2)
plot_shift("SLA_cm2g", log=TRUE,  "#2166ac",
           expression(paste("SLA (",cm^2," ",g^-1,")")))
dev.off()

svg("Results/shift_WSG.svg", width = 2, height = 2)
plot_shift("WSG_gcm3", log=FALSE, "#b2182b",
           expression(paste("WSG (g ",cm^-3,")")))
dev.off()

svg("Results/shift_TFW.svg", width = 2, height = 2)
plot_shift("TFW_um", log=FALSE, "black",
           "TFW (μm)")
dev.off()

svg("Results/shift_dh.svg", width = 2, height = 2)
plot_shift("dh_um", log=TRUE, "#b2182b",
           expression(paste(d[h]," (μm)")))
dev.off()

svg("Results/shift_VD.svg", width = 2, height = 2)
plot_shift("VD_nomm2", log=TRUE, "black",
           expression(paste("VD (vessels ",mm^-2,")")))
dev.off()

svg("Results/shift_DApit.svg", width = 2, height = 2)
plot_shift("DApit_um", log=TRUE, "#2166ac",
           expression(paste(DA[pit]," (μm)")))
dev.off()

# ========== Fig 2: TRAIT SHIFTS BY SITE ==========

# Load model outputs for shift in trait values by site
#   Accounting for phenology and growth form
LA <- MCMCsummary(readRDS("Data/compare_shift_LA.RDS"))
LDMC <- MCMCsummary(readRDS("Data/compare_shift_LDMC.RDS"))
SLA <- MCMCsummary(readRDS("Data/compare_shift_SLA.RDS"))
WSG <- MCMCsummary(readRDS("Data/compare_shift_WSG.RDS"))
TFW <- MCMCsummary(readRDS("Data/compare_shift_TFW.RDS"))
dh <- MCMCsummary(readRDS("Data/compare_shift_dh.RDS"))
VD <- MCMCsummary(readRDS("Data/compare_shift_VD.RDS"))
DApit <- MCMCsummary(readRDS("Data/compare_shift_DApit.RDS"))

# Extract mean and 95% CI of each site for each trait
site_effects <- function(Trait){
  df <- MCMCvis::MCMCsummary(readRDS(
    paste("Data/compare_shift_",Trait,".RDS", sep="")))[1:4,] %>%
    as.matrix() %>% as.data.frame() %>%
    rownames_to_column("site") %>%
    mutate(Trait = Trait)
  return(df)
}

site_shifts <- bind_rows(site_effects("LA"),
                         site_effects("LDMC"),
                         site_effects("SLA"),
                         site_effects("WSG"),
                         site_effects("TFW"),
                         site_effects("dh"),
                         site_effects("VD"),
                         mutate(site_effects("DApit"),
                                Trait = "DA pit")) %>%
  mutate(Trait = factor(Trait,levels=c("LA","LDMC","SLA",
                                       "WSG","TFW","dh",
                                       "VD","DA pit")),
         Site = case_when(site == "a1_Jabi" ~ "Jabiru",
                          site == "a2_Colo" ~ "Colorados",
                          site == "a3_Coto" ~ "Cotove",
                          site == "a4_Tayr" ~ "Tayrona"),
         Site = factor(Site,levels=c("Jabiru","Colorados",
                                     "Cotove","Tayrona"))) %>%
  group_by(Trait) %>%
  arrange(Site) %>%
  mutate(order = as.factor(row_number()))

# Plot results
svg("Results/shift_by_site.svg", width = 7.5, height = 4)
ggplot(site_shifts, aes( x=Trait)) +
  geom_hline(yintercept=0, color="gray", linetype="longdash") +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`, y = mean,
                      alpha = order,
                      color = Site, fill = Site, shape = Site),
                  size = 0.5,
                  position =  position_dodge(width=0.8)) +
  scale_alpha_manual(values=c(1,1,1,1), guide = "none") +
  scale_shape_manual(values=c(23,22,21,24)) +
  scale_color_manual(values=c("#2166ac","#92c5de","#f4a582","#b2182b")) +
  scale_fill_manual(values=c("#2166ac","#92c5de","#f4a582","#b2182b")) +
  theme_classic() +
  coord_cartesian(ylim=c(min(site_shifts$`2.5%`),
                         max(site_shifts$`97.5%`) + .5)) +
  ylab("Shift in standardized trait value") +
  theme(legend.position = "top",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25),
        axis.text = element_text(colour = "black", size=10),
        axis.title = element_text(colour = "black", size=12),
        axis.title.x = element_blank(),
        line = element_line(size = 0.25))
dev.off()

# ========== Fig 3: TRAIT PCA ==========

# Plot loadings function
plot_loadings <- function(loadings_summary){
  ggplot(loadings_summary) +
    geom_hline(yintercept = 0, color="gray90") +
    geom_vline(xintercept = 0, color="gray90") +
    geom_segment(color="#DDAA33", arrow = arrow(length=unit(.05, 'in')),
                 linewidth=0.5,
                 aes(x=0, y=0, xend=Comp1.Median, yend=Comp2.Median)) +
    geom_linerange(aes(ymin=Comp2.Q02.5, ymax=Comp2.Q97.5,
                       x=Comp1.Median), linewidth=0.1) +
    geom_linerange(aes(xmin=Comp1.Q02.5, xmax=Comp1.Q97.5,
                       y=Comp2.Median), linewidth=0.1) +
    geom_linerange(aes(ymin=Comp2.Q25, ymax=Comp2.Q75,
                       x=Comp1.Median), linewidth=0.5) +
    geom_linerange(aes(xmin=Comp1.Q25, xmax=Comp1.Q75,
                       y=Comp2.Median), linewidth=0.5) +
    theme_classic() +
    coord_equal(ratio=1, xlim=c(-1,1), ylim=c(-1,1)) +
    theme(axis.text = element_text(colour = "black", size=9),
          axis.title = element_text(colour = "black", size=10),
          line = element_line(linewidth = 0.25)) +
    # geom_text(aes(x=Comp1.Median+.1, y=Comp2.Median+.1,
    #               label=Trait), size=3) +
    xlab("Component 1") + ylab("Component 2")
}

# Adult loadings
svg("Results/pca_loadings_adult.svg", width = 2.5, height = 2.5)
plot_loadings(read.csv("Data/pca_loadings_adult.csv"))
dev.off()

# Seedling loadings
svg("Results/pca_loadings_sdlg.svg", width = 2.5, height = 2.5)
plot_loadings(read.csv("Data/pca_loadings_sdlg.csv"))
dev.off()

# Combined loadings
svg("Results/pca_loadings_both.svg", width = 2.5, height = 2.5)
plot_loadings(read.csv("Data/pca_loadings_both.csv"))
dev.off()

# Plot scores function
plot_scores <- function(scores, loadings, palette){
  range_x <- range(scores$Comp.1)
  range_y <- range(scores$Comp.2)
  ylimits <- if(diff(range_y) > diff(range_x)) {range_y} else {
    mean(range_y) + c(-1, +1) * diff(range_x)/2
  }
  xlimits <- if(diff(range_x) > diff(range_y)) {range_x} else {
    mean(range_x) + c(-1, +1) * diff(range_y)/2
  }
  
  ggplot() +
    geom_point(data = scores, aes(x = Comp.1, color=stage,
                                     y = Comp.2),
               shape = "circle", size = 1, alpha=.5) +
    geom_segment(data = loadings, color="black", linewidth=0.25,
                 aes(x=0, y=0, xend=Comp1.Median*5, yend=Comp2.Median*5),
                 arrow = arrow(length=unit(.05, 'in'))) +
    scale_color_manual(values=palette) +
    # geom_text(data = loadings, size=9/.pt,
    #           aes(x=Comp1.Median*6, y=Comp2.Median*6, label=Trait)) +
    theme_classic() +
    coord_equal(ratio=1, ylim=ylimits, xlim=xlimits) +
    theme(axis.text = element_text(colour = "black", size=9),
          axis.title = element_text(colour = "black", size=10),
          legend.position = "none",
          line = element_line(linewidth = 0.25)) +
    xlab("Component 1") + ylab("Component 2")
}

# Adult scores
svg("Results/pca_scores_adult.svg", width = 2.5, height = 2.5)
plot_scores(read.csv("Data/pca_scores_adult.csv"),
            read.csv("Data/pca_loadings_adult.csv"),
            "#F8766D")
dev.off()

# Seedling scores
svg("Results/pca_scores_sdlg.svg", width = 2.5, height = 2.5)
plot_scores(read.csv("Data/pca_scores_sdlg.csv"),
            read.csv("Data/pca_loadings_sdlg.csv"),
            "#00BFC4")
dev.off()

# Combined scores
svg("Results/pca_scores_both.svg", width = 2.5, height = 2.5)
plot_scores(read.csv("Data/pca_scores_both.csv"),
            read.csv("Data/pca_loadings_both.csv"),
            c("#F8766D","#00BFC4"))
dev.off()

# ========== Fig 4: TRAIT CORRELATIONS ==========

my_theme <- list(geom_point(shape = "circle", size = .5,
                            alpha = .5),
                 theme_classic(),
                 theme(axis.text = element_text(colour = "black", size=10),
                       axis.title = element_text(colour = "black", size=12),
                       line = element_line(linewidth = 0.25),
                       legend.position = "none"))

# DApit vs. dh
svg("Results/DApit_vs_dh.svg", width = 2.5, height = 2.5)
ggplot(traits_data %>% arrange(stage)) +
  aes(x = DApit_um, y = dh_um, color = stage) +
  geom_smooth(method="lm", alpha = 0.4, linewidth=.25) +
  my_theme +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(x = expression(paste(DA[pit]," (μm)")),
       y = expression(paste(d[h]," (μm)")))
dev.off()

# DApit vs. VD
svg("Results/DApit_vs_VD.svg", width = 2.5, height = 2.5)
ggplot(traits_data %>% arrange(stage)) +
  aes(x = DApit_um, y = VD_nomm2, color = stage) +
  geom_smooth(method="lm", alpha = 0.4, linewidth=.25) +
  my_theme +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(x = expression(paste(DA[pit]," (μm)")),
       y = expression(paste("VD (vessels ",mm^-2,")")))
dev.off()

# VD vs. dh
svg("Results/VD_vs_dh.svg", width = 2.5, height = 2.5)
ggplot(traits_data %>% arrange(stage)) +
  aes(x = VD_nomm2, y = dh_um, color = stage) +
  geom_smooth(method="lm", alpha = 0.4, linewidth=.25) +
  my_theme +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(x = expression(paste("VD (vessels ",mm^-2,")")),
       y = expression(paste(d[h]," (μm)")))
dev.off()

# LDMC vs. SLA
svg("Results/LDMC_vs_SLA.svg", width = 2.5, height = 2.5)
ggplot(traits_data %>% arrange(stage)) +
  aes(x = LDMC_gg, y = SLA_cm2g, color = stage) +
  geom_smooth(method="lm", alpha = 0.4, linewidth=.25) +
  my_theme +
  scale_y_continuous(trans = "log10") +
  labs(x = expression(paste("LDMC (g ",g^-1,")")),
       y = expression(paste("SLA (",cm^2," ",g^-1,")")))
dev.off()

# ========== Fig S1: TRAIT SHIFTS BY SITE ==========

# Load model outputs for shift in trait values by site
#   Accounting for phenology and growth form
LA <- MCMCsummary(readRDS("Data/compare_shift2_LA.RDS"))
LDMC <- MCMCsummary(readRDS("Data/compare_shift2_LDMC.RDS"))
SLA <- MCMCsummary(readRDS("Data/compare_shift2_SLA.RDS"))
WSG <- MCMCsummary(readRDS("Data/compare_shift2_WSG.RDS"))
TFW <- MCMCsummary(readRDS("Data/compare_shift2_TFW.RDS"))
dh <- MCMCsummary(readRDS("Data/compare_shift2_dh.RDS"))
VD <- MCMCsummary(readRDS("Data/compare_shift2_VD.RDS"))
DApit <- MCMCsummary(readRDS("Data/compare_shift2_DApit.RDS"))

# Extract mean and 95% CI of each site for each trait
site_effects <- function(Trait){
  df <- MCMCvis::MCMCsummary(readRDS(
    paste("Data/compare_shift2_",Trait,".RDS", sep="")))[1:4,] %>%
    as.matrix() %>% as.data.frame() %>%
    rownames_to_column("site") %>%
    mutate(Trait = Trait)
  return(df)
}

site_shifts <- bind_rows(site_effects("LA"),
                         site_effects("LDMC"),
                         site_effects("SLA"),
                         site_effects("WSG"),
                         site_effects("TFW"),
                         site_effects("dh"),
                         site_effects("VD"),
                         mutate(site_effects("DApit"),
                                Trait = "DA pit")) %>%
  mutate(Trait = factor(Trait,levels=c("LA","LDMC","SLA",
                                       "WSG","TFW","dh",
                                       "VD","DA pit")),
         Site = case_when(site == "a1_Jabi" ~ "Jabiru",
                          site == "a2_Colo" ~ "Colorados",
                          site == "a3_Coto" ~ "Cotove",
                          site == "a4_Tayr" ~ "Tayrona"),
         Site = factor(Site,levels=c("Jabiru","Colorados",
                                     "Cotove","Tayrona"))) %>%
  group_by(Trait) %>%
  arrange(Site) %>%
  mutate(order = as.factor(row_number()))

# Plot results
svg("Results/shift2_by_site.svg", width = 7.5, height = 4)
ggplot(site_shifts, aes( x=Trait)) +
  geom_hline(yintercept=0, color="gray", linetype="longdash") +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`, y = mean,
                      alpha = order,
                      color = Site, fill = Site, shape = Site),
                  size = 0.5,
                  position =  position_dodge(width=0.8)) +
  scale_alpha_manual(values=c(1,1,1,1), guide = "none") +
  scale_shape_manual(values=c(23,22,21,24)) +
  scale_color_manual(values=c("#2166ac","#92c5de","#f4a582","#b2182b")) +
  scale_fill_manual(values=c("#2166ac","#92c5de","#f4a582","#b2182b")) +
  theme_classic() +
  coord_cartesian(ylim=c(min(site_shifts$`2.5%`),
                         max(site_shifts$`97.5%`) + .5)) +
  ylab("Shift in standardized trait value") +
  theme(legend.position = "top",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25),
        axis.text = element_text(colour = "black", size=10),
        axis.title = element_text(colour = "black", size=12),
        axis.title.x = element_blank(),
        line = element_line(size = 0.25))
dev.off()
