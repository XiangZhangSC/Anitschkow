#
# This script describes univariate lognormal regression analyses for 
# Lp(a) and lipoprotein subclasses particle concentrations, adjusting
# usage of lipid lowering drugs
#
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(rstan)
library(ggplot2)
library(multidplyr)
library(purrr)
library(broom)
library(tidybayes)
library(stringr)

##########################################################
# 1. Load data
##########################################################

Nightingale <- read_rds("Anitschkow_metabolomics.rds")
Anitschkow <- read_rds("Anischkow_meta.rds")

Nightingale1 <- Nightingale %>% 
  select(sampleid:y)

Anitschkow1 <- Anitschkow %>% 
  rename(Treatment = `Randomized treatment Group`) %>% 
  select(sampleid:Ezetimibe) %>% 
  distinct()

############################################################
# 2. Data preparation
############################################################

Nightingale.nested <- Anitschkow1 %>% 
  left_join(Nightingale1, by = "sampleid") %>% 
  select(-sampleid) %>% 
  spread(visit, y) %>% 
  group_by(what) %>% 
  nest()

anitschkow_dat <- function(my_dat) {
  
  # convert from mol/L to nmol/L (lipoprotein particle concentration)
  y1 <- my_dat$V1 * 10^9
  y2 <- my_dat$V2 * 10^9
  
  # missing observations
  which_y1_mis <- which(is.na(y1))
  N_y1_mis <- length(which_y1_mis)
  
  which_y2_mis <- which(is.na(y2))
  N_y2_mis <- length(which_y2_mis)
  
  # observations equal to 0
  which_y1_cen <- which(y1 == 0)
  N_y1_cen <- length(which_y1_cen)
  y1[which_y1_cen] <- NA
  
  which_y2_cen <- which(y2 == 0)
  N_y2_cen <- length(which_y2_cen)
  y2[which_y2_cen] <- NA
  
  # observed
  which_y1_obs <- which(!is.na(y1))
  N_y1_obs <- length(which_y1_obs)
  
  which_y2_obs <- which(!is.na(y2))
  N_y2_obs <- length(which_y2_obs)


  N <- nrow(my_dat)
  y1_obs <- y1[which_y1_obs]
  y2_obs <- y2[which_y2_obs]
  
  P <- ifelse(my_dat$Treatment == "AMG145", 1L, 0L)
  S <- my_dat$Statine_use
  E <- my_dat$Ezetimibe
  
  my_stan_dat <- list(
    N = N, 
    N_y1_obs = N_y1_obs, 
    N_y1_mis = N_y1_mis, 
    N_y1_cen = N_y1_cen, 
    N_y2_obs = N_y2_obs, 
    N_y2_mis = N_y2_mis, 
    N_y2_cen = N_y2_cen, 
    which_y1_obs = as.array(which_y1_obs), 
    which_y1_mis = as.array(which_y1_mis), 
    which_y1_cen = as.array(which_y1_cen), 
    which_y2_obs = as.array(which_y2_obs), 
    which_y2_mis = as.array(which_y2_mis), 
    which_y2_cen = as.array(which_y2_cen),
    y1_obs = y1_obs, 
    y2_obs = y2_obs,
    P = P, 
    S = S, 
    E = E
  )
  return(my_stan_dat)
}

# nest the Nightinglae metabolomics data per lipoprotein subclass
lipoproteins <- Nightingale.nested %>% 
  filter(what %in% c("XXL-VLDL-P", "XL-VLDL-P", 'L-VLDL-P', "M-VLDL-P", "S-VLDL-P", "XS-VLDL-P", 
                     "IDL-P", "L-LDL-P", 'M-LDL-P', "S-LDL-P", 
                     "XL-HDL-P", "L-HDL-P", "M-HDL-P", "S-HDL-P")) %>% 
  mutate(stan_dat = map(data, anitschkow_dat))

######################################################################
# 3. Run univariate lognormal regression analyses
######################################################################

# load the univariate lognormal regression model
# with additional covariates (statin and ezetimibe)
AK4 <- stan_model("univariate_lognormal_regression_model_adj_lipid_lowering.stan")

xiang_cluster <- create_cluster(cores = 10) %>% 
  cluster_library("purrr") %>% 
  cluster_library("rstan") %>% 
  cluster_assign_value("AK4", AK4)

lipoproteins.stan <- lipoproteins %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(stan_output = map(stan_dat, ~sampling(AK4, data = .x, chains = 4))) %>% 
  collect() %>% 
  as_tibble()

#write_rds(lipoproteins.stan, "stan_lipoproteins_evolocumab_statin_enzitmibe.rds")
#lipoproteins.stan <- read_rds('stan_lipoproteins_evolocumab_statin_enzitmibe.rds')

########################################################################################
# 4. Visualization
########################################################################################

# Figure S1: Lipoprotein subclasses changes adjusting for lipid lowering drug usage
post_delta <- lipoproteins.stan %>% 
  filter(!what %in% c("XXL-VLDL-P", "XL-VLDL-P", "L-VLDL-P")) %>% 
  mutate(result = map(stan_output, ~tidy(.x, pars = c("delta_S0_E0", 
                                                      "delta_S1_E0", 
                                                      "delta_S1_E1"), 
                                         conf.int = TRUE, 
                                         conf.method = "HPDinterval"))) %>% 
  unnest(result)

post_delta1 <- post_delta %>% 
  mutate(what = factor(what, levels = c("M-VLDL-P", 
                                        "S-VLDL-P", 
                                        "XS-VLDL-P", 
                                        "IDL-P", 
                                        "L-LDL-P", 
                                        "M-LDL-P", 
                                        "S-LDL-P", 
                                        "XL-HDL-P", 
                                        "L-HDL-P", 
                                        "M-HDL-P", 
                                        "S-HDL-P")), 
         term = factor(term, levels = c("delta_S0_E0", 
                                        "delta_S1_E0", 
                                        "delta_S1_E1"),
                       labels = c("No medication", 
                                  "Statin",
                                  "Statin and Ezetimibe")))

ggplot(post_delta1) + 
  geom_pointrange(aes(x = what, y = estimate, ymin = conf.low, ymax = conf.high), shape = 1) + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") + 
  labs(x = NULL, y = "Mean difference between treatment and control with 95% credible interval") + 
  coord_flip() + 
  facet_wrap(~term, ncol = 1) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 14), 
        axis.text.x = element_text(size = 12), 
        strip.text = element_text(size = 12))

ggsave("evolocumab_effects_on_lipoproteins.png", 
       width = 200, 
       height = 300, 
       units = "mm", 
       dpi = 300)
 