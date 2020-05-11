#
# This script descirbed the multivariate multi level model
# which is used to evaluate relationship between Lp(a)
# reduction and other lipoprotein subclasses changes
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
library(forcats)

# load the multivariate multi level model
AK3 <- stan_model("multivariate_multi_level_regression_model.stan")

#######################################################
# 1. Import data
#######################################################

# load metabolomics data and clinical data
Nightingale <- read_rds("Anitschkow_metabolomics.rds")
Anitschkow <- read_rds("Anischkow_meta.rds")

Anitschkow <- Anitschkow %>% 
  rename(Treatment = `Randomized treatment Group`)

########################################################
# 2. Data preparation
########################################################

dat.lpa <- Anitschkow %>% 
  filter(what == "Lipoprotein(a) nmol/l", Treatment == "AMG145") %>% 
  select(-what, -unit) %>% 
  rename(LPa = y)

Nightingale.nested <- dat.lpa %>% 
  left_join(Nightingale, by = "sampleid") %>% 
  select(patient, visit, Treatment, LPa, what, y, miss, cens) %>% 
  mutate(y = y * 10^9) %>% 
  group_by(what) %>% 
  nest() %>% 
  filter(what %in% c("XXL-VLDL-P", "XL-VLDL-P", 'L-VLDL-P', "M-VLDL-P", "S-VLDL-P", "XS-VLDL-P", 
                     "IDL-P", "L-LDL-P", 'M-LDL-P', "S-LDL-P", 
                     "XL-HDL-P", "L-HDL-P", "M-HDL-P", "S-HDL-P"))

# convert the nested lipoprotein data into stan friendly data
anitschkow_dat <- function(my_dat) {
  subj_dat <- my_dat %>% 
    select(patient, Treatment) %>% 
    distinct()
  
  N <- nrow(my_dat)
  L <- my_dat$LPa
  N_mis <- sum(my_dat$miss)
  N_cen <- sum(my_dat$cens)
  N_obs <- N - N_mis - N_cen
  which_y_mis <- which(my_dat$miss == 1)
  which_y_cen <- which(my_dat$cens == 1)
  which_y_obs <- which(my_dat$miss == 0 & my_dat$cens == 0)
  y_obs <- my_dat$y[which_y_obs]
  N_subj <- nrow(subj_dat)
  subj_id <- as.integer(as.factor(my_dat$patient))
  V2 <- ifelse(my_dat$visit == "V2", 1L, 0L)
  
  my_stan_dat <- list(
    N = N, 
    L = L, 
    N_mis = N_mis, 
    N_cen = N_cen, 
    N_obs = N_obs, 
    which_y_mis = as.array(which_y_mis), 
    which_y_cen = as.array(which_y_cen), 
    which_y_obs = as.array(which_y_obs), 
    y_obs = y_obs, 
    N_subj = N_subj, 
    subj_id = subj_id, 
    V2 = V2
  )
  return(my_stan_dat)
}

lipoproteins <- Nightingale.nested %>% 
  mutate(stan_dat = map(data, anitschkow_dat))

##########################################################
# 3. Run multivariate multi level regression analyses
##########################################################

xiang_cluster <- create_cluster(cores = 10) %>% 
  cluster_library("purrr") %>% 
  cluster_library("rstan") %>% 
  cluster_assign_value("AK3", AK3)

lipoproteins.stan <- lipoproteins %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(stan_output = map(stan_dat, ~sampling(AK3, data = .x, chains = 4, control = list(max_treedepth = 13, adapt_delta = 0.99)))) %>% 
  collect() %>% 
  as_tibble()

#write_rds(lipoproteins.stan, "stan_output_lipoproteins_Lpa.rds")

##########################################################################
# 4. Visualization
##########################################################################

# Figure 2: Medium VLDL particle concentration reduction vs. baseline Lp(a)
AK3.post <- lipoproteins.stan$stan_output[[which(lipoproteins.stan$what == "M-VLDL-P")]]


delta_L <- tidy(AK3.post, pars = "delta_L", estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")
delta_y <- tidy(AK3.post, pars = "delta_y", estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")
aL_subj <- tidy(AK3.post, pars = "aL_subj", estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval")

delta_y$Lpa_baseline <- aL_subj$estimate
delta_y$Lpa_baseline_low <- aL_subj$conf.low
delta_y$Lpa_baseline_high <- aL_subj$conf.high

ggplot(delta_y, aes(Lpa_baseline, estimate)) + 
  geom_pointrange(shape = 1, aes(ymin = conf.low, ymax = conf.high)) + 
  geom_errorbarh(aes(xmin = Lpa_baseline_low, xmax = Lpa_baseline_high)) + 
  geom_hline(yintercept = -0.5, color = "steelblue", linetype = "dashed", size = 1) + 
  labs(x = "Baseline lipoprotein(a) nmol/l", 
       y = "% change in medium VLDL particle") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14))


# Figure S4: Lp(a) reduction vs lipoprotein particle concentration changes
anitschkow_posterior <- function(x, which_par) {
  x %>% 
    mutate(tidy_stan_output = map(stan_output, ~tidy(.x, pars = which_par, estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval"))) %>% 
    unnest(tidy_stan_output)
}

# posterior Pearson correlation coefficient between Lp(a) reduction
# and changes in lipoprotein subclasses
post_rho <- anitschkow_posterior(lipoproteins.stan, "Rho_subj") %>% 
  filter(term == "Rho_subj[2,4]") %>% 
  mutate(what = factor(what, levels = c("XXL-VLDL-P", "XL-VLDL-P", 'L-VLDL-P', "M-VLDL-P", "S-VLDL-P", "XS-VLDL-P", 
                                        "IDL-P", "L-LDL-P", 'M-LDL-P', "S-LDL-P", 
                                        "XL-HDL-P", "L-HDL-P", "M-HDL-P", "S-HDL-P"))) %>% 
  arrange(what)

ggplot(post_rho, aes(fct_rev(what), estimate)) + 
  geom_pointrange(shape = 1, aes(ymin = conf.low, ymax = conf.high)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = NULL, y = "Pearson correlation coefficient") + 
  coord_flip() + 
  theme_classic() + 
  theme(axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14))

post_rho <- anitschkow_posterior(lipoproteins.stan, "Rho_subj") %>% 
  filter(term == "Rho_subj[1,4]") %>% 
  mutate(what = factor(what, levels = c("XXL-VLDL-P", "XL-VLDL-P", 'L-VLDL-P', "M-VLDL-P", "S-VLDL-P", "XS-VLDL-P", 
                                        "IDL-P", "L-LDL-P", 'M-LDL-P', "S-LDL-P", 
                                        "XL-HDL-P", "L-HDL-P", "M-HDL-P", "S-HDL-P"))) %>% 
  arrange(what)
