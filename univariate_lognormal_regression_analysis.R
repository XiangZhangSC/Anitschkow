#
# This script dscribed the univariate lognormal regression analyses
# for lipoportein subclasses and LP(a), comparing evolocumab and 
# placebo group
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
library(GGally)
library(forcats)

############################################
# Import data
############################################

# load metabolomics data and clinical data
Nightingale <- read_rds("Anitschkow_metabolomics.rds")
Anitschkow <- read_rds("Anischkow_meta.rds")


Anitschkow <- Anitschkow %>% 
  rename(Treatment = `Randomized treatment Group`)

sample_visit <- Anitschkow %>% 
  select(sampleid:Treatment, Statine_use, Ezetimibe) %>% 
  distinct()

############################################
# 1. Analysis-ready data preparation
############################################

# reshpae clinical data (LP[a], LDLc and TG) for
# lognormal regression model for assessing the effects
# of evolocumab
anitschkow_biochemical_dat <- function(df, which_biochemical_compound) {
  df %>% 
    filter(what == which_biochemical_compound) %>% 
    select(patient, visit, Treatment, Statine_use, Ezetimibe, what, y) %>% 
    spread(visit, y) %>% 
    mutate(v1_mis = ifelse(is.na(V1), 1L, 0L), 
           v1_cen = ifelse(V1 == 0, 1L, 0L), 
           v2_mis = ifelse(is.na(V2), 1L, 0L), 
           v2_cen = ifelse(V2 == 0, 1L, 0L))
}

d.lpa <- anitschkow_biochemical_dat(Anitschkow, "Lipoprotein(a) nmol/l")
d.ldlc <- anitschkow_biochemical_dat(Anitschkow, "LDL cholesterol mmol/l")
d.tg <- anitschkow_biochemical_dat(Anitschkow, "Triglyceride mmol/l")

# nest the Nightingale metabolomics data based on the metabolites
Nightingale.nested <- sample_visit %>% 
  left_join(Nightingale, by = "sampleid") %>% 
  select(patient:unit) %>% 
  spread(visit, y) %>% 
  mutate(v1_mis = ifelse(is.na(V1), 1L, 0L), 
         v1_cen = ifelse(!is.na(V1) & V1 == 0, 1L, 0L), 
         v2_mis = ifelse(is.na(V2), 1L, 0L), 
         v2_cen = ifelse(!is.na(V2) & V2 == 0, 1L, 0L)) %>% 
  group_by(what) %>% 
  nest()

# convert the nested NIghtingale metabolomics data into
# stan friendly data
anitschkow_dat <- function(df, scalar) {
  
  metabolite.v1 <- df$V1 * scalar
  metabolite.v2 <- df$V2 * scalar
  
  E <- df$Ezetimibe
  S <- df$Statine_use
  
  N <- nrow(df)
  N_y1_mis <- sum(df$v1_mis)
  N_y1_cen <- sum(df$v1_cen)
  N_y1_obs <- N - N_y1_mis - N_y1_cen
  
  N_y2_mis <- sum(df$v2_mis)
  N_y2_cen <- sum(df$v2_cen)
  N_y2_obs <- N - N_y2_mis - N_y2_cen
  
  which_y1_mis <- which(df$v1_mis == 1)
  which_y1_cen <- which(df$v1_cen == 1)
  which_y1_obs <- which(df$v1_mis == 0 & df$v1_cen == 0)
  
  metabolite.v1[which_y1_cen] <- NA
  
  which_y2_mis <- which(df$v2_mis == 1)
  which_y2_cen <- which(df$v2_cen == 1)
  which_y2_obs <- which(df$v2_mis == 0 & df$v2_cen == 0)
  
  metabolite.v2[which_y2_cen] <- NA
  
  N_group <- length(unique(df$Treatment))
  group_id <- as.integer(factor(df$Treatment, levels = c("PLACEBO", "AMG145")))
  
  y1_obs <- metabolite.v1[which_y1_obs]
  y2_obs <- metabolite.v2[which_y2_obs]
  
  AK.dat <- list(
    N = N, 
    N_y1_mis = N_y1_mis, 
    N_y1_cen = N_y1_cen, 
    N_y1_obs = N_y1_obs, 
    N_y2_mis = N_y2_mis, 
    N_y2_cen = N_y2_cen, 
    N_y2_obs = N_y2_obs, 
    which_y1_obs = as.array(which_y1_obs), 
    which_y1_mis = as.array(which_y1_mis), 
    which_y1_cen = as.array(which_y1_cen), 
    which_y2_obs = as.array(which_y2_obs), 
    which_y2_mis = as.array(which_y2_mis), 
    which_y2_cen = as.array(which_y2_cen), 
    N_group = N_group, 
    y1_obs = y1_obs, 
    y2_obs = y2_obs, 
    group_id = group_id, 
    E = E, 
    S = S
  )
  
  return(AK.dat)
}

dat.lpa.stan <- anitschkow_dat(d.lpa,1)
dat.ldlc.stan <- anitschkow_dat(d.ldlc,1)
dat.tg.stan <- anitschkow_dat(d.tg,1)

# generate stan friendly lipoprotein particle concentration data
lipoproteins <- Nightingale.nested %>% 
  filter(what %in% c("XXL-VLDL-P", "XL-VLDL-P", 'L-VLDL-P', "M-VLDL-P", "S-VLDL-P", "XS-VLDL-P", 
                     "IDL-P", "L-LDL-P", 'M-LDL-P', "S-LDL-P", 
                     "XL-HDL-P", "L-HDL-P", "M-HDL-P", "S-HDL-P")) %>% 
  mutate(stan_dat = map(data, ~anitschkow_dat(.x, 10^9)))

# generate stan friendly lipoprotein lipid composition data
my_lipoproteins <- c("XXL-VLDL", "XL-VLDL", 'L-VLDL', "M-VLDL", "S-VLDL", "XS-VLDL", 
                     "IDL", "L-LDL", 'M-LDL', "S-LDL", 
                     "XL-HDL", "L-HDL", "M-HDL", "S-HDL")

my_tg <- paste(my_lipoproteins, "TG", sep = "-")
my_ce <- paste(my_lipoproteins, "CE", sep = "-")
my_fc <- paste(my_lipoproteins, "FC", sep = "-")
my_pl <- paste(my_lipoproteins, "PL", sep = "-")

lipoprotein_lipids <- Nightingale.nested %>% 
  filter(what %in% c(my_tg, my_ce, my_fc, my_pl)) %>% 
  mutate(stan_dat = map(data, ~anitschkow_dat(.x, 1)))

# generate stan friendly amino acid data
AAs <- Nightingale.nested %>% 
  filter(what %in% c("Ala", "Gln", "His", "Ile", "Leu", "Val", "Phe", "Tyr")) %>% 
  mutate(stan_dat = map(data, ~anitschkow_dat(.x, 10^3)))

# generate stan friendly ketone body and glycolysis related metabolites
polar_metabolites <- Nightingale.nested %>% 
  filter(what %in% c("AcAce", "Ace", "Alb", "bOHBut", "Cit", "Crea", "Glc", "Gp", "Lac")) %>% 
  mutate(stan_dat = map(data, ~anitschkow_dat(.x, 1)))

# generate stan friendly fatty acids data
ffas <- Nightingale.nested %>% 
  filter(what %in% c("MUFA", "DHA", "FAw3", "FAw6", "LA", "PUFA", "SFA", "TotFA")) %>% 
  mutate(stan_dat = map(data, ~anitschkow_dat(.x, 1)))

############################################
# 2. Run lognormal regression analyses
############################################

# load the lognormal regression model
AK <- stan_model("univariate_lognormal_regression.stan")

# Lipoprotein(a)
LPa.post <- sampling(AK, data = dat.lpa.stan, chains = 4, cores = 4L)
tidy(LPa.post, conf.int = TRUE, conf.method = "HPDinterval")

# LDL-C
LDLc.post <- sampling(AK, data = dat.ldlc.stan, chains = 4, cores = 4L)
tidy(LDLc.post, conf.int = TRUE, conf.method = "HPDinterval")

# TG
TG.post <- sampling(AK, data = dat.tg.stan, chains = 4, cores = 4L)
tidy(TG.post, conf.int = TRUE, conf.method = "HPDinterval")

# define the cluster
xiang_cluster <- create_cluster(cores = 10) %>% 
  cluster_library("purrr") %>% 
  cluster_library("rstan") %>% 
  cluster_assign_value("AK", AK)

# regression analyses for lipoprotein particle concentrations
lipoproteins.stan <- lipoproteins %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(stan_output = map(stan_dat, ~sampling(AK, data = .x, chains = 4, control = list(max_treedepth = 11, adapt_delta = 0.95)))) %>% 
  collect() %>% 
  as_tibble()

# regression analyses for lipoprotein lipid composition
lipids.stan <- lipoprotein_lipids %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(stan_output = map(stan_dat, ~sampling(AK, data = .x, chains = 4, control = list(max_treedepth = 11, adapt_delta = 0.95)))) %>% 
  collect() %>% 
  as_tibble()

# regression analyses for amino acids
AAs.stan <- AAs %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(stan_output = map(stan_dat, ~sampling(AK, data = .x, chains = 4, control = list(max_treedepth = 11, adapt_delta = 0.95)))) %>% 
  collect() %>% 
  as_tibble()

# regression analyses for ketone bodies and glycolysis related metabolites
polar.stan <- polar_metabolites %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(stan_output = map(stan_dat, ~sampling(AK, data = .x, chains = 4, control = list(max_treedepth = 11, adapt_delta = 0.95)))) %>% 
  collect() %>% 
  as_tibble()

# regression analyses for fatty acids
ffas.stan <- ffas %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(stan_output = map(stan_dat, ~sampling(AK, data = .x, chains = 4, control = list(max_treedepth = 11, adapt_delta = 0.95)))) %>% 
  collect() %>% 
  as_tibble()

#write_rds(lipoproteins.stan, "stan_output_lipoproteins.rds")

############################################
# 3. Visualizing outcomes of regression analyses
############################################

# Every metabolite will be represented by its
# posterior mean + 95% credible interval
anitschkow_posterior <- function(x, which_par) {
  x %>% 
    mutate(tidy_stan_output = map(stan_output, ~tidy(.x, pars = which_par, estimate.method = "mean", conf.int = TRUE, conf.level = 0.95, conf.method = "HPDinterval"))) %>% 
    unnest(tidy_stan_output)
}

# visualization for lipoproteins (either particle concentration or lipid compositions)
anitschkow_lipoprotein <- function(df, what_name = "Mean difference between treatment and control with 95% credible interval") {
  df %>% 
    mutate(what = factor(what, c("XXL-VLDL-P", "XL-VLDL-P", 'L-VLDL-P', "M-VLDL-P", "S-VLDL-P", "XS-VLDL-P", 
                                 "IDL-P", "L-LDL-P", 'M-LDL-P', "S-LDL-P", 
                                 "XL-HDL-P", "L-HDL-P", "M-HDL-P", "S-HDL-P"))) %>% 
    ggplot(aes(what, estimate)) + 
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high), shape = 1) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
    labs(x = NULL, y = what_name) + 
    coord_flip() + 
    theme_bw()
}

# Figure S2: lipoprotein lipid composition
post.delta <- anitschkow_posterior(lipids.stan, which_par = "delta") %>% 
  select(-data, -stan_dat, -stan_output) %>% 
  mutate(lipoprotein = str_sub(what, 1, -4), 
         lipid = str_sub(what, -2, -1), 
         lipoprotein = factor(lipoprotein, levels = my_lipoproteins))

ggplot(post.delta, aes(lipoprotein, estimate)) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), shape = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  facet_wrap(~lipid, ncol = 4) + 
  labs(x = NULL, y = "Mean difference between treatment and control with 95% credible interval") + 
  coord_flip() + 
  theme_bw()

ggsave("Treatment_effect_lipoprotein_lipid.png", width = 180, height = 180, units = "mm", dpi = 300)

# Figure S3: Fatty acids 
ffa.df <- tibble(
  what = c("DHA", "PUFA", "SFA", "LA", "MUFA", "FAw6", "FAw3", "TotFA"), 
  metabolite = c("Docosahexaenoic acid", 
                 "Polyunsaturated fatty acids", 
                 "Saturated fatty acids", 
                 "Linoleic acid", 
                 "Monounsaturated fatty acids", 
                 "Omega-6", 
                 "Omega-3", 
                 "Total fatty acids")
)

anitschkow_posterior(ffas.stan, which_par = "delta") %>% 
  select(-data, -stan_dat, -stan_output) %>% 
  left_join(ffa.df, by = "what") %>% 
  ggplot(aes(metabolite, estimate)) + 
  geom_pointrange(shape = 1, aes(ymin = conf.low, ymax = conf.high)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = NULL, y = "Mean difference between treatment and control with 95% credible interval") + 
  coord_flip() + 
  theme_bw()

ggsave("Treatment_effect_FFAs.png", width = 180, height = 100, units = "mm", dpi = 300)

# Figure S5: TG in lipoproteins
Nightingale %>% 
  filter(what %in% c("XXL-VLDL-TG", "XL-VLDL-TG", 'L-VLDL-TG', "M-VLDL-TG", "S-VLDL-TG", "XS-VLDL-TG", 
                     "IDL-TG", "L-LDL-TG", 'M-LDL-TG', "S-LDL-TG", 
                     "XL-HDL-TG", "L-HDL-TG", "M-HDL-TG", "S-HDL-TG")) %>% 
  left_join(sample_visit, by = "sampleid") %>% 
  mutate(what = fct_inorder(what), Treatment = factor(Treatment, 
                                                      levels = c("PLACEBO", "AMG145"), 
                                                      labels = c("Placebo", "Evolocumab")), 
         visit = factor(visit, 
                        levels = c("V1", "V2"), 
                        labels = c("Baseline", "16 Weeks after the treatment"))) %>% 
  ggplot(aes(what, y)) + 
  geom_point(aes(color = Treatment), position = position_jitterdodge()) + 
  scale_color_brewer(palette = "Set2") + 
  facet_wrap(~visit, ncol = 1) + 
  labs(x = NULL, y = "Triglycerides (mmol/l)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "bottom", 
        legend.title = element_blank())

ggsave("Triglycerides_in_lipoproteins.png", width = 240, height = 240, units = "mm", dpi = 300)

# Figure S6: lipoprotein particle concentrations
Nightingale %>% 
  filter(what %in% c("XXL-VLDL-P", "XL-VLDL-P", 'L-VLDL-P', "M-VLDL-P", "S-VLDL-P", "XS-VLDL-P")) %>% 
  left_join(sample_visit, by = "sampleid") %>% 
  mutate(what = fct_inorder(what), Treatment = factor(Treatment, 
                                                      levels = c("PLACEBO", "AMG145"), 
                                                      labels = c("Placebo", "Evolocumab")), 
         visit = factor(visit, 
                        levels = c("V1", "V2"), 
                        labels = c("Baseline", "16 Weeks after the treatment"))) %>% 
  ggplot(aes(what, y)) + 
  geom_point(aes(color = Treatment), position = position_jitterdodge()) + 
  scale_color_brewer(palette = "Set2") + 
  facet_wrap(~visit, ncol = 1) + 
  labs(x = NULL, y = "Particle concentration (mol/l)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "bottom", 
        legend.title = element_blank())

ggsave("VLDL_particle_concentrations.png", width = 180, height = 240, units = "mm", dpi = 300)


