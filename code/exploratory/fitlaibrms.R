library(tidyverse)
library(brms)

lai_depth <- read_csv('~/google_drive/ForestLight/data/data_forplotting/LAI_Depth (1).csv') %>%
  mutate(log_LAI = log10(LAI),
         log_depth = log10(Depth)) %>%
  filter(!Species %in% 'Cecropia')

pfd_lai <- read_csv('~/google_drive/ForestLight/data/data_forplotting/PFD_LAI.csv') %>%
  mutate(log_PFD = log(PFD/100)) %>%
  filter(!Species %in% 'Cecropia')

options(mc.cores = 3)

# Fit Bayesian random-intercept, random-slope model
lai_fit <- brm(log_LAI ~ log_depth + (log_depth | Species), data = lai_depth, family = 'gaussian', chains = 3, iter = 10000, warmup = 9000, seed = 33)

summary(lai_fit)
bayes_R2(lai_fit)

pfd_fit <- brm(log_PFD ~ LAI + (LAI | Species), data = pfd_lai, family = 'gaussian', chains = 3, iter = 20000, warmup = 18000, seed = 66)

summary(pfd_fit)
bayes_R2(pfd_fit)

curve(1 - exp(-0.5 * x), from = 0.1, to = 15, log = 'x', ylim = c(0,1))
curve(1 - exp(-0.5 * (1.58 * x ^ 0.35)), from = 0.1, to = 15, log = 'x', ylim = c(0,1), col = 'red', add  =TRUE)
curve(1 - exp(-0.5 * (2.29 * x ^ 0.37)), from = 0.1, to = 15, log = 'x', col = 'green', add = TRUE)

library(lme4)
lmer(log(PFD/100) ~ LAI + (LAI | Species), data = pfd_lai)
lm(log_PFD ~ LAI, data = pfd_lai)
