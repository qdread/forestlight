# Test new density 3 fit

library(tidyverse)
library(forestscaling)
library(rstan)
library(bayesplot)

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

options(mc.cores = 2)
rstan_options(auto_write = TRUE)

# Stan model
oldmod <- stan_model('~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/density3.stan')
newmod <- stan_model('~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/density3_decreasingslopes.stan')

d3_pars <- c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high')

# Data objects
dat_fg3 <- with(alltreedat[[3]] %>% filter(fg == 3), list(x = dbh_corr, y = production, N = length(dbh_corr), x_min = min(dbh_corr), x_max = max(dbh_corr)))

N <- 5000
set.seed(333)
dat_subsample <- with(alltreedat[[3]] %>% filter(fg == 3) %>% sample_n(N), list(x = dbh_corr, y = production, N = length(dbh_corr), x_min = min(dbh_corr), x_max = max(dbh_corr)))

str(dat_subsample)

# Fit the model!
fit_sub <- sampling(newmod, data = dat_subsample, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222, pars = 'log_lik', include = FALSE)
#fit_sub_old <- sampling(oldmod, data = dat_subsample, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222)

summary(fit_sub, pars = d3_pars)
mcmc_trace(as.array(fit_sub), pars = d3_pars)
