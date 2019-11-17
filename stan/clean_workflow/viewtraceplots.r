# Check which chains got stuck on the models that didn't converge.

fp <- '~/forestlight/stanoutput'

library(rstan)
library(bayesplot)

d3_pars <- c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high')
p1_pars <- c('beta0', 'beta1')
p2_pars <- c('beta0', 'beta1_low', 'beta1_high', 'x0', 'delta')

prefixes <- c('fit_density3_production_fg3_1995_', 'fit_density3_production_alltree_1995_', 'fit_production1_volumescaling_unclassified_1995_', 'fit_production2_diamgrowthscaling_fg4_1995_')

files <- lapply(prefixes, paste0, 1:3, '.csv')

fits <- as.list(rep(NA, 4))

fits[[1]] <- read_stan_csv(file.path(fp, files[[1]]))
fits[[2]] <- read_stan_csv(file.path(fp, files[[2]]))
fits[[3]] <- read_stan_csv(file.path(fp, files[[3]]))
fits[[4]] <- read_stan_csv(file.path(fp, files[[4]]))



mcmc_trace(as.array(fits[[1]]), pars = d3_pars)
mcmc_trace(as.array(fits[[2]]), pars = d3_pars)
mcmc_trace(as.array(fits[[3]]), pars = p1_pars)
mcmc_trace(as.array(fits[[4]]), pars = p2_pars)
