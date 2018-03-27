# For each CMDstan model fit, get the following
# Credible intervals of parameters
# Predictive intervals to make graphs
# Information criteria

extract_all_fit <- function(dens_model, prod_model, fg, year, xmin, n) {
  require(rstan)
  require(loo)
  
  
  # Load CSVs as stanfit object
  files <- paste0('~/forestlight/stanoutput/fit_', dens_model, 'x', prod_model, '_', fg, '_', year, '_', 1:3, '.csv')
  fit <- read_stan_csv(files)
  
  prod_model <- ifelse(prod_model == 'power', 'powerlaw', 'powerlawexp')
  
  # Get credible intervals of parameters
  pareto_par <- c('alpha')
  weibull_par <- c('m', 'n')
  powerlaw_par <- c('beta0', 'beta1')
  powerlawexp_par <- c('beta0', 'beta1', 'a', 'b', 'c')
  
  if (dens_model == 'pareto') get_pars <- pareto_par else get_pars <- weibull_par
  if (prod_model == 'powerlaw') get_pars <- c(get_pars, powerlaw_par) else get_pars <- c(get_pars, powerlawexp_par)

  summ_fit <- summary(fit)
  param_cis <- cbind(data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg), summ_fit[[1]][get_pars, ])

  param_cis <- cbind(parameter = dimnames(param_cis)[[1]], param_cis)
  names(param_cis)[9:13] <- c('q025', 'q25', 'q50', 'q75', 'q975')

  # Get predictive intervals
  pred_interval <- dens_prod_ci(fit, dbh_pred, dens_form = dens_model, prod_form = prod_model, x_min = xmin, n_indiv = n)
  pred_interval <- data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg, pred_interval)

  # Get WAIC and LOOIC
  ll_dens <- extract_log_lik(fit, 'log_lik_dens')
  ll_prod <- extract_log_lik(fit, 'log_lik_prod')
  waic_dens <- waic(ll_dens)
  waic_prod <- waic(ll_prod)
  loo_dens <- loo(ll_dens)
  loo_prod <- loo(ll_prod)
  
  list(waic_dens = waic_dens, 
       waic_prod = waic_prod, 
       loo_dens = loo_dens, 
       loo_prod = loo_prod, 
       param_cis = param_cis, 
       pred_interval = pred_interval)
}

source('~/forestlight/stancode/extract_ci_stan.r')

library(purrr)
library(dplyr)

mod_df <- expand.grid(dens_model = c('pareto', 'weibull'),
                      prod_model = c('power', 'exp'),
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = seq(1990, 2010, 5), 
                      stringsAsFactors = FALSE)

min_n <- read.csv('~/forestlight/stancode/min_n.csv', stringsAsFactors = FALSE)

mod_df <- left_join(mod_df, min_n)

dbh_pred <- exp(seq(log(1.2), log(315), length.out = 50))

fit_info <- pmap(mod_df, extract_all_fit)