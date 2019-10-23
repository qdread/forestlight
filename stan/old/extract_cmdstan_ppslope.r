# For each CMDstan model fit, get the following
# Credible intervals of parameters
# Predictive intervals to make graphs
# Information criteria

# Edit 12 April: Do this in parallel because it is too slow
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

extract_all_fit <- function(dens_model, prod_model, fg, year, xmin, n, total_production) {
  require(rstan)

  # Load CSVs as stanfit object
  fp <- '~/forestlight/stanoutput'
  files <- paste0('midsizefit_slope_', dens_model, 'x', prod_model, '_', fg, '_', year, '_', 1:3, '.csv')
  fit <- read_stan_csv(file.path(fp,files))
  
  prod_model <- ifelse(prod_model == 'power', 'powerlaw', 'powerlawexp')
  
  # Get credible intervals of parameters
  pareto_par <- c('alpha')
  weibull_par <- c('m', 'n')
  powerlaw_par <- c('beta0', 'beta1', 'ee_slope')
  powerlawexp_par <- c('beta0', 'beta1', 'a', 'b', 'c')
  
  if (dens_model == 'pareto') get_pars <- pareto_par else get_pars <- weibull_par
  if (prod_model == 'powerlaw') get_pars <- c(get_pars, powerlaw_par) else get_pars <- c(get_pars, powerlawexp_par)

  summ_fit <- summary(fit)
  param_cis <- cbind(data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg), summ_fit[[1]][get_pars, ])

  param_cis <- cbind(parameter = dimnames(param_cis)[[1]], param_cis)
  names(param_cis)[9:13] <- c('q025', 'q25', 'q50', 'q75', 'q975')

  # Get predictive intervals
  pred_interval <- dens_prod_ci(fit, dbh_pred, dens_form = dens_model, prod_form = prod_model, total_prod = total_production, x_min = xmin, n_indiv = n)
  pred_interval <- data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg, pred_interval)
  
  list(param_cis = param_cis, 
       pred_interval = pred_interval)
}

source('~/forestlight/stancode/extract_ci_stan.r')

library(purrr)
library(dplyr)

mod_df <- expand.grid(dens_model = c('pareto'),
                      prod_model = c('power'),
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = seq(1990, 2010, 5), 
                      stringsAsFactors = FALSE)

min_n <- read.csv('~/forestlight/stancode/min_n_midsize.csv', stringsAsFactors = FALSE)

total_prod <- read.csv('~/forestlight/production_total_midsize.csv', stringsAsFactors = FALSE)

mod_df <- mod_df %>%
  left_join(total_prod) %>%
  rename(total_production = production) %>%
  left_join(min_n)

dbh_pred <- exp(seq(log(1.2), log(315), length.out = 101))

fit_info <- extract_all_fit(dens_model = mod_df$dens_model[task],
                            prod_model = mod_df$prod_model[task],
                            fg = mod_df$fg[task],
                            year = mod_df$year[task],
                            xmin = mod_df$xmin[task],
                            n = mod_df$n[task],
                            total_production = mod_df$total_production[task])

save(fit_info, file = paste0('~/forestlight/stanoutput/fitinfo/midsizeinfo_slope_',task,'.r'))