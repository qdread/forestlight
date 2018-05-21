# For each CMDstan model fit, get the fitted values of log slope at all predicted value points, and their credible intervals

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

extract_logslope <- function(dens_model, prod_model, fg, year, xmin, n, total_production) {
  require(rstan)

  # Load CSVs as stanfit object
  fp <- '~/forestlight/stanoutput'
  files <- paste0('fit_', dens_model, 'x', prod_model, '_', fg, '_', year, '_', 1:3, '.csv')
  if (fg == 'alltree') files <- paste0('ss', files) # Use the 25K subset for all trees.
  fit <- read_stan_csv(file.path(fp,files))
  
  prod_model <- ifelse(prod_model == 'power', 'powerlaw', 'powerlawexp')
  
  # Get fitted slopes
  fitted_slopes <- fitted_slope_ci(fit, dbh_pred, dens_form = dens_model, prod_form = prod_model, total_prod = total_production, x_min = xmin, n_indiv = n)
  data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg, fitted_slopes)
  
}

source('~/forestlight/stancode/get_fitted_slope_quant.r')

library(purrr)
library(dplyr)

mod_df <- expand.grid(dens_model = c('pareto', 'weibull'),
                      prod_model = c('power', 'exp'),
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = seq(1990, 2010, 5), 
                      stringsAsFactors = FALSE)

min_n <- read.csv('~/forestlight/stancode/min_n.csv', stringsAsFactors = FALSE)

total_prod <- read.csv('~/forestlight/production_total.csv', stringsAsFactors = FALSE)

mod_df <- mod_df %>%
  left_join(total_prod) %>%
  rename(total_production = production) %>%
  left_join(min_n)

dbh_pred <- exp(seq(log(1.2), log(315), length.out = 101))

fit_info <- extract_logslope(dens_model = mod_df$dens_model[task],
                            prod_model = mod_df$prod_model[task],
                            fg = mod_df$fg[task],
                            year = mod_df$year[task],
                            xmin = mod_df$xmin[task],
                            n = mod_df$n[task],
                            total_production = mod_df$total_production[task])

save(fit_info, file = paste0('~/forestlight/stanoutput/fitinfo/logslope_',task,'.r'))