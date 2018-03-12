# Credible intervals for all fits for all years.

fp <- '~/forestlight/stanoutput'
fnames <- dir(fp, pattern = 'fit_')

source('~/forestlight/stancode/extract_ci_stan.r')
library(rstan)
library(purrr)

z <- expand.grid(year = c(1990, 1995, 2000, 2005, 2010),
                 fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified', 'alltree'),
                 model_name = c('ppow', 'wpow', 'pexp', 'wexp'),
                 stringsAsFactors = FALSE)

min_n <- read.csv('~/forestlight/stancode/min_n.csv', stringsAsFactors = FALSE)
dbh_pred <- seq(1.2, 315, length.out = 50)

all_cis <- pmap(z, function(year, fg, model_name) {
  fn <- file.path(fp, paste0('fit_', model_name, '_', fg, '_', year, '.r'))
  prod_model <- ifelse(model_name %in% c('ppow','wpow'), 'powerlaw', 'powerlawexp')
  dens_model <- ifelse(model_name %in% c('ppow', 'pexp'), 'pareto', 'weibull')
  x_min_i <- min_n$xmin[min_n$year == year & min_n$fg == fg]
  n_i <- min_n$n[min_n$year == year & min_n$fg == fg]
  load(fn)
  ci <- dens_prod_ci(fit, dbh_pred, dens_form = dens_model, prod_form = prod_model, x_min = x_min_i, n_indiv = n_i)
  data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg, ci)
})

ci_df <- do.call(rbind, all_cis)

