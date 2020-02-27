# Fit models for 3-Seg density & 1-Seg production for each year.
# Get Bayesian estimates, but with subsample of data so that it fits quickly.

library(tidyverse)
library(rstan)
options(mc.cores = 2)
rstan_options(auto_write = TRUE)

# Load data

load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r') # doesn't include imputed values

# Compile models

mod_dens3 <- stan_model('stan/clean_workflow/model_scripts/density3.stan')
mod_prod1 <- stan_model('stan/clean_workflow/model_scripts/production1.stan')

# Create data objects
N <- 10000
set.seed(111)

get_stan_data <- function(dat) with(dat, list(N = nrow(dat), x = dat$dbh_corr, y = dat$production, x_min = 1, x_max = max(dat$dbh_corr)))

allstandat_byfg <- map(alltreedat, function(x) {
  x %>%
    filter(!recruit) %>%
    mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
    filter(!is.na(dbh_corr)) %>%
    group_by(fg) %>%
    group_map(~ get_stan_data(if (nrow(.) <= N) . else sample_n(., N)))
})

allstandat_alltrees <- map(alltreedat, function(x) {
  x %>%
    filter(!recruit) %>%
    mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
    filter(!is.na(dbh_corr)) %>%
    sample_n(N) %>%
    get_stan_data()
})

# Fit with Bayesian method

n_chains <- 2
n_iter <- 4000
n_warmup <- 3000

prodfit_alltrees <- map(allstandat_alltrees[-1], ~ sampling(mod_prod1, data = ., pars = 'log_lik', include = FALSE, chains = n_chains, iter = n_iter, warmup = n_warmup))
densfit_alltrees <- map(allstandat_alltrees[-1], ~ sampling(mod_dens3, data = ., pars = 'log_lik', include = FALSE, chains = n_chains, iter = n_iter, warmup = n_warmup))

prodfit_byfg <- map(allstandat_byfg[-1], function(yeardat) map(yeardat, ~ sampling(mod_prod1, data = ., pars = 'log_lik', include = FALSE, chains = n_chains, iter = n_iter, warmup = n_warmup)))
densfit_byfg <- map(allstandat_byfg[-1], function(yeardat) map(yeardat, ~ sampling(mod_dens3, data = ., pars = 'log_lik', include = FALSE, chains = n_chains, iter = n_iter, warmup = n_warmup)))

# Extract slopes from model fits (point estimates only)
get_pars <- function(fit) {
  summ <- summary(fit)$summary
  summ[-nrow(summ), "50%"]
}

prodpars_alltrees <- map(prodfit_alltrees, get_pars)
denspars_alltrees <- map(densfit_alltrees, get_pars)

prodpars_byfg <- map(prodfit_byfg, function(yearfit) map(yearfit, get_pars))
denspars_byfg <- map(densfit_byfg, function(yearfit) map(yearfit, get_pars))

# Combine production and density slopes into a single data frame.
ys <- seq(1990, 2010, 5)
fgs <- c(paste0('fg', 1:5), 'unclassified')
prodpars_alltrees_df <- map2_dfr(prodpars_alltrees, ys, ~ data.frame(year = .y, fg = 'alltree', parameter = names(.x), value = .x))
denspars_alltrees_df <- map2_dfr(denspars_alltrees, ys, ~ data.frame(year = .y, fg = 'alltree', parameter = names(.x), value = .x))

prodpars_byfg_df <- map2_dfr(prodpars_byfg, ys, ~ data.frame(year = .y, map2_dfr(.x, fgs, function(dat, fg) data.frame(fg = fg, parameter = names(dat), value = dat))))
denspars_byfg_df <- map2_dfr(denspars_byfg, ys, ~ data.frame(year = .y, map2_dfr(.x, fgs, function(dat, fg) data.frame(fg = fg, parameter = names(dat), value = dat))))

# Write outputs temporarily so that models do not need to be run again
write_csv(prodpars_byfg_df, '~/google_drive/ForestLight/data/data_piecewisefits/allyearfits_prodpars.csv')
write_csv(denspars_byfg_df, '~/google_drive/ForestLight/data/data_piecewisefits/allyearfits_denspars.csv')

# Get total production slopes using cutoffs from the data.