# Fit models for 3-Seg density & 1-Seg production for each year.
# Just get the point estimates so that it works more quickly.

library(tidyverse)
library(rstan)

# Load data

load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r') # doesn't include imputed values

# Compile models

mod_dens3 <- stan_model('stan/clean_workflow/model_scripts/density3.stan')
mod_prod1 <- stan_model('stan/clean_workflow/model_scripts/production1.stan')

# Create data objects

get_stan_data <- function(dat) with(dat, list(N = nrow(dat), x = dat$dbh_corr, y = dat$production, x_min = 1, x_max = max(dat$dbh_corr)))

allstandat_byfg <- map(alltreedat, function(x) {
  x %>%
    filter(!recruit) %>%
    mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
    filter(!is.na(dbh_corr)) %>%
    group_by(fg) %>%
    group_map(~ get_stan_data(.))
})

allstandat_alltrees <- map(alltreedat, function(x) {
  x %>%
    filter(!recruit) %>%
    mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
    filter(!is.na(dbh_corr)) %>%
    get_stan_data()
})


# Fit with MLE
# This will take a few minutes on the local machine but nothing outrageous.
dens_inits <- list(tau_low = 5.0, tau_high = 25.0, alpha_low = 0.1, alpha_mid = 2.0, alpha_high = 5.0)

prodfit_alltrees <- map(allstandat_alltrees[-1], ~ optimizing(mod_prod1, data = ., verbose = TRUE))
densfit_alltrees <- map(allstandat_alltrees[-1], ~ optimizing(mod_dens3, data = ., verbose = TRUE, init = dens_inits))

prodfit_byfg <- map(allstandat_byfg[-1], function(yeardat) map(yeardat, ~ optimizing(mod_prod1, data = ., verbose = TRUE)))
densfit_byfg <- map(allstandat_byfg[-1], function(yeardat) map(yeardat, ~ optimizing(mod_dens3, data = ., verbose = TRUE, init = dens_inits[3:5])))

# Extract slopes from model fits
get_pars <- function(fit) fit$par[!grepl('log_lik', names(fit$par))]

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

# Get total production slopes