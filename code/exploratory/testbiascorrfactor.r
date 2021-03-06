library(tidyverse)
library(rstan)
library(forestscaling)
options(mc.cores = 2)
rstan_options(auto_write = TRUE)

# Code to test Sprugel 1983 correction factor

log_correction_factor <- function(log_y, log_y_pred, n_pars, base = exp(1)) {
  std_err_est <- (sum((log_y - log_y_pred)^2) / (length(log_y) - n_pars)) ^ 0.5
  base ^ ((std_err_est^2)/2)
}

# Load some data to test
load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')

# Compile stan models
mod_p1 <- stan_model(file = '~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/production1.stan', model_name = 'p1segment')
mod_p2 <- stan_model(file = '~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/production2.stan', model_name = 'p2segment')
mod_d2 <- stan_model(file = '~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/density2.stan', model_name = 'd2segment')
mod_d3 <- stan_model(file = '~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/density3.stan', model_name = 'd3segment')

# Make a toy dataset to fit the model
set.seed(4)
N <- 10000
idx <- sample(nrow(alltreedat[[3]]), N)
dat <- with(alltreedat[[3]][idx, ], list(x = dbh_corr, y = production, N = N, x_min = 1, x_max = 286))

# Fit both models
fit_p1 <- sampling(mod_p1, dat, pars = 'log_lik', include = FALSE, chains = 2, iter = 2000, warmup = 1000, seed = 99)
fit_p2 <- sampling(mod_p2, dat, pars = 'log_lik', include = FALSE, chains = 2, iter = 2000, warmup = 1000, seed = 99)
fit_d2 <- sampling(mod_d2, dat, pars = 'log_lik', include = FALSE, chains = 2, iter = 2000, warmup = 1000, seed = 99)
fit_d3 <- sampling(mod_d3, dat, pars = 'log_lik', include = FALSE, chains = 2, iter = 2000, warmup = 1000, seed = 99)

# Check convergence
summary(fit_p1)$summary
summary(fit_p2)$summary
summary(fit_d2)$summary
summary(fit_d3)$summary

###### Get predicted values from all iterations and use this to calculate the correction factor
# extract pars
pars_p1 <- extract(fit_p1, c('beta0', 'beta1'))
pars_p2 <- extract(fit_p2, c('x0', 'beta0', 'beta1_low', 'beta1_high', 'delta'))
pars_d2 <- extract(fit_d2, c('alpha_low', 'alpha_high', 'tau'))
pars_d3 <- extract(fit_d3, c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high'))

hinge_fn <- function(x = dat$x, x0, beta0, beta1_low, beta1_high, delta) {
  xdiff <- log(x) - log(x0)
  exp( log(beta0) + beta1_low * xdiff + (beta1_high - beta1_low) * delta * log(1 + exp(Brobdingnag::as.brob(xdiff / delta))) )
}

# use pars to get fitted values
fittedval_p1 <- mapply(function(x = dat$x, beta0, beta1) beta0 * x^beta1, beta0 = pars_p1$beta0, beta1 = pars_p1$beta1)
fittedval_p2 <- mapply(hinge_fn, x0 = pars_p2$x0, beta0 = pars_p2$beta0, beta1_low = pars_p2$beta1_low, beta1_high = pars_p2$beta1_high, delta = pars_p2$delta)

# Calculate correction factor for each iteration
cf_p1 <- apply(fittedval_p1, 2, function(ypred) log_correction_factor(log_y = log(dat$y), log_y_pred = log(ypred), n_pars = 2))
cf_p2 <- apply(fittedval_p2, 2, function(ypred) log_correction_factor(log_y = log(dat$y), log_y_pred = log(ypred), n_pars = 5))

# Correct the fitted values by multiplying them by the correction factor
correctedval_p1 <- sweep(fittedval_p1, 2, cf_p1, `*`)
correctedval_p2 <- sweep(fittedval_p2, 2, cf_p2, `*`)

# Fitted values at spaced out intervals
dbh_pred <- logseq(1, 286, 50)

hinge_fn <- function(x = dbh_pred, x0, beta0, beta1_low, beta1_high, delta) {
  xdiff <- log(x) - log(x0)
  exp( log(beta0) + beta1_low * xdiff + (beta1_high - beta1_low) * delta * log(1 + exp(Brobdingnag::as.brob(xdiff / delta))) )
}

fitted_p1_plot <- mapply(function(x = dbh_pred, beta0, beta1) beta0 * x^beta1, beta0 = pars_p1$beta0, beta1 = pars_p1$beta1)
fitted_p2_plot <- mapply(hinge_fn, x0 = pars_p2$x0, beta0 = pars_p2$beta0, beta1_low = pars_p2$beta1_low, beta1_high = pars_p2$beta1_high, delta = pars_p2$delta)

pdf_2 <- function(x = dbh_pred, xmin = 1, alpha_low, alpha_high, tau) {
  C_con <- tau ^ (alpha_low - alpha_high)
  C_norm <- ( (C_con / alpha_low) * (xmin ^ (-alpha_low) - tau ^ (-alpha_low)) + ( tau ^ (-alpha_high) ) / alpha_high ) ^ -1
 
  prob <- case_when(
    x < tau ~ C_con * C_norm * ( x ^ - (alpha_low + 1) ),
    x >= tau ~ C_norm * ( x ^ - (alpha_high + 1) )
  )
  return(prob)
}

pdf_3 <- function(x = dbh_pred, xmin = 1, alpha_low, alpha_mid, alpha_high, tau_low, tau_high) {

    C_con_low <- tau_low ^ (alpha_low - alpha_mid)
    C_con_high <- tau_high ^ (alpha_high - alpha_mid)
    C_norm <- ( (C_con_low / alpha_low) * (xmin ^ -alpha_low - tau_low ^ -alpha_low) + (1 / alpha_mid) * (tau_low ^ -alpha_mid - tau_high ^ -alpha_mid) + (C_con_high / alpha_high) * (tau_high ^ -alpha_high) ) ^ -1
    
    prob <- case_when(
      x < tau_low ~ C_con_low * C_norm * ( x ^ - (alpha_low + 1) ),
      x >= tau_low & x <= tau_high ~ C_norm * ( x ^ - (alpha_mid + 1) ),
      x > tau_high ~ C_con_high * C_norm * ( x ^ - (alpha_high + 1) ) 
    )

    return(prob)
}


fitted_d2_plot <- mapply(pdf_2, alpha_low = pars_d2$alpha_low, alpha_high = pars_d2$alpha_high, tau = pars_d2$tau)
fitted_d3_plot <- mapply(pdf_3, alpha_low = pars_d3$alpha_low, alpha_mid = pars_d3$alpha_mid, alpha_high = pars_d3$alpha_high, tau_low = pars_d3$tau_low, tau_high = pars_d3$tau_high)


fitted_totalprod <- fitted_d3_plot * fitted_p1_plot * nrow(alltreedat[[3]])
fitted_totalprod_corrected <- sweep(fitted_totalprod, 2, cf_p1, `*`)

quantile_corrected <- apply(fitted_totalprod_corrected, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t %>% as.data.frame %>% setNames(c('qlow','med','qhi'))
quantile_uncorrected <- apply(fitted_totalprod, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t %>% as.data.frame %>% setNames(c('qlow','med','qhi'))

# Calculate old school correction factor too.
integrals <- apply(fitted_totalprod, 2, function(y) pracma::trapz(x=dbh_pred,y=y))
totalprod <- sum(alltreedat[[3]]$production)
fitted_totalprod_corrected_oldschool <- sweep(fitted_totalprod, 2, totalprod/integrals, `*`)
quantile_corrected_oldschool <- apply(fitted_totalprod_corrected_oldschool, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t %>% as.data.frame %>% setNames(c('qlow','med','qhi'))


quantile_dens <- apply(fitted_d2_plot * nrow(alltreedat[[3]]), 1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t %>% as.data.frame %>% setNames(c('qlow','med','qhi'))
quantile_dens3 <- apply(fitted_d3_plot * nrow(alltreedat[[3]]), 1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t %>% as.data.frame %>% setNames(c('qlow','med','qhi'))

# Load binned data
bintp <- read_csv('~/google_drive/ForestLight/data/data_binned/totalproductionbin_byyear.csv')
all95 <- bintp %>% filter(fg=='all',year==1995)

plot(all95$bin_midpoint, all95$bin_value, log = 'xy')
lines(dbh_pred, quantile_corrected$med, col = 'forestgreen')
lines(dbh_pred, quantile_uncorrected$med, col = 'goldenrod')
lines(dbh_pred, quantile_corrected_oldschool$med, col = 'red')

bindens <- read_csv('~/google_drive/ForestLight/data/data_binned/densitybin_byyear.csv')
all95dens <- bindens %>% filter(fg=='all',year==1995)
plot(all95dens$bin_midpoint, all95dens$bin_value, log = 'xy')
lines(dbh_pred, quantile_dens$med, col = 'red')
lines(dbh_pred, quantile_dens3$med, col = 'blue')
