# Test ratio fit
# Sample down to 10000.

library(tidyverse)
library(forestscaling)
library(rstan)

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.R'))

options(mc.cores = 2)
rstan_options(auto_write = TRUE)

# Stan model
p1mod <- stan_model('~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/production1.stan')
d3mod <- stan_model('~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/density3_decreasingslopes.stan')

# Weibull model for light per area
weibull_model <- '
  data {
  	int<lower=0> N;
    vector<lower=0>[N] x;
    real<lower=0> UL; // Lower truncation limit
    real<lower=0> LL; // Upper truncation limit
  }
  
  parameters {
    real<lower=0> m;
    real<lower=0> n;
  }
  
  model {
    // Priors: Weibull density
    m ~ lognormal(1, 1);
    n ~ lognormal(1, 1);
    
    // Likelihood: Weibull density
    for (i in 1:N) x[i] ~ weibull(m, n) T[LL,UL];
  }'

wmod <- stan_model(model_code = weibull_model)

# Data for fg's 1-4
N <- 15000 # for FG 3.

set.seed(77)

#range(alltreedat[[3]]$light_received_byarea, na.rm = TRUE)

dat_fg <- map(1:4, ~ with(alltreedat[[3]] %>% filter(fg == .x) %>% sample_n(min(nrow(.), N)), list(x = dbh_corr, y = production, N = length(dbh_corr), x_min = min(dbh_corr), x_max = max(dbh_corr))))
dat_fg_lightarea <- map(1:4, ~ with(alltreedat[[3]] %>% filter(fg == .x, !is.na(light_received_byarea)) %>% sample_n(min(nrow(.), N)), list(x = light_received_byarea, y = production, N = length(dbh_corr), LL = 1, UL = 413, x_min = min(light_received_byarea), x_max = max(light_received_byarea))))

# Fits
density_fit_fg <- map2(dat_fg, 101:104, ~ sampling(d3mod, data = .x, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = .y, pars = 'log_lik', include = FALSE))
production_fit_fg <- map2(dat_fg, 201:204, ~ sampling(p1mod, data = .x, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = .y, pars = 'log_lik', include = FALSE))
density_fit_fg_lightarea <- map2(dat_fg_lightarea, 401:404, ~ sampling(wmod, data = .x, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = .y))
production_fit_fg_lightarea <- map2(dat_fg_lightarea, 701:704, ~ sampling(p1mod, data = .x, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = .y, pars = 'log_lik', include = FALSE))

save(density_fit_fg, production_fit_fg, density_fit_fg_lightarea, production_fit_fg_lightarea, file = '~/Dropbox/Q/projects/forestlight/fitsforratio.RData')

# check convergence
summ_dens <- map(density_fit_fg, summary)
summ_prod <- map(production_fit_fg, summary)
summ_dens_la <- map(density_fit_fg_lightarea, summary)
summ_prod_la <- map(production_fit_fg_lightarea, summary)

# Generate predicted values for each fit
dbh_pred <- logseq(1,286,101)
la_pred <- logseq(1,412,101)
parnames_dens <- c('tau_low', 'tau_high', 'alpha_low', 'alpha_mid', 'alpha_high')
parnames_prod <- c('beta0', 'beta1')

dens_pred_fg <- map(density_fit_fg, function(fit) {
  pars_fg <- extract(fit, parnames_dens) %>% bind_cols
  pmap(pars_fg, pdf_3part, x = dbh_pred, xmin = 1)
})
prod_pred_fg <- map(production_fit_fg, function(fit) {
  pars_fg <- extract(fit, parnames_prod) %>% bind_cols
  pmap(pars_fg, powerlaw_log, x = dbh_pred)
})
dens_pred_fg_lightarea <- map(density_fit_fg_lightarea, function(fit) {
  pars_fg <- extract(fit, c('m','n')) %>% bind_cols
  pmap(pars_fg, function(m, n) truncdist::dtrunc(x = la_pred, spec = 'weibull', a = 1, b = 413, shape = m, scale = n))
})
prod_pred_fg_lightarea <- map(production_fit_fg_lightarea, function(fit) {
  pars_fg <- extract(fit, parnames_prod) %>% bind_cols
  pmap(pars_fg, powerlaw_log, x = la_pred)
})

# Get total production including correction factors.

corr_factor <- function(y, y_fit, n_pars) {
  y_fit <- do.call(cbind, y_fit)
  # Get residuals by subtracting log y from linear predictor
  resids <- -1 * sweep(log(y_fit), 2, log(y))
  
  # Calculate variances and ratio
  pred_var <- apply(log(y_fit), 1, var)
  resid_var <- apply(resids, 1, var)
  r2s <- pred_var / (pred_var + resid_var)
  
  # Sum of squared residuals
  ssq_resid <- apply(resids^2, 1, sum)
  # Standard error of estimates
  sse <- (ssq_resid / (length(y) - n_pars))^0.5
  # Correction factors
  exp((sse^2)/2)
}

# We need all fitted values for production, not just the 101 values, to calculate correction factor.
prod_pred_all <- map2(dat_fg, production_fit_fg, function(dat, fit) {
  pars_fg <- extract(fit, parnames_prod) %>% bind_cols
  pmap(pars_fg, powerlaw_log, x = dat$x)
})
prod_pred_all_lightarea <- map2(dat_fg_lightarea, production_fit_fg_lightarea, function(dat, fit) {
  pars_fg <- extract(fit, parnames_prod) %>% bind_cols
  pmap(pars_fg, powerlaw_log, x = dat$x)
})

prod_cf_fg <- map2(dat_fg, prod_pred_all, ~ corr_factor(y = .x$y, y_fit = .y, n_pars = 2))
prod_cf_fg_lightarea <- map2(dat_fg_lightarea, prod_pred_all_lightarea, ~ corr_factor(y = .x$y, y_fit = .y, n_pars = 2))

# Multiply density x production then multiply by correction factor.
totalprod_pred_fg <- map2(dens_pred_fg, prod_pred_fg, ~ do.call(cbind, .x) * do.call(cbind, .y)) %>%
  map2(prod_cf_fg, ~ sweep(.x, 2, .y, `*`))
totalprod_pred_fg_lightarea <- map2(dens_pred_fg_lightarea, prod_pred_fg_lightarea, ~ do.call(cbind, .x) * do.call(cbind, .y)) %>%
  map2(prod_cf_fg_lightarea, ~ sweep(.x, 2, .y, `*`))

# Take the ratio of the predicted values from each sampling iteration
dens_ratio13 <- map2(dens_pred_fg[[1]], dens_pred_fg[[3]], `/`)
dens_ratio24 <- map2(dens_pred_fg[[2]], dens_pred_fg[[4]], `/`)
dens_ratio13_lightarea <- map2(dens_pred_fg_lightarea[[1]], dens_pred_fg_lightarea[[3]], `/`)
dens_ratio24_lightarea <- map2(dens_pred_fg_lightarea[[2]], dens_pred_fg_lightarea[[4]], `/`)

totalprod_ratio13 <- totalprod_pred_fg[[1]] / totalprod_pred_fg[[3]]
totalprod_ratio24 <- totalprod_pred_fg[[2]] / totalprod_pred_fg[[4]]
totalprod_ratio13_lightarea <- totalprod_pred_fg_lightarea[[1]] / totalprod_pred_fg_lightarea[[3]]
totalprod_ratio24_lightarea <- totalprod_pred_fg_lightarea[[2]] / totalprod_pred_fg_lightarea[[4]]

# Generate credible intervals from the ratios
# Correct for the "subsetting" of fg3
subset_corr_factor <- N/sum(alltreedat[[3]]$fg %in% 3)

# By diameter
dens_ratio_ci13 <- do.call(cbind, dens_ratio13) %>%
  `*`(subset_corr_factor) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(dbh = dbh_pred) 

dens_ratio_ci24 <- do.call(cbind, dens_ratio24) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(dbh = dbh_pred)

totalprod_ratio_ci13 <- totalprod_ratio13 %>%
  `*`(subset_corr_factor) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(dbh = dbh_pred) 

totalprod_ratio_ci24 <- totalprod_ratio24 %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(dbh = dbh_pred)


# By light per area
subset_corr_factor_la <- N/sum(alltreedat[[3]]$fg %in% 3 & !is.na(alltreedat[[3]]$light_received_byarea))

dens_ratio_ci13_lightarea <- do.call(cbind, dens_ratio13_lightarea) %>%
  `*`(subset_corr_factor_la) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(dbh = dbh_pred) 

dens_ratio_ci24_lightarea <- do.call(cbind, dens_ratio24_lightarea) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(dbh = dbh_pred)

totalprod_ratio_ci13_lightarea <- totalprod_ratio13_lightarea %>%
  `*`(subset_corr_factor_la) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(dbh = dbh_pred) 

totalprod_ratio_ci24_lightarea <- totalprod_ratio24_lightarea %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(dbh = dbh_pred)

# Load observed ratio
obs_ratio13 <- read_csv(file.path(gdrive_path, 'data/data_binned/fastslow_stats_bydiam_byyear.csv')) %>%
  filter(n_individuals > 10, year == 1995)
obs_ratio24 <- read_csv(file.path(gdrive_path, 'data/data_binned/breeder_stats_bydiam_byyear.csv')) %>%
  filter(n_individuals > 10, year == 1995)
obs_ratio13_lightarea <- read_csv(file.path(gdrive_path, 'data/data_binned/fastslow_stats_bylight_byyear.csv')) %>%
  filter(n_individuals > 10, year == 1995)
obs_ratio24_lightarea <- read_csv(file.path(gdrive_path, 'data/data_binned/breeder_stats_bylight_byyear.csv')) %>%
  filter(n_individuals > 10, year == 1995)


# Plot predicted and observed ratio 1:3 and 2:4: density & production x diameter & light

ggplot() +
  geom_ribbon(data = dens_ratio_ci24, aes(x = dbh, ymin = low, ymax = hi), alpha = 0.5) +
  geom_ribbon(data = dens_ratio_ci13 %>% filter(between(dbh, 1, 100)), aes(x = dbh, ymin = low, ymax = hi), alpha = 0.5) +
  geom_line(data = dens_ratio_ci24, aes(x = dbh, y = med)) +
  geom_line(data = dens_ratio_ci13 %>% filter(between(dbh, 1, 100)), aes(x = dbh, y = med)) +
  geom_point(data = obs_ratio24, aes(x = bin_midpoint, y = breeder_density_ratio), size = 4, pch = 21) +
  geom_point(data = obs_ratio13, aes(x = bin_midpoint, y = fastslow_density_ratio), size = 4, pch = 21) +
  scale_x_log10(limits = c(1, 100)) + scale_y_log10(name = 'abundance ratio', limits = c(0.01, 100)) +
  theme_minimal()
  
ggplot() +
  geom_ribbon(data = totalprod_ratio_ci24, aes(x = dbh, ymin = low, ymax = hi), alpha = 0.5) +
  geom_ribbon(data = totalprod_ratio_ci13 %>% filter(between(dbh, 1, 100)), aes(x = dbh, ymin = low, ymax = hi), alpha = 0.5) +
  geom_line(data = totalprod_ratio_ci24, aes(x = dbh, y = med)) +
  geom_line(data = totalprod_ratio_ci13 %>% filter(between(dbh, 1, 100)), aes(x = dbh, y = med)) +
  geom_point(data = obs_ratio24, aes(x = bin_midpoint, y = breeder_production_ratio), size = 4, pch = 21) +
  geom_point(data = obs_ratio13, aes(x = bin_midpoint, y = fastslow_production_ratio), size = 4, pch = 21) +
  scale_x_log10(limits = c(1, 100)) + scale_y_log10(name = 'production ratio', limits = c(0.01, 100)) +
  theme_minimal()

ggplot() +
  geom_ribbon(data = dens_ratio_ci24_lightarea, aes(x = dbh, ymin = low, ymax = hi), alpha = 0.5) +
  geom_ribbon(data = dens_ratio_ci13_lightarea, aes(x = dbh, ymin = low, ymax = hi), alpha = 0.5) +
  geom_line(data = dens_ratio_ci24_lightarea, aes(x = dbh, y = med)) +
  geom_line(data = dens_ratio_ci13_lightarea, aes(x = dbh, y = med)) +
  geom_point(data = obs_ratio24_lightarea, aes(x = bin_midpoint, y = breeder_density_ratio), size = 4, pch = 21) +
  geom_point(data = obs_ratio13_lightarea, aes(x = bin_midpoint, y = fastslow_density_ratio), size = 4, pch = 21) +
  scale_x_log10(limits = c(1, 413)) + scale_y_log10(name = 'abundance ratio', limits = c(0.01, 100)) +
  theme_minimal()

ggplot() +
  geom_ribbon(data = totalprod_ratio_ci24_lightarea, aes(x = dbh, ymin = low, ymax = hi), alpha = 0.5) +
  geom_ribbon(data = totalprod_ratio_ci13_lightarea %>% filter(between(dbh, 1, 100)), aes(x = dbh, ymin = low, ymax = hi), alpha = 0.5) +
  geom_line(data = totalprod_ratio_ci24_lightarea, aes(x = dbh, y = med)) +
  geom_line(data = totalprod_ratio_ci13_lightarea %>% filter(between(dbh, 1, 100)), aes(x = dbh, y = med)) +
  geom_point(data = obs_ratio24, aes(x = bin_midpoint, y = breeder_production_ratio), size = 4, pch = 21) +
  geom_point(data = obs_ratio13, aes(x = bin_midpoint, y = fastslow_production_ratio), size = 4, pch = 21) +
  scale_x_log10(limits = c(1, 413)) + scale_y_log10(name = 'production ratio', limits = c(0.01, 100)) +
  theme_minimal()
