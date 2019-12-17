# Test ratio fit

library(tidyverse)
library(forestscaling)
library(rstan)

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.R'))

options(mc.cores = 2)
rstan_options(auto_write = TRUE)

# Stan model
d1mod <- stan_model('~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/density1.stan')
d2mod <- stan_model('~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/density2.stan')
d3mod <- stan_model('~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/density3_decreasingslopes.stan')

# Data for fg's 1-4

dat_fg <- map(1:4, ~ with(alltreedat[[3]] %>% filter(fg == .x), list(x = dbh_corr, y = production, N = length(dbh_corr), x_min = min(dbh_corr), x_max = max(dbh_corr))))

fit_fg1 <- sampling(d3mod, data = dat_fg[[1]], chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222, pars = 'log_lik', include = FALSE)
fit_fg2 <- sampling(d3mod, data = dat_fg[[2]], chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222, pars = 'log_lik', include = FALSE)
fit_fg3 <- sampling(d3mod, data = dat_fg[[3]], chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222, pars = 'log_lik', include = FALSE)
fit_fg4 <- sampling(d3mod, data = dat_fg[[4]], chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222, pars = 'log_lik', include = FALSE)

fit1_fg2 <- sampling(d1mod, data = dat_fg2, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222)
fit1_fg4 <- sampling(d1mod, data = dat_fg4, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222)

# Generate predicted values for each fit
dbh_pred <- logseq(1,286,101)
parnames <- c('tau_low', 'tau_high', 'alpha_low', 'alpha_mid', 'alpha_high')

pred_fg <- map(list(fit_fg1, fit_fg2, fit_fg3, fit_fg4), function(fit) {
  pars_fg <- extract(fit, parnames) %>% bind_cols
  pmap(pars_fg, pdf_3part, x = dbh_pred, xmin = 1)
})

# Take the ratio of the predicted values from each sampling iteration
pred_ratio13 <- map2(pred_fg[[1]], pred_fg[[3]], '/')
pred_ratio24 <- map2(pred_fg[[2]], pred_fg[[4]], `/`)

# Generate credible intervals from the ratios
pred_ratio_ci13 <- do.call(cbind, pred_ratio13) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(dbh = dbh_pred) 

pred_ratio_ci24 <- do.call(cbind, pred_ratio24) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(dbh = dbh_pred) %>%
  filter(dbh <= 25)

# Load observed ratio
obs_ratio <- read_csv(file.path(gdrive_path, 'data/data_binned/breeder_stats_bydiam_byyear.csv')) %>%
  filter(n_individuals > 10, year == 1995)

# Plot predicted and observed ratio

ggplot() +
  geom_ribbon(data = pred_ratio_ci, aes(x = dbh, ymin = low, ymax = hi), alpha = 0.5) +
  geom_line(data = pred_ratio_ci, aes(x = dbh, y = med)) +
  geom_point(data = obs_ratio, aes(x = bin_midpoint, y = breeder_density_ratio), size = 4, pch = 21) +
  scale_x_log10(limits = c(1, 25)) + scale_y_log10(name = 'breeder:pioneer abundance ratio', limits = c(0.1, 100)) +
  theme_minimal()
  