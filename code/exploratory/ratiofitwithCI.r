# Test ratio fit

library(tidyverse)
library(forestscaling)
library(rstan)

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

options(mc.cores = 2)
rstan_options(auto_write = TRUE)

# Stan model
d1mod <- stan_model('~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/density1.stan')
d2mod <- stan_model('~/Documents/GitHub/forestlight/stan/clean_workflow/model_scripts/density2.stan')

# Data for two functional groups, fg2 and fg4

dat_fg2 <- with(alltreedat[[3]] %>% filter(fg == 2), list(x = dbh_corr, y = production, N = length(dbh_corr), x_min = min(dbh_corr), x_max = max(dbh_corr)))

dat_fg4 <- with(alltreedat[[3]] %>% filter(fg == 4), list(x = dbh_corr, y = production, N = length(dbh_corr), x_min = min(dbh_corr), x_max = max(dbh_corr)))

fit_fg2 <- sampling(d2mod, data = dat_fg2, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222)
fit_fg4 <- sampling(d2mod, data = dat_fg4, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222)

fit1_fg2 <- sampling(d1mod, data = dat_fg2, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222)
fit1_fg4 <- sampling(d1mod, data = dat_fg4, chains = 2, iter = 2000, warmup = 1000, thin = 1, seed = 222)

# Generate predicted values for each fit
dbh_pred <- logseq(1,286,101)
parnames <- c('tau', 'alpha_low', 'alpha_high')

pars_fg2 <- extract(fit_fg2, parnames) %>% bind_cols
pred_fg2 <- pmap(pars_fg2, pdf_2part, x = dbh_pred, xmin = 1)

pars_fg4 <- extract(fit_fg4, parnames) %>% bind_cols
pred_fg4 <- pmap(pars_fg4, pdf_2part, x = dbh_pred, xmin = 1)
###
pars_fg2 <- extract(fit1_fg2, 'alpha') %>% bind_cols
pred_fg2 <- pmap(pars_fg2, pdf_pareto, x = dbh_pred, xmin = 1)

pars_fg4 <- extract(fit1_fg4, 'alpha') %>% bind_cols
pred_fg4 <- pmap(pars_fg4, pdf_pareto, x = dbh_pred, xmin = 1)

# Take the ratio of the predicted values from each sampling iteration
pred_ratio <- map2(pred_fg2, pred_fg4, `/`)

# Generate credible intervals from the ratios
pred_ratio_ci <- do.call(cbind, pred_ratio) %>%
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
  