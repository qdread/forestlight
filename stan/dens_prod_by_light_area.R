# Fit models for truncated power law at 10 W m-2 for FG1,2,3,4, individual production vs light/area and density vs light/area

library(tidyverse)
library(rstan)
library(forestscaling)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))

# Load data
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r')) # doesn't include imputed values

# What is the truncation point
#ggplot(alltree_light_95, aes(x = light_received_byarea)) + geom_histogram() + scale_x_log10() + scale_y_log10() + theme_bw()
# It's clearly ~ 10
table(round(alltree_light_95$light_received_byarea))[1:20] # Use 7 as the number since it is the mode of the distribution.
lightbin <- logbin(alltree_light_95$light_received_byarea, n = 20)
ggplot(lightbin, aes(x = bin_midpoint, y = bin_value)) + geom_point() + scale_x_log10() + scale_y_log10() # Still looks like 7.

# Create data objects
x_min <- 7

get_stan_data <- function(dat, x_min) with(dat, list(N = nrow(dat), x = dat$light_received_byarea, y = dat$production, x_min = x_min))

stan_data_list <- alltree_light_95 %>%
  filter(!recruit) %>%
  filter(!fg %in% 5, !is.na(fg), light_received_byarea >= x_min) %>%
  mutate(fg = paste0('fg', fg)) %>%
  group_by(fg) %>%
  group_map(~ get_stan_data(., x_min))

# Compile models
# Density 1 model will estimate x_min from data too.
mod_dens1 <- stan_model('stan/clean_workflow/model_scripts/density1_simplified.stan')
mod_prod1 <- stan_model('stan/clean_workflow/model_scripts/production1_nologlik.stan')

# Fit models

# Fitting options
n_chains <- 3
n_iter <- 6000
n_warmup <- 5000

prodfit_alltrees <- map(stan_data_list, ~ sampling(mod_prod1, data = ., chains = n_chains, iter = n_iter, warmup = n_warmup))
densfit_alltrees <- map(stan_data_list, ~ sampling(mod_dens1, data = ., chains = n_chains, iter = n_iter, warmup = n_warmup))

# Pull out the summaries to make sure models look OK
prodfit_summaries <- map(prodfit_alltrees, ~ summary(.)$summary)
densfit_summaries <- map(densfit_alltrees, ~ summary(.)$summary)

# Save fits
save(prodfit_alltrees, densfit_alltrees, prodfit_summaries, densfit_summaries, file = '~/Dropbox/Q/projects/forestlight/fits_bylight_forratio.RData')


#### Extract model output to get the fitted values, slopes, etc.
load('~/Dropbox/Q/projects/forestlight/fits_bylight_forratio.RData')

# source the extra extraction functions that aren't in the package
source('~/Documents/GitHub/forestscalingworkflow/R_functions/model_output_extraction_functions.r')
source('~/Documents/GitHub/forestlight/stan/get_ratio_slopes_fromfit.R')


# Get the statistics on the ratio trends.
la_pred <- logseq(1,412,101)
parnames_prod <- c('beta0', 'beta1')

# Predicted values for each fit.
dens_pred_fg_lightarea <- map(densfit_alltrees, function(fit) {
  pars_fg <- extract(fit, c('alpha')) %>% bind_cols
  pmap(pars_fg, pdf_pareto, x = la_pred, xmin = x_min)
})
prod_pred_fg_lightarea <- map(prodfit_alltrees, function(fit) {
  pars_fg <- extract(fit, parnames_prod) %>% bind_cols
  pmap(pars_fg, powerlaw_log, x = la_pred)
})

# Get total production including correction factors.

corr_factor <- function(y, y_fit, n_pars) {
  y_fit <- do.call(cbind, y_fit)
  # Get residuals by subtracting log y from linear predictor
  resids <- -1 * sweep(log(y_fit), 1, log(y))
  
  # Sum of squared residuals
  ssq_resid <- apply(resids^2, 2, sum)
  # Standard error of estimates
  sse <- (ssq_resid / (length(y) - n_pars))^0.5
  # Correction factors
  exp((sse^2)/2)
}

# We need all fitted values for production, not just the 101 values, to calculate correction factor.
prod_pred_all_lightarea <- map2(stan_data_list, prodfit_alltrees, function(dat, fit) {
  pars_fg <- extract(fit, parnames_prod) %>% bind_cols
  pmap(pars_fg, powerlaw_log, x = dat$x)
})

prod_cf_fg_lightarea <- map2(stan_data_list, prod_pred_all_lightarea, ~ corr_factor(y = .x$y, y_fit = .y, n_pars = 2))

totalprod_pred_fg_lightarea <- map2(dens_pred_fg_lightarea, prod_pred_fg_lightarea, ~ do.call(cbind, .x) * do.call(cbind, .y)) %>%
  map2(prod_cf_fg_lightarea, ~ sweep(.x, 2, .y, `*`))

# Take the ratio of the predicted values from each sampling iteration
# Multiply by number of individuals
dens_ratio13_lightarea <- map2(dens_pred_fg_lightarea[[1]], dens_pred_fg_lightarea[[3]], ~ (.x * stan_data_list[[1]]$N) / (.y * stan_data_list[[3]]$N) )
dens_ratio24_lightarea <- map2(dens_pred_fg_lightarea[[2]], dens_pred_fg_lightarea[[4]], ~ (.x * stan_data_list[[2]]$N) / (.y * stan_data_list[[4]]$N))

totalprod_ratio13_lightarea <- totalprod_pred_fg_lightarea[[1]] / totalprod_pred_fg_lightarea[[3]] * (stan_data_list[[1]]$N / stan_data_list[[3]]$N)
totalprod_ratio24_lightarea <- totalprod_pred_fg_lightarea[[2]] / totalprod_pred_fg_lightarea[[4]] * (stan_data_list[[2]]$N / stan_data_list[[4]]$N)

# Generate credible intervals from the ratios
dens_ratio_ci13_lightarea <- do.call(cbind, dens_ratio13_lightarea) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(light_area = la_pred) 

dens_ratio_ci24_lightarea <- do.call(cbind, dens_ratio24_lightarea) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(light_area = la_pred) 

totalprod_ratio_ci13_lightarea <- totalprod_ratio13_lightarea %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(light_area = la_pred) 

totalprod_ratio_ci24_lightarea <- totalprod_ratio24_lightarea %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  t %>% as.data.frame %>% setNames(c('low','med','hi')) %>%
  mutate(light_area = la_pred) 

# Save fitted values and credible intervals
ratio_fitted_lightarea <- map2_dfr(list(data.frame(ratio = 'fast:slow', variable = 'density'),
                                    data.frame(ratio = 'pioneer:breeder', variable = 'density'),
                                    data.frame(ratio = 'fast:slow', variable = 'total production'),
                                    data.frame(ratio = 'pioneer:breeder', variable = 'total production')),
                                   list(dens_ratio_ci13_lightarea, dens_ratio_ci24_lightarea, totalprod_ratio_ci13_lightarea, totalprod_ratio_ci24_lightarea),
                                   ~data.frame(.x, .y)) %>%
  setNames(c('ratio','variable','q025','q50','q975','light_area'))

write_csv(ratio_fitted_lightarea, file.path(gdrive_path, 'data/data_piecewisefits/ratio_fittedvalues_lightarea.csv'))


# Plotting ----------------------------------------------------------------


obs_ratio13_lightarea <- read_csv(file.path(gdrive_path, 'data/data_binned/fastslow_stats_bylight_byyear.csv')) %>%
  filter(n_individuals > 20, year == 1995)
obs_ratio24_lightarea <- read_csv(file.path(gdrive_path, 'data/data_binned/breeder_stats_bylight_byyear.csv')) %>%
  filter(n_individuals > 20, year == 1995)

ggplot() +
  geom_ribbon(data = dens_ratio_ci24_lightarea, aes(x = light_area, ymin = low, ymax = hi), alpha = 0.5) +
  geom_ribbon(data = dens_ratio_ci13_lightarea, aes(x = light_area, ymin = low, ymax = hi), alpha = 0.5) +
  geom_line(data = dens_ratio_ci24_lightarea, aes(x = light_area, y = med)) +
  geom_line(data = dens_ratio_ci13_lightarea, aes(x = light_area, y = med)) +
  geom_point(data = obs_ratio24_lightarea, aes(x = bin_midpoint, y = breeder_density_ratio), size = 4, pch = 21) +
  geom_point(data = obs_ratio13_lightarea, aes(x = bin_midpoint, y = fastslow_density_ratio), size = 4, pch = 21) +
  scale_x_log10(limits = c(1, 413)) + scale_y_log10(name = 'abundance ratio', limits = c(0.01, 200)) +
  theme_minimal()


ggplot() +
  geom_ribbon(data = totalprod_ratio_ci24_lightarea %>% filter(light_area>7), aes(x = light_area, ymin = low, ymax = hi), alpha = 0.5) +
  geom_ribbon(data = totalprod_ratio_ci13_lightarea %>% filter(light_area>7), aes(x = light_area, ymin = low, ymax = hi), alpha = 0.5) +
  geom_line(data = totalprod_ratio_ci24_lightarea %>% filter(light_area>7), aes(x = light_area, y = med)) +
  geom_line(data = totalprod_ratio_ci13_lightarea %>% filter(light_area>7), aes(x = light_area, y = med)) +
  geom_point(data = obs_ratio24_lightarea, aes(x = bin_midpoint, y = breeder_production_ratio), size = 4, pch = 21) +
  geom_point(data = obs_ratio13_lightarea, aes(x = bin_midpoint, y = fastslow_production_ratio), size = 4, pch = 21) +
  scale_x_log10(limits = c(3, 413), breaks = c(3,30,300)) + scale_y_log10(name = 'production ratio', limits = c(0.01, 300)) +
  theme_minimal()
