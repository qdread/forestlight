# Code to fit binned values of richness, normalized by dividing by bin width, to diameter
# QDR / Forestlight / 27 July 2020


# Load data ---------------------------------------------------------------

# See lines 1-50 of richness_ratio_trend.R
# This loads the data and sets up the binned values with n=20

library(rstan)
options(mc.cores = 3)

# Set up binned values for "all" ------------------------------------------

bin_all <- data.frame(fg = 'all', bin = 1:20) %>%
  left_join(bin_bounds %>% mutate(bin = 1:20))

bin_all_light <- data.frame(fg = 'all', bin = 1:20) %>%
  left_join(bin_bounds_light %>% mutate(bin = 1:20))

bin_all <- bin_all %>%
  cbind(pmap_dfr(bin_all, function(bin_min, bin_max, ...) {
    sp_ids <- as.character(dat$sp[dat$dbh_corr >= bin_min & dat$dbh_corr < bin_max])
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min))

bin_all_light <- bin_all_light %>%
  cbind(pmap_dfr(bin_all_light, function(bin_min, bin_max, ...) {
    sp_ids <- as.character(dat_light$sp[dat_light$light_received_byarea >= bin_min & dat_light$light_received_byarea < bin_max])
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min))

# Load model --------------------------------------------------------------

# Model is log(richness by bin width) ~ log(diameter), or log(light per area), with 3 segments, and FG as random effect
# Also fit with all, no mixed model (no random effect)

dat_all <- with(bin_all %>% filter(richness > 0), list(x = bin_midpoint, y = richness_by_bin_width, N = nrow(bin_all)))
dat_all_light <- with(bin_all_light %>% filter(richness > 0), list(x = bin_midpoint, y = richness_by_bin_width, N = nrow(bin_all_light)))

bin_x_fg_use <- bin_x_fg %>% filter(richness > 0, !fg %in% 'unclassified')
dat_fg <- with(bin_x_fg_use, list(x = bin_midpoint, y = richness_by_bin_width, fg = as.numeric(factor(fg)), N = nrow(bin_x_fg_use), M = 5))
bin_x_fg_light_use <- bin_x_fg %>% filter(richness > 0, !fg %in% 'unclassified')
dat_fg_light <- with(bin_x_fg_light_use, list(x = bin_midpoint, y = richness_by_bin_width, fg = as.numeric(factor(fg)), N = nrow(bin_x_fg_light_use), M = 5))

# Fit models --------------------------------------------------------------

mod_2seg_linear <- stan_model(file.path(github_path, 'forestlight/stan/richness_2segment_linearmodel.stan'))
mod_3seg_linear <- stan_model(file.path(github_path, 'forestlight/stan/richness_3segment_linearmodel.stan'))

fit_2seg_linear_all <- sampling(mod_2seg_linear, data = dat_all, chains = 3, iter = 5000, warmup = 4000, seed = 222)
fit_3seg_linear_all <- sampling(mod_3seg_linear, data = dat_all, chains = 3, iter = 5000, warmup = 4000, seed = 111)

fit_2seg_linear_all_light <- sampling(mod_2seg_linear, data = dat_all_light, chains = 3, iter = 5000, warmup = 4000, seed = 222)
fit_3seg_linear_all_light <- sampling(mod_3seg_linear, data = dat_all_light, chains = 3, iter = 5000, warmup = 4000, seed = 111)


# Mixed models: load model and data ---------------------------------------

mod_2seg_mixed <- stan_model(file.path(github_path, 'forestlight/stan/richness_2segment_mixedmodel.stan'))
mod_3seg_mixed <- stan_model(file.path(github_path, 'forestlight/stan/richness_3segment_mixedmodel.stan'))

fit_2seg_mixed_all <- sampling(mod_2seg_mixed, data = dat_fg, chains = 3, iter = 5000, warmup = 4000, seed = 333)
# More iterations needed to converge 3 segment mixed model.
fit_3seg_mixed_all <- sampling(mod_3seg_mixed, data = dat_fg, chains = 3, iter = 5000, warmup = 4000, seed = 444)

fit_2seg_mixed_all_light <- sampling(mod_2seg_mixed, data = dat_fg_light, chains = 3, iter = 5000, warmup = 4000, seed = 555)
fit_3seg_mixed_all_light <- sampling(mod_3seg_mixed, data = dat_fg_light, chains = 3, iter = 5000, warmup = 4000, seed = 666)
