# Script to clean up parameter CI CSV, removing duplicates and including slopes explicitly.

library(tidyverse)

params <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta')
density_params <- c('alpha', 'alpha_low', 'alpha_high', 'tau', 'alpha_mid', 'tau_low', 'tau_high')

# Remove duplicates.
production_param_df <- params %>% filter(parameter %in% production_params, dens_model == 3)
density_param_df <- params %>% filter(parameter %in% density_params, prod_model == 2)

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter')
density_params_names <- c('slope', 'slope small trees', 'slope large trees', 'cutoff small to large', 'slope middle trees', 'cutoff small to middle', 'cutoff middle to large')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# Convert parameters in density data frame to slopes.
density_param_df <- density_param_df %>%
  mutate_at(vars(mean:q975), ~ case_when(parameter %in% c('alpha', 'alpha_mid', 'alpha_high') ~ -(.x+1),
                                         parameter %in% c('alpha_low') ~ .x-1,
                                         TRUE ~ .x))

# correct to alternative names
production_param_df <- production_param_df %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(prod_model == 1 ~ 'one segment',
                           prod_model == 2 ~ 'two segment')) %>%
  select(-dens_model, -prod_model, -parameter) %>%
  select(year, fg, model, parameter_description, everything())

density_param_df <- density_param_df %>%
  mutate(parameter_description = density_params_names[match(parameter, density_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(dens_model == 1 ~ 'one segment',
                           dens_model == 2 ~ 'two segment',
                           dens_model == 3 ~ 'three segment')) %>%
  select(-dens_model, -prod_model, -parameter) %>%
  select(year, fg, model, parameter_description, everything())

write.csv(production_param_df, '~/google_drive/ForestLight/data/data_piecewisefits/clean_parameters_production.csv', row.names = FALSE)
write.csv(density_param_df, '~/google_drive/ForestLight/data/data_piecewisefits/clean_parameters_density.csv', row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  filter(dens_model == 3) %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                    prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2)

# Clean output of number of individuals and species by fg
n_ind <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/tally_indiv_by_fg.csv', stringsAsFactors = FALSE)
n_spp <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/tally_spp_by_fg.csv', stringsAsFactors = FALSE)

n_ind <- n_ind %>% 
  filter(year == 1995) %>%
  left_join(n_spp) %>%
  mutate(fg = fg_full_names[match(fg, fgs)])

write.csv(r2df, '~/google_drive/ForestLight/data/data_piecewisefits/clean_rsquared_production.csv', row.names = FALSE)
write.csv(n_ind, '~/google_drive/ForestLight/data/data_piecewisefits/clean_N_by_functionalgroup.csv', row.names = FALSE)
