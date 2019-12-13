# Script to clean up all new piecewise output and put in tidy CSVs
# QDR / Forestlight / 14 June 2019

# Modified 20 June: add light params and R2s.
# Modified 10 Sep: include sigma in production fits.

fp_out <- '~/google_drive/ForestLight/data/clean_summary_tables'

library(tidyverse)

params <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

density_param_df <- params %>% filter(variable == 'density')
production_param_df <- params %>% filter(variable == 'production')

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta','sigma')
density_params <- c('alpha', 'alpha_low', 'alpha_high', 'tau', 'alpha_mid', 'tau_low', 'tau_high')

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter', 'sigma')
density_params_names <- c('slope', 'slope small trees', 'slope large trees', 'cutoff small to large', 'slope middle trees', 'cutoff small to middle', 'cutoff middle to large')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# Convert parameters in density data frame to slopes.
density_param_df <- density_param_df %>%
  mutate_at(vars(matches('^mean|^q')), ~ case_when(parameter %in% c('alpha', 'alpha_mid', 'alpha_high') ~ -(.x+1),
                                         parameter %in% c('alpha_low') ~ .x-1,
                                         TRUE ~ .x))

# correct to alternative names
production_param_df <- production_param_df %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')

density_param_df <- density_param_df %>%
  mutate(parameter_description = density_params_names[match(parameter, density_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment',
                           model == 3 ~ 'three segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')

write.csv(production_param_df, file.path(fp_out, 'clean_parameters_production.csv'), row.names = FALSE)
write.csv(density_param_df, file.path(fp_out, 'clean_parameters_density.csv'), row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                    prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_production.csv'), row.names = FALSE)

# Clean output of informatiion criteria
ics <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/piecewise_ics_by_fg.csv', stringsAsFactors = FALSE)

ics <- ics %>%
  filter(criterion == 'WAIC') %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment',
                      dens_model == 1 ~ 'one segment',
                      dens_model == 2 ~ 'two segment',
                      dens_model == 3 ~ 'three segment'),
    fg = fg_full_names[match(fg, fgs)]
  ) %>%
  rename(WAIC = IC_value,
         stderr_WAIC = IC_stderr) %>%
  select(year, variable, fg, model, WAIC, stderr_WAIC) %>%
  filter(!fg %in% 'unclassified')

write.csv(ics, file.path(fp_out, 'clean_WAIC_by_functionalgroup.csv'), row.names = FALSE)


# Individual light output -------------------------------------------------

library(tidyverse)

params <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/light_piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta','sigma')

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter','sigma')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# correct to alternative names
indivlight_param_df <- params %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')


write.csv(indivlight_param_df, file.path(fp_out, 'clean_parameters_individuallight.csv'), row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/light_piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_individuallight.csv'), row.names = FALSE)

# Clean output of information criteria
ics <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/light_piecewise_ics_by_fg.csv', stringsAsFactors = FALSE)

ics <- ics %>%
  filter(criterion == 'WAIC') %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'
                      ),
    fg = fg_full_names[match(fg, fgs)]
  ) %>%
  rename(WAIC = IC_value,
         stderr_WAIC = IC_stderr) %>%
  select(year, variable, fg, model, WAIC, stderr_WAIC) %>%
  filter(!fg %in% 'unclassified')

write.csv(ics, file.path(fp_out, 'clean_WAIC_by_functionalgroup_individuallight.csv'), row.names = FALSE)

# Volume output -------------------------------------------------

library(tidyverse)

params <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/volume_piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta','sigma')

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter','sigma')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# correct to alternative names
indivlight_param_df <- params %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')


write.csv(indivlight_param_df, file.path(fp_out, 'clean_parameters_crownvolume.csv'), row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/volume_piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_crownvolume.csv'), row.names = FALSE)

# Clean output of information criteria
ics <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/volume_piecewise_ics_by_fg.csv', stringsAsFactors = FALSE)

ics <- ics %>%
  filter(criterion == 'WAIC') %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'
    ),
    fg = fg_full_names[match(fg, fgs)]
  ) %>%
  rename(WAIC = IC_value,
         stderr_WAIC = IC_stderr) %>%
  select(year, variable, fg, model, WAIC, stderr_WAIC) %>%
  filter(!fg %in% 'unclassified')

write.csv(ics, file.path(fp_out, 'clean_WAIC_by_functionalgroup_crownvolume.csv'), row.names = FALSE)


# Diameter growth output --------------------------------------------------

library(tidyverse)

params <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/diamgrowth_piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)

production_params <- c('beta0','beta1','beta1_low','beta1_high','x0','delta','sigma')

# Alternative names.
production_params_names <- c('intercept', 'slope', 'slope small trees', 'slope large trees', 'cutoff', 'smoothing parameter', 'sigma')

# Full names of functional groups
fg_full_names <- c('fast', 'large pioneer', 'slow', 'small breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

# correct to alternative names
indivdiamgrowth_param_df <- params %>%
  mutate(parameter_description = production_params_names[match(parameter, production_params)],
         fg = fg_full_names[match(fg, fgs)],
         model = case_when(model == 1 ~ 'one segment',
                           model == 2 ~ 'two segment')) %>%
  select(-parameter) %>%
  select(year, fg, model, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')


write.csv(indivdiamgrowth_param_df, file.path(fp_out, 'clean_parameters_individualdiametergrowth.csv'), row.names = FALSE)

# Clean output of r-squared for production fits
r2df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/diamgrowth_piecewise_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'),
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(year, fg, model, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_individualdiametergrowth.csv'), row.names = FALSE)

# Clean output of information criteria
ics <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/diamgrowth_piecewise_ics_by_fg.csv', stringsAsFactors = FALSE)

ics <- ics %>%
  filter(criterion == 'WAIC') %>%
  mutate(
    model = case_when(prod_model == 1 ~ 'one segment',
                      prod_model == 2 ~ 'two segment'
    ),
    fg = fg_full_names[match(fg, fgs)]
  ) %>%
  rename(WAIC = IC_value,
         stderr_WAIC = IC_stderr) %>%
  select(year, variable, fg, model, WAIC, stderr_WAIC) %>%
  filter(!fg %in% 'unclassified')

write.csv(ics, file.path(fp_out, 'clean_WAIC_by_functionalgroup_individualdiametergrowth.csv'), row.names = FALSE)


# Mortality regression ----------------------------------------------------

mort_pars <- read_csv('~/google_drive/ForestLight/data/data_piecewisefits/mortality_paramci_by_fg.csv')

mort_pars <- filter(mort_pars, !parameter %in% 'lp__', !grepl('fg', parameter))

mort_pars$parameter <- c('global mean intercept', 'intercept std. dev', 'global mean slope', 'slope std.dev', rep('intercept', 5), rep('slope', 5))

mort_pars$fg <- c(rep('--', 4), rep(fg_full_names[1:5], 2))

mort_pars <- mort_pars %>% select(fg, parameter, everything())

names(mort_pars)[6:10] <- c('q025', 'q25', 'q50', 'q75', 'q975')

write_csv(mort_pars, file.path(fp_out, 'clean_parameters_mortality.csv'))


# Production per area by light per area -----------------------------------

# parameters
la_pars <- read_csv('~/google_drive/ForestLight/data/data_piecewisefits/lightbyarea_paramci_by_fg.csv')

la_par_description <- c('A','b','k','light per area at maximum slope of curve', 'growth per area at maximum slope of curve', 'maximum slope of curve in log-log space')

la_pars_df <- la_pars %>%
  mutate(parameter_description = rep(la_par_description,7),
         fg = fg_full_names[match(fg, fgs)]) %>%
  select(-parameter) %>%
  select(year, fg, parameter_description, everything()) %>%
  filter(!fg %in% 'unclassified')

write_csv(la_pars_df, file.path(fp_out, 'clean_parameters_growthperareabylightperarea.csv'))

# r squared values
r2df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/lightbyarea_r2_by_fg.csv', stringsAsFactors = FALSE)

r2df <- r2df %>%
  mutate(
    fg = fg_full_names[match(fg, fgs)],
    r2 = paste0(round(q50, 3), ' [', round(q025, 3), ',', round(q975, 3), ']')) %>%
  select(fg, r2) %>%
  filter(!fg %in% 'unclassified')

write.csv(r2df, file.path(fp_out, 'clean_rsquared_growthperareabylightperarea.csv'), row.names = FALSE)
