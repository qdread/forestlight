# Create piecewise data for plotting.
# Modified 02 Jul 2019: Calculate new correction factors to get the proper normalizations (based on integral with bounded limits)

library(dplyr)
library(pracma)
gdrive_path <- '~/google_drive/ForestLight'
github_path <- '~/Documents/GitHub/forestlight'
fp_obs <- file.path(gdrive_path, 'data/data_forplotting_aug2018')

# Load raw data
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.R'))

source(file.path(github_path, 'code/allfunctions27july.r'))
# Load binned data
obs_totalprod <- read.csv(file.path(fp_obs, 'obs_totalprod.csv'), stringsAsFactors = FALSE)

binedgedata <- read.csv(file.path(fp_obs,'obs_dens.csv'),stringsAsFactors = FALSE) %>% filter(fg == 'all', year == 1995) 
area_core <- 42.84

totallightbins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = light_received, edges = binedgedata))
totalvolbins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = crownvolume, edges = binedgedata))

totallightbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$light_received, edges = binedgedata)) %>%
  ungroup %>%
  rbind(data.frame(fg = 'all', totallightbins_all)) %>%
  mutate(bin_value = bin_value / area_core, year = 1995)
totalvolbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$crownvolume, edges = binedgedata)) %>%
  ungroup %>%
  rbind(data.frame(fg = 'all', totalvolbins_all)) %>%
  mutate(bin_value = bin_value / area_core, year = 1995)


# Calculate limits of integral for each FG
lim_totalprod <- obs_totalprod %>%
  filter(year == 1995, bin_count >= 10) %>%
  group_by(fg) %>%
  summarize(lim_min = min(bin_min), lim_max = max(bin_max))
lim_totallight <-  totallightbins_fg %>%
  filter(year == 1995, bin_count >= 10) %>%
  group_by(fg) %>%
  summarize(lim_min = min(bin_min), lim_max = max(bin_max))
lim_totalvol <-  totalvolbins_fg %>%
  filter(year == 1995, bin_count >= 10) %>%
  group_by(fg) %>%
  summarize(lim_min = min(bin_min), lim_max = max(bin_max))

# Calculate total production, light, and volume within the limits
prod_totals_fg <- alltreedat[[3]] %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
  left_join(lim_totalprod) %>%
  group_by(fg) %>%
  filter(dbh_corr <= lim_max) %>%
  summarize(total = sum(production))

prod_total_all <- alltreedat[[3]] %>%
  filter(dbh_corr <= lim_totalprod$lim_max[lim_totalprod$fg == 'all']) %>%
  summarize(total = sum(production))

prod_totals_fg <- rbind(prod_totals_fg, data.frame(fg = 'all', total = prod_total_all$total))

light_totals_fg <- alltree_light_95 %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
  left_join(lim_totallight) %>%
  group_by(fg) %>%
  filter(dbh_corr <= lim_max) %>%
  summarize(total = sum(light_received))

light_total_all <- alltree_light_95 %>%
  filter(dbh_corr <= lim_totallight$lim_max[lim_totallight$fg == 'all']) %>%
  summarize(total = sum(light_received))

light_totals_fg <- rbind(light_totals_fg, data.frame(fg = 'all', total = light_total_all$total))

volume_totals_fg <- alltreedat[[3]] %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
  left_join(lim_totalvol) %>%
  group_by(fg) %>%
  filter(dbh_corr <= lim_max) %>%
  summarize(total = sum(crownvolume))

volume_total_all <- alltreedat[[3]] %>%
  filter(dbh_corr <= lim_totalvol$lim_max[lim_totalvol$fg == 'all']) %>%
  summarize(total = sum(crownvolume))

volume_totals_fg <- rbind(volume_totals_fg, data.frame(fg = 'all', total = volume_total_all$total))


ci_df <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/newpiecewise_ci_by_fg.csv'), stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

pred_dens <- ci_df %>%
  filter(variable == 'density') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), ~ ./area_core)

fitted_indivprod <- ci_df %>%
  filter(variable == 'production_fitted') %>%
  select(-variable)

fitted_totalprod <- ci_df %>%
  filter(variable == 'total_production_fitted') %>%
  select(-variable) 

pred_indivprod <- ci_df %>%
  filter(variable == 'production') %>%
  select(-variable)

pred_totalprod <- ci_df %>%
  filter(variable == 'total_production') %>%
  select(-variable) 

# Get integrals for total production and total production fitted, then multiply the fitted values by the normalization constant
# Norm constant is the summed bin total production / the integral
# Then divide by area
totalprod_integrals <- fitted_totalprod %>%
  left_join(lim_totalprod) %>%
  group_by(fg, dens_model, prod_model) %>%
  filter(dbh <= lim_max) %>%
  summarize(integral = trapz(x = dbh, y = q50)) %>%
  left_join(prod_totals_fg) %>%
  mutate(constant = total/integral)

fitted_totalprod <- fitted_totalprod %>%
  left_join(totalprod_integrals) %>%
  mutate_at(vars(starts_with('q')), ~ . * (constant/area_core))

pred_totalprod <- pred_totalprod %>%
  left_join(totalprod_integrals) %>%
  mutate_at(vars(starts_with('q')), ~ . * (constant/area_core))

fp <- '~/google_drive/ForestLight/data/data_piecewisefits'

write.csv(pred_dens, file.path(fp, 'pred_dens.csv'), row.names = FALSE)
write.csv(pred_indivprod, file.path(fp, 'pred_indivprod.csv'), row.names = FALSE)
write.csv(pred_totalprod, file.path(fp, 'pred_totalprod.csv'), row.names = FALSE)
write.csv(fitted_indivprod, file.path(fp, 'fitted_indivprod.csv'), row.names = FALSE)
write.csv(fitted_totalprod, file.path(fp, 'fitted_totalprod.csv'), row.names = FALSE)



#################
# For individual + total light scaling and total volume scaling
# Added 20 June 2019

library(dplyr)

ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/totallightscaling/light_piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

fitted_indivlight <- ci_df %>%
  filter(variable == 'incoming_light_fitted') %>%
  select(-variable)

fitted_totallight <- ci_df %>%
  filter(variable == 'total_incoming_light_fitted') %>%
  select(-variable) 

pred_indivlight <- ci_df %>%
  filter(variable == 'incoming_light') %>%
  select(-variable)

pred_totallight <- ci_df %>%
  filter(variable == 'total_incoming_light') %>%
  select(-variable) 

# Get integral and multiply fitted and predicted total light by the new constant
# Then divide by area
totallight_integrals <- fitted_totallight %>%
  left_join(lim_totallight) %>%
  group_by(fg, dens_model, prod_model) %>%
  filter(dbh <= lim_max) %>%
  summarize(integral = trapz(x = dbh, y = q50)) %>%
  left_join(light_totals_fg) %>%
  mutate(constant = total/integral)

fitted_totallight <- fitted_totallight %>%
  left_join(totallight_integrals) %>%
  mutate_at(vars(starts_with('q')), ~ . * (constant/area_core))

pred_totallight <- pred_totallight %>%
  left_join(totallight_integrals) %>%
  mutate_at(vars(starts_with('q')), ~ . * (constant/area_core))


fp <- '~/google_drive/ForestLight/data/data_piecewisefits'

write.csv(pred_indivlight, file.path(fp, 'pred_indivlight.csv'), row.names = FALSE)
write.csv(pred_totallight, file.path(fp, 'pred_totallight.csv'), row.names = FALSE)
write.csv(fitted_indivlight, file.path(fp, 'fitted_indivlight.csv'), row.names = FALSE)
write.csv(fitted_totallight, file.path(fp, 'fitted_totallight.csv'), row.names = FALSE)

## Volume
ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/volume_piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

fitted_totalvol <- ci_df %>%
  select(-variable) 

# Get integral and multiply fitted total volume by the new constant
# Then divide by area
totalvol_integrals <- fitted_totalvol %>%
  left_join(lim_totalvol) %>%
  group_by(fg, dens_model) %>%
  filter(dbh <= lim_max) %>%
  summarize(integral = trapz(x = dbh, y = q50)) %>%
  left_join(volume_totals_fg) %>%
  mutate(constant = total/integral)

fitted_totalvol <- fitted_totalvol %>%
  left_join(totalvol_integrals) %>%
  mutate_at(vars(starts_with('q')), ~ . * (constant/area_core))

fp <- '~/google_drive/ForestLight/data/data_piecewisefits/totallightscaling'

write.csv(fitted_totalvol, file.path(fp, 'fitted_totalvol.csv'), row.names = FALSE)
