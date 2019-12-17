# Create piecewise data for plotting.
# Modified 13 Dec 2019: Create both observed and predicted data in this script. 
# Modified 13 Dec 2019: Use the new up to date correction factor based on the true correction of Jensen's Inequality.
# Modified 02 Jul 2019: Calculate new correction factors to get the proper normalizations (based on integral with bounded limits)

library(tidyverse)
library(forestscaling)

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users',Sys.info()['user'],'Google Drive/ForestLight'))

fp_out <- file.path(gdrive_path, 'data/data_forplotting')

years <- seq(1990, 2010, 5)
fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")

# Load all the by-year binned data.
load(file.path(gdrive_path, 'data/data_binned/bin_object_singleyear.RData'))

# Load raw data
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.R'))

binedgedata <- densitybin_byyear %>% filter(fg == 'all', year == 1995) 
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

indivlightbins_all <- alltree_light_95 %>%
  mutate(indivlight_bin = cut(dbh_corr, breaks = c(binedgedata$bin_min[1], binedgedata$bin_max), labels = binedgedata$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(indivlight_bin) %>%
  do(c(n = nrow(.), quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))
indivlightbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  mutate(indivlight_bin =cut(dbh_corr, breaks = c(binedgedata$bin_min[1], binedgedata$bin_max), labels = binedgedata$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, indivlight_bin) %>%
  do(c(n = nrow(.), quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))
indivlightbins_fg <- data.frame(fg = 'all', indivlightbins_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(indivlightbins_fg)) %>%
  mutate(indivlight_bin = as.numeric(as.character(indivlight_bin))) %>%
  rename(bin_midpoint = indivlight_bin)

write.csv(indivlightbins_fg, file.path(fp_out, 'obs_indivlight.csv'), row.names = FALSE)
write.csv(totallightbins_fg, file.path(fp_out, 'obs_totallight.csv'), row.names = FALSE)
write.csv(totalvolbins_fg, file.path(fp_out, 'obs_totalvol.csv'), row.names = FALSE)

# Observed individual production
obs_indivprod_df <- map(alltreedat[-1],
                         function(x) fakebin_across_years(dat_values = x$production, dat_classes = x$dbh_corr, edges = binedgedata, n_census = 1))
obs_indivprod_df <- cbind(fg = 'all', year = rep(years, each = nrow(binedgedata)), do.call(rbind, obs_indivprod_df))

obs_indivprod_fg <- replicate(length(fg_names), list())

for (i in 1:length(fg_names)) {
  for (j in 1:length(years)) {
    obs_indivprod_fg[[i]][[j]] <-
      cbind(fg = fg_names[i], year = years[j],
            fakebin_across_years(dat_values = fgdat[[i]][[j+1]]$production,
                                 dat_classes = fgdat[[i]][[j+1]]$dbh_corr,
                                 edges = binedgedata,
                                 n_census = 1))
  }
}

obs_indivprod_fg <- do.call(rbind, map(obs_indivprod_fg, function(x) do.call(rbind, x)))
obs_indivprod_df <- rbind(obs_indivprod_df, obs_indivprod_fg)

# Load correction factor for total production
cf_production <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/piecewise_cf_by_fg.csv')) %>% 
  select(prod_model, fg, q50) %>% 
  rename(corr_factor = q50) %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))


# Convert to individuals and production per hectare.
obs_dens <- densitybin_byyear %>%
  mutate(bin_value = bin_value / area_core)
obs_totalprod <- totalproductionbin_byyear %>%
  mutate(bin_value = bin_value / area_core)


ci_df <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/piecewise_ci_by_fg.csv'), stringsAsFactors = FALSE)
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

fitted_totalprod <- fitted_totalprod %>%
  left_join(cf_production) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))

pred_totalprod <- pred_totalprod %>%
  left_join(cf_production) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))

write.csv(obs_dens, file.path(fp_out, 'obs_dens.csv'), row.names = FALSE)
write.csv(obs_totalprod, file.path(fp_out, 'obs_totalprod.csv'), row.names = FALSE)
write.csv(obs_indivprod_df, file.path(fp_out, 'obs_indivprod.csv'), row.names = FALSE)
write.csv(pred_dens, file.path(fp_out, 'pred_dens.csv'), row.names = FALSE)
write.csv(pred_indivprod, file.path(fp_out, 'pred_indivprod.csv'), row.names = FALSE)
write.csv(pred_totalprod, file.path(fp_out, 'pred_totalprod.csv'), row.names = FALSE)
write.csv(fitted_indivprod, file.path(fp_out, 'fitted_indivprod.csv'), row.names = FALSE)
write.csv(fitted_totalprod, file.path(fp_out, 'fitted_totalprod.csv'), row.names = FALSE)



#################
# For individual + total light scaling and total volume scaling
# Added 20 June 2019

library(dplyr)

ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/light_piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

cf_totallight<- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/light_piecewise_cf_by_fg.csv')) %>% 
  select(prod_model, fg, q50) %>% 
  rename(corr_factor = q50) %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))

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

fitted_totallight <- fitted_totallight %>%
  left_join(cf_totallight) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))

pred_totallight <- pred_totallight %>%
  left_join(cf_totallight) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))


write.csv(pred_indivlight, file.path(fp_out, 'pred_indivlight.csv'), row.names = FALSE)
write.csv(pred_totallight, file.path(fp_out, 'pred_totallight.csv'), row.names = FALSE)
write.csv(fitted_indivlight, file.path(fp_out, 'fitted_indivlight.csv'), row.names = FALSE)
write.csv(fitted_totallight, file.path(fp_out, 'fitted_totallight.csv'), row.names = FALSE)

## Volume
ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/volume_piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

cf_volume <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/volume_piecewise_cf_by_fg.csv')) %>% 
  select(prod_model, fg, q50) %>% 
  rename(corr_factor = q50) %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))

fitted_totalvol <- ci_df %>%
  filter(variable == 'total_crown_volume_fitted') %>%
  select(-variable) 

fitted_totalvol <- fitted_totalvol %>%
  left_join(cf_volume) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))

write.csv(fitted_totalvol, file.path(fp_out, 'fitted_totalvol.csv'), row.names = FALSE)


# Data for growth per area vs light per area ------------------------------

# In previous workflow this code was in createlightplotdata_final.r

dat90 <- alltree_light_90 %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

dat95 <- alltree_light_95 %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

numbins <- 20
light_bins_9095 <- logbin(x = c(dat90$light_area, dat95$light_area), n = numbins)

prod_light_bin_90 <- fakebin_across_years(dat_values = dat90$production_area,
                                          dat_classes = dat90$light_area,
                                          edges = light_bins_9095,
                                          n_census = 1)
prod_light_bin_95 <- fakebin_across_years(dat_values = dat95$production_area,
                                          dat_classes = dat95$light_area,
                                          edges = light_bins_9095,
                                          n_census = 1)

prod_light_bin_fg_90 <- dat90 %>%
  group_by(fg) %>%
  do(bin = fakebin_across_years(dat_values = .$production_area,
                                dat_classes = .$light_area,
                                edges = light_bins_9095,
                                n_census = 1))

prod_light_bin_fg_90 <- cbind(fg = rep(c('fg1','fg2','fg3','fg4','fg5','unclassified'), each = numbins),
                              do.call(rbind, prod_light_bin_fg_90$bin)) %>%
  filter(complete.cases(.))

prod_light_bin_fg_95 <- dat95 %>%
  group_by(fg) %>%
  do(bin = fakebin_across_years(dat_values = .$production_area,
                                dat_classes = .$light_area,
                                edges = light_bins_9095,
                                n_census = 1))

prod_light_bin_fg_95 <- cbind(fg = rep(c('fg1','fg2','fg3','fg4','fg5','unclassified'), each = numbins),
                              do.call(rbind, prod_light_bin_fg_95$bin)) %>%
  filter(complete.cases(.))

prod_light_bin_all <- rbind(
  cbind(year = 1990, fg = 'alltree', prod_light_bin_90),
  cbind(year = 1990, prod_light_bin_fg_90),
  cbind(year = 1995, fg = 'alltree', prod_light_bin_95),
  cbind(year = 1995, prod_light_bin_fg_95)
)

write.csv(prod_light_bin_all, file = file.path(fp_out, 'obs_light_binned.csv'), row.names = FALSE)

write.csv(rbind(
  cbind(year = 1990, dat90),
  cbind(year = 1995, dat95)
), 
file = file.path(fp_out, 'obs_light_raw.csv'), row.names = FALSE)

# Create predicted light by area from existing file

system2("cp", args = paste('~/google_drive/ForestLight/data/data_piecewisefits/lightbyarea_predci_by_fg.csv', file.path(fp_out, 'pred_light.csv')))


# Create bin data for MS figure 5 -----------------------------------------

alltree_light_95 <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), as.character(NA)))

dbhbin1995 <- with(alltree_light_95, logbin(x = dbh_corr, n = 20))

lightperareafakebin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  do(quantile(.$light_received/.$crownarea, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightperareafakebin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  do(quantile(.$light_received/.$crownarea, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightperareafakebin_fg <- data.frame(fg = 'all', lightperareafakebin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightperareafakebin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

lightpervolfakebin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  do(quantile(.$light_received/.$crownvolume, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightpervolfakebin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  do(quantile(.$light_received/.$crownvolume, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightpervolfakebin_fg <- data.frame(fg = 'all', lightpervolfakebin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightpervolfakebin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

unscaledlightbydbhfakebin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  do(quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

unscaledlightbydbhfakebin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  do(quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

unscaledlightbydbhfakebin_fg <- data.frame(fg = 'all', unscaledlightbydbhfakebin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(unscaledlightbydbhfakebin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

write.csv(lightperareafakebin_fg, file.path(fp_out, 'lightperareafakebin_fg.csv'), row.names = FALSE)
write.csv(lightpervolfakebin_fg, file.path(fp_out, 'lightpervolfakebin_fg.csv'), row.names = FALSE)
write.csv(unscaledlightbydbhfakebin_fg, file.path(fp_out, 'unscaledlightbydbhfakebin_fg.csv'), row.names = FALSE)


# Diameter growth plots ---------------------------------------------------


# Observed individual diameter growth
obs_indivdiamgrowth_df <- map(alltreedat[-1],
                        function(x) fakebin_across_years(dat_values = x$diam_growth_rate, dat_classes = x$dbh_corr, edges = binedgedata, n_census = 1))
obs_indivdiamgrowth_df <- cbind(fg = 'all', year = rep(years, each = nrow(binedgedata)), do.call(rbind, obs_indivdiamgrowth_df))

obs_indivdiamgrowth_fg <- replicate(length(fg_names), list())

for (i in 1:length(fg_names)) {
  for (j in 1:length(years)) {
    obs_indivdiamgrowth_fg[[i]][[j]] <-
      cbind(fg = fg_names[i], year = years[j],
            fakebin_across_years(dat_values = fgdat[[i]][[j+1]]$diam_growth_rate,
                                 dat_classes = fgdat[[i]][[j+1]]$dbh_corr,
                                 edges = binedgedata,
                                 n_census = 1))
  }
}

obs_indivdiamgrowth_fg <- do.call(rbind, map(obs_indivdiamgrowth_fg, function(x) do.call(rbind, x)))
obs_indivdiamgrowth_df <- rbind(obs_indivdiamgrowth_df, obs_indivdiamgrowth_fg)

# Fitted
fittedvals <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/diamgrowth_piecewise_ci_by_fg.csv'), stringsAsFactors = FALSE) %>%
  filter(variable == 'diameter_growth_fitted') %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))

write.csv(obs_indivdiamgrowth_df, file = file.path(fp_out, 'obs_indivdiamgrowth.csv'), row.names = FALSE)
write.csv(fittedvals, file = file.path(fp_out, 'fitted_indivdiamgrowth.csv'), row.names = FALSE)
