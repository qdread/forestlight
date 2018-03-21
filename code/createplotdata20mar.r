# Combine modeled and observed data into a single data frame for plotting

# Put bins by year and by functional group into single data frames

years <- seq(1990, 2010, 5)
fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")

library(purrr)

# Observed density
obs_dens_df <- cbind(fg = 'all', year = rep(years, each = numbins), do.call(rbind, dbhbin_all_byyear))
obs_dens_fg <- map(dbhbin_fg_byyear, function(x) map2(x, years, function(z, yz) cbind(year = yz, z)))
obs_dens_fg <- map2(obs_dens_fg, fg_names, function(x, nx) map(x, function(x) cbind(fg = nx, x)))
obs_dens_fg <- do.call(rbind, map(obs_dens_fg, function(x) do.call(rbind, x)))
obs_dens_df <- rbind(obs_dens_df, obs_dens_fg)

# Observed total production
obs_totalprod_df <- cbind(fg = 'all', year = rep(years, each = numbins), do.call(rbind, totalprodbin_alltree_byyear))
obs_totalprod_fg <- map(totalprodbin_fg_byyear, function(x) map2(x, years, function(z, yz) cbind(year = yz, z)))
obs_totalprod_fg <- map2(obs_totalprod_fg, fg_names, function(x, nx) map(x, function(x) cbind(fg = nx, x)))
obs_totalprod_fg <- do.call(rbind, map(obs_totalprod_fg, function(x) do.call(rbind, x)))
obs_totalprod_df <- rbind(obs_totalprod_df, obs_totalprod_fg)

# Observed individual production
# "Fake bin" for each functional group and each year separately.

all_prod_1995 <- fakebin_across_years(dat_values = alltreedat[[3]]$production, dat_classes = alltreedat[[3]]$dbh_corr, edges = all_dens_1995, n_census = 1)

obs_indivprod_df <- map2(alltreedat[-1], dbhbin_all_byyear, 
     function(x, y) fakebin_across_years(dat_values = x$production, dat_classes = x$dbh_corr, edges = y, n_census = 1))
obs_indivprod_df <- cbind(fg = 'all', year = rep(years, each = numbins), do.call(rbind, obs_indivprod_df))

obs_indivprod_fg <- replicate(length(fg_names), list())

for (i in 1:length(fg_names)) {
  for (j in 1:length(years)) {
    obs_indivprod_fg[[i]][[j]] <-
      cbind(fg = fg_names[i], year = years[j],
      fakebin_across_years(dat_values = fgdat[[i]][[j+1]]$production,
                           dat_classes = fgdat[[i]][[j+1]]$dbh_corr,
                           edges = dbhbin_fg_byyear[[i]][[j]],
                           n_census = 1))
  }
}

obs_indivprod_fg <- do.call(rbind, map(obs_indivprod_fg, function(x) do.call(rbind, x)))
obs_indivprod_df <- rbind(obs_indivprod_df, obs_indivprod_fg)

# Divide the appropriate columns by area
area_core <- 42.84
obs_dens_df <- mutate(obs_dens_df, bin_value = bin_value/area_core)
obs_totalprod_df <- mutate(obs_totalprod_df, bin_value = bin_value/area_core)

# Write the data
fp <- 'C:/Users/Q/google_drive/ForestLight/data/data_20mar2018'

write.csv(obs_dens_df, file.path(fp, 'obs_dens.csv'), row.names = FALSE)
write.csv(obs_indivprod_df, file.path(fp, 'obs_indivprod.csv'), row.names = FALSE)
write.csv(obs_totalprod_df, file.path(fp, 'obs_totalprod.csv'), row.names = FALSE)

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

pred_dens_df <- ci_df %>%
  filter(variable == 'density') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

pred_indivprod_df <- ci_df %>%
  filter(variable == 'production') %>%
  select(-variable)

pred_totalprod_df <- ci_df %>%
  filter(variable == 'total_production') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

write.csv(pred_dens_df, file.path(fp, 'pred_dens.csv'), row.names = FALSE)
write.csv(pred_indivprod_df, file.path(fp, 'pred_indivprod.csv'), row.names = FALSE)
write.csv(pred_totalprod_df, file.path(fp, 'pred_totalprod.csv'), row.names = FALSE)
