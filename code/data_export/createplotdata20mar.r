# Combine modeled and observed data into separate data frames for plotting
# Edit 17 Apr: Do midsize trees as well.

# Observed data -----------------------------------------------------------



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


# Modeled data -----------------------------------------------------------

library(dplyr)

ci_df <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

pred_dens <- ci_df %>%
  filter(variable == 'density') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

pred_indivprod <- ci_df %>%
  filter(variable == 'production') %>%
  select(-variable)

pred_totalprod <- ci_df %>%
  filter(variable == 'total_production') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

fp <- 'C:/Users/Q/google_drive/ForestLight/data/data_forplotting_12apr2018'

write.csv(pred_dens, file.path(fp, 'pred_dens.csv'), row.names = FALSE)
write.csv(pred_indivprod, file.path(fp, 'pred_indivprod.csv'), row.names = FALSE)
write.csv(pred_totalprod, file.path(fp, 'pred_totalprod.csv'), row.names = FALSE)

# Observed data midsize ---------------------------------------------------

# Put bins by year and by functional group into single data frames

years <- seq(1990, 2010, 5)
fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")
size_limits <- c(5, 50)
numbins <- 20

library(purrr)
library(dplyr)
source('code/allfunctions27july.r')
source('code/fakebin.r')

# Bin dbh and total production with the size limits.
allmidsize <- map(alltreedat[2:6], function(x) filter(x, between(dbh_corr, size_limits[1], size_limits[2])))
alldbh <- do.call(c, map(allmidsize, 'dbh_corr'))
alldbhbin <- logbin(x=alldbh, n=numbins)
dbhbin_all_byyear <- lapply(allmidsize, function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = alldbhbin))
totalprodbin_alltree_byyear <- lapply(allmidsize, function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = alldbhbin))

# Make a version of alltreedat without the unclassified trees
allmidsize_classified <- lapply(allmidsize, function(x) subset(x, !is.na(fg)))
# Bin classified trees. (log binning of density)
alldbh_classified <- unlist(lapply(allmidsize_classified, '[', , 'dbh_corr'))
dbhbin_allclassified <- logbin(x = alldbh_classified, y = NULL, n = numbins)

fgmidsize <- map(fgdat, function(x) map(x, function(y) filter(y, between(dbh_corr, size_limits[1], size_limits[2]))))

dbhbin_fg_midsize <- list()

for (i in 1:6) {
  dbhbin_fg_midsize[[i]] <- lapply(fgmidsize[[i]][2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_allclassified))
}

totalprodbin_fg_midsize <- list()

for (i in 1:6) {
  totalprodbin_fg_midsize[[i]] <- lapply(fgmidsize[[i]][2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_allclassified))
}


# Observed density
obs_dens_df <- cbind(fg = 'all', year = rep(years, each = numbins), do.call(rbind, dbhbin_all_byyear))
obs_dens_fg <- map(dbhbin_fg_midsize, function(x) map2(x, years, function(z, yz) cbind(year = yz, z)))
obs_dens_fg <- map2(obs_dens_fg, fg_names, function(x, nx) map(x, function(x) cbind(fg = nx, x)))
obs_dens_fg <- do.call(rbind, map(obs_dens_fg, function(x) do.call(rbind, x)))
obs_dens_df <- rbind(obs_dens_df, obs_dens_fg)

# Observed total production
obs_totalprod_df <- cbind(fg = 'all', year = rep(years, each = numbins), do.call(rbind, totalprodbin_alltree_byyear))
obs_totalprod_fg <- map(totalprodbin_fg_midsize, function(x) map2(x, years, function(z, yz) cbind(year = yz, z)))
obs_totalprod_fg <- map2(obs_totalprod_fg, fg_names, function(x, nx) map(x, function(x) cbind(fg = nx, x)))
obs_totalprod_fg <- do.call(rbind, map(obs_totalprod_fg, function(x) do.call(rbind, x)))
obs_totalprod_df <- rbind(obs_totalprod_df, obs_totalprod_fg)

# Observed individual production
# "Fake bin" for each functional group and each year separately.

obs_indivprod_df <- map2(allmidsize, dbhbin_all_byyear, 
                         function(x, y) fakebin_across_years(dat_values = x$production, dat_classes = x$dbh_corr, edges = y, n_census = 1))
obs_indivprod_df <- cbind(fg = 'all', year = rep(years, each = numbins), do.call(rbind, obs_indivprod_df))

obs_indivprod_fg <- replicate(length(fg_names), list())

for (i in 1:length(fg_names)) {
  for (j in 1:length(years)) {
    obs_indivprod_fg[[i]][[j]] <-
      cbind(fg = fg_names[i], year = years[j],
            fakebin_across_years(dat_values = fgmidsize[[i]][[j+1]]$production,
                                 dat_classes = fgmidsize[[i]][[j+1]]$dbh_corr,
                                 edges = dbhbin_fg_midsize[[i]][[j]],
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
fp <- 'C:/Users/Q/google_drive/ForestLight/data/data_forplotting_12apr2018'

write.csv(obs_dens_df, file.path(fp, 'obs_dens_midsize.csv'), row.names = FALSE)
write.csv(obs_indivprod_df, file.path(fp, 'obs_indivprod_midsize.csv'), row.names = FALSE)
write.csv(obs_totalprod_df, file.path(fp, 'obs_totalprod_midsize.csv'), row.names = FALSE)

# Modeled data midsize ---------------------------------------------------

library(dplyr)

ci_df <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/ci_by_fg_midsizetrees.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

pred_dens <- ci_df %>%
  filter(variable == 'density') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

pred_indivprod <- ci_df %>%
  filter(variable == 'production') %>%
  select(-variable)

pred_totalprod <- ci_df %>%
  filter(variable == 'total_production') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

fp <- 'C:/Users/Q/google_drive/ForestLight/data/data_forplotting_12apr2018'

write.csv(pred_dens, file.path(fp, 'pred_dens_midsize.csv'), row.names = FALSE)
write.csv(pred_indivprod, file.path(fp, 'pred_indivprod_midsize.csv'), row.names = FALSE)
write.csv(pred_totalprod, file.path(fp, 'pred_totalprod_midsize.csv'), row.names = FALSE)
