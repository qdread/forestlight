# Create piecewise data for plotting.
# Modified 13 Jan 2020: Use imputed production totals (1 segment production) to create total production bins (but not individual)
# Modified 09 Jan 2020: Fix incorrect number of individuals in some of the observed data.
# Modified 13 Dec 2019: Create both observed and predicted data in this script. 
# Modified 13 Dec 2019: Use the new up to date correction factor based on the true correction of Jensen's Inequality.
# Modified 02 Jul 2019: Calculate new correction factors to get the proper normalizations (based on integral with bounded limits)
p_out <- 'data/data_forplotting'

years <- seq(1990, 2010, 5)
fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")

# Load all the by-year binned data.

load(file.path(github_path, 'forestscalingworkflow/data/data_binned/bin_object_singleyear.RData'))

# Load raw data
load(file.path(gdrive_path, 'data/rawdataobj_withimputedproduction.RData'))

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

# Observed individual production
obs_indivprod_df <- map(alltreedat[-1],
                        function(x) with(x %>% filter(!recruit), cloudbin_across_years(dat_values = production, dat_classes = dbh_corr, edges = binedgedata, n_census = 1)))
obs_indivprod_df <- cbind(fg = 'all', year = rep(years, each = nrow(binedgedata)), do.call(rbind, obs_indivprod_df))

obs_indivprod_fg <- replicate(length(fg_names), list())

for (i in 1:length(fg_names)) {
  for (j in 1:length(years)) {
    obs_indivprod_fg[[i]][[j]] <-
      cbind(fg = fg_names[i], year = years[j],
            cloudbin_across_years(dat_values = fgdat[[i]][[j+1]] %>% filter(!recruit) %>% pull(production),
                                  dat_classes = fgdat[[i]][[j+1]] %>% filter(!recruit) %>% pull(dbh_corr),
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

# Added 09 Jan 2020: replace incorrect bin counts in observed production data with the bin counts from observed density.
joined_counts <- obs_indivprod_df %>% 
  select(fg, year, bin_midpoint) %>%
  left_join(obs_dens %>% select(fg, year, bin_midpoint, bin_count))

obs_indivprod_df$mean_n_individuals <- joined_counts$bin_count


#################
# For individual + total light scaling and total volume scaling
# Added 20 June 2019

library(dplyr)

ci_df <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/light_piecewise_cf_by_fg.csv'), stringsAsFactors = FALSE)
area_core <- 42.84
cf_totallight <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/light_piecewise_cf_by_fg.csv')) %>% 
  #cf_totallight<- read_csv('finalcsvs/light_piecewise_cf_by_fg.csv') %>% 
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

## Volume
ci_df <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/volume_piecewise_cf_by_fg.csv'))  
  
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


# Data for growth per area vs light per area ------------------------------

# In previous workflow this code was in createlightplotdata_final.r

dat90 <- alltree_light_90 %>%
  filter(!recruit) %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

dat95 <- alltree_light_95 %>%
  filter(!recruit) %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

numbins <- 20

light_bins_9095 <- logbin(x = c(dat90$light_area, dat95$light_area), n = numbins)

prod_light_bin_90 <- cloudbin_across_years(dat_values = dat90$production_area,
                                           dat_classes = dat90$light_area,
                                           edges = light_bins_9095,
                                           n_census = 1)
prod_light_bin_95 <- cloudbin_across_years(dat_values = dat95$production_area,
                                           dat_classes = dat95$light_area,
                                           edges = light_bins_9095,
                                           n_census = 1)

prod_light_bin_fg_90 <- dat90 %>%
  group_by(fg) %>%
  do(bin = cloudbin_across_years(dat_values = .$production_area,
                                 dat_classes = .$light_area,
                                 edges = light_bins_9095,
                                 n_census = 1))

prod_light_bin_fg_90 <- cbind(fg = rep(c('fg1','fg2','fg3','fg4','fg5','unclassified'), each = numbins),
                              do.call(rbind, prod_light_bin_fg_90$bin)) %>%
  filter(complete.cases(.))

prod_light_bin_fg_95 <- dat95 %>%
  group_by(fg) %>%
  do(bin = cloudbin_across_years(dat_values = .$production_area,
                                 dat_classes = .$light_area,
                                 edges = light_bins_9095,
                                 n_census = 1))

prod_light_bin_fg_95 <- cbind(fg = rep(c('fg1','fg2','fg3','fg4','fg5','unclassified'), each = numbins),
                              do.call(rbind, prod_light_bin_fg_95$bin)) %>%
  filter(complete.cases(.))

# Added 09 Jan 2020: Replace incorrect number of individuals in binned light data frame with correct ones by group.
# There is no existing bin count so must be created manually.
nlight90 <- dat90 %>%
  group_by(fg) %>%
  mutate(bin = cut(light_area, breaks = c(light_bins_9095$bin_min, light_bins_9095$bin_max[numbins]))) %>%
  group_by(fg, bin) %>%
  summarize(correct_n = n()) %>%
  mutate(bin_midpoint = light_bins_9095$bin_midpoint[as.numeric(bin)])
nlight95 <- dat95 %>%
  group_by(fg) %>%
  mutate(bin = cut(light_area, breaks = c(light_bins_9095$bin_min, light_bins_9095$bin_max[numbins]))) %>%
  group_by(fg, bin) %>%
  summarize(correct_n = n()) %>%
  mutate(bin_midpoint = light_bins_9095$bin_midpoint[as.numeric(bin)])

# Join 1990 and 1995 correct counts with the binned values.
joined_counts_light90 <- prod_light_bin_fg_90 %>%
  select(fg, bin_midpoint) %>%
  left_join(nlight90)
joined_counts_light95 <- prod_light_bin_fg_95 %>%
  select(fg, bin_midpoint) %>%
  left_join(nlight95)

prod_light_bin_fg_90$mean_n_individuals <- joined_counts_light90$correct_n
prod_light_bin_fg_95$mean_n_individuals <- joined_counts_light95$correct_n

prod_light_bin_all <- rbind(
  cbind(year = 1990, fg = 'alltree', prod_light_bin_90),
  cbind(year = 1990, prod_light_bin_fg_90),
  cbind(year = 1995, fg = 'alltree', prod_light_bin_95),
  cbind(year = 1995, prod_light_bin_fg_95)
)



# Create bin data for MS figure 5 -----------------------------------------

# Edited 01 Apr 2020 to add (geometric) means, fix NA bin.

light_bin_stats <- function(indivs) {
  qs <- quantile(indivs, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
    setNames(c('q025','q25','q50','q75','q975'))
  ci_width <- qnorm(0.975) * sd(log(indivs))/sqrt(length(indivs))
  data.frame(t(qs), 
             mean = exp(mean(log(indivs))),
             sd = exp(sd(log(indivs))),
             ci_min = exp(mean(log(indivs)) - ci_width),
             ci_max = exp(mean(log(indivs)) + ci_width))
}

alltree_light_95 <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), as.character(NA)))

dbhbin1995 <- with(alltree_light_95, logbin(x = dbh_corr, n = 20))

lightperareacloudbin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min, dbhbin1995$bin_max[numbins] + 1), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  group_modify(~ light_bin_stats(.$light_received/.$crownarea))

lightperareacloudbin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min, dbhbin1995$bin_max[numbins] + 1), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  group_modify(~ light_bin_stats(.$light_received/.$crownarea))

lightperareacloudbin_fg <- data.frame(fg = 'all', lightperareacloudbin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightperareacloudbin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

lightpervolcloudbin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min, dbhbin1995$bin_max[numbins] + 1), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  group_modify(~ light_bin_stats(.$light_received_byvolume))

lightpervolcloudbin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min, dbhbin1995$bin_max[numbins] + 1), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  group_modify(~ light_bin_stats(.$light_received_byvolume))

lightpervolcloudbin_fg <- data.frame(fg = 'all', lightpervolcloudbin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightpervolcloudbin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

unscaledlightbydbhcloudbin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min, dbhbin1995$bin_max[numbins] + 1), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  group_modify(~ light_bin_stats(.$light_received))

unscaledlightbydbhcloudbin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min, dbhbin1995$bin_max[numbins] + 1), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  group_modify(~ light_bin_stats(.$light_received))

unscaledlightbydbhcloudbin_fg <- data.frame(fg = 'all', unscaledlightbydbhcloudbin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(unscaledlightbydbhcloudbin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))


# Diameter growth plots ---------------------------------------------------


# Observed individual diameter growth
obs_indivdiamgrowth_df <- map(alltreedat[-1],
                              function(x) with(x %>% filter(!recruit), cloudbin_across_years(dat_values = diam_growth_rate, dat_classes = dbh_corr, edges = binedgedata, n_census = 1)))
obs_indivdiamgrowth_df <- cbind(fg = 'all', year = rep(years, each = nrow(binedgedata)), do.call(rbind, obs_indivdiamgrowth_df))

obs_indivdiamgrowth_fg <- replicate(length(fg_names), list())

for (i in 1:length(fg_names)) {
  for (j in 1:length(years)) {
    obs_indivdiamgrowth_fg[[i]][[j]] <-
      cbind(fg = fg_names[i], year = years[j],
            cloudbin_across_years(dat_values = fgdat[[i]][[j+1]] %>% filter(!recruit) %>% pull(diam_growth_rate),
                                  dat_classes = fgdat[[i]][[j+1]] %>% filter(!recruit) %>% pull(dbh_corr),
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

# Added 09 Jan 2020: Correct bad numbers of individuals in observed diameter growth by replacing with values from density
obs_indivdiamgrowth_df$mean_n_individuals <- joined_counts$bin_count



#-----------------------------------------

#john tomfoolery, --light per area calculations
# smallest 
dbhbin1995[1,]
lightperareacloudbin_all[1,]
alltree_light_95_filt <- alltree_light_95 %>%
  filter(dbh_corr >= 1.000000, dbh_corr < 1.310926) %>%
  mutate(median_light_received_byarea = median(light_received_byarea)) %>%
  select(dbh_corr, median_light_received_byarea)
alltree_light_95_filt[1:3,]
# median = 7.436764 
range(alltree_light_95$dbh_corr)

#second_smallest
dbhbin1995[2,]
lightperareacloudbin_all[2,]
alltree_light_95_filt2 <- alltree_light_95 %>%
  filter(dbh_corr >= 1.310926, dbh_corr < 1.718526) %>%
  mutate(median_light_received_byarea = median(light_received_byarea)) %>%
  select(dbh_corr, median_light_received_byarea)
alltree_light_95_filt2[1:3,]
#median = 8.81549

#third_smallest
dbhbin1995[3,]
lightperareacloudbin_all[2,]
alltree_light_95_filt3 <- alltree_light_95 %>%
  filter(dbh_corr >= 1.718526, dbh_corr < 2.25286) %>%
  mutate(median_light_received_byarea = median(light_received_byarea)) %>%
  select(dbh_corr, median_light_received_byarea)
alltree_light_95_filt3[1:3,]
# 12.19964

#7th_smallest
dbhbin1995[7,]
lightperareacloudbin_all[7,]
alltree_light_95_filt7<- alltree_light_95 %>%
  filter(dbh_corr >= 5.075376, dbh_corr < 6.65344) %>%
  mutate(median_light_received_byarea = median(light_received_byarea)) %>%
  select(dbh_corr, median_light_received_byarea)
alltree_light_95_filt7[1:3,]

dbhbin1995[12,]
lightperareacloudbin_all[12,]
alltree_light_95_filt12<- alltree_light_95 %>%
  filter(dbh_corr >= 19.64981, dbh_corr < 25.75944) %>%
  mutate(median_light_received_byarea = median(light_received_byarea)) %>%
  select(dbh_corr, median_light_received_byarea)
alltree_light_95_filt12[1:3,]
# 12.19964
#combine
dbhbin1995[16,]
alltree_light_95_filt16 <- alltree_light_95 %>%
  filter(dbh_corr >= 58.0324, dbh_corr < 76.07616) %>%
  mutate(median_light_received_byarea = median(light_received_byarea)) %>%
  select(dbh_corr, median_light_received_byarea)
alltree_light_95_filt16[1:3,]
# 12.19964
#combine
new_df <- data.frame(x = c(1.15546278812523, 1.514726, 1.985693, 5.86440808053007,22.70463,67.05428 ), 
                     y = c(7.436764,8.81549,  12.19964, 37.62251, 166.7423, 306.9545))

#plot light per area with bins, but don't bother with individual data
exl <- expression(atop('Light per Crown Area', paste('(W m'^-2, ')')))
exv <- expression(atop('Light per Crown Volume', paste('(W m'^-3, ')')))
exd <- 'Stem Diameter (cm)'

ggplot() +
   #geom_point(alpha = 0.01, data = alltree_light_95, 
    #                       aes(x = dbh_corr, y = light_received_byarea), color = 'chartreuse3') +
  geom_pointrange(data = lightperareacloudbin_fg %>% filter(fg %in% 'all', dbh_bin < 156), 
                  aes(x = dbh_bin, y = q50, ymin = mean, ymax = mean)) +
  geom_ribbon(data = fitted_lightcloudbin_fg %>% filter(fit == 'light per area', dbh < 156), 
              aes(x = dbh, ymin = q025, ymax = q975), alpha = 0.4) +
  geom_line(data = fitted_lightcloudbin_fg %>% filter(fit == 'light per area', dbh < 156),
            aes(x = dbh, y = q50)) +
  geom_point(data = new_df, aes(x = x, y = y), col = "red",shape = 21, fill = NA, size = 4) + 
  scale_x_log10(limits = c(0.8, 200), name = exd) +
  scale_y_log10(name = exl, limits = c(0.8, 1500)) +
  theme_plant()


# total light
#john tomfoolery
#smallest
dbhbin1995[1,]
lightperareacloudbin_all[1,]
alltree_light_95_filt_a <- alltree_light_95 %>%
  filter(dbh_corr >= 1.000000, dbh_corr <= 1.310926) %>%
  mutate(median_light_received = median(light_received)) %>%
  select(dbh_corr, median_light_received) 
alltree_light_95_filt_a[1:3,]
# median = 6.676022
range(alltree_light_95$dbh_corr)

#second_smallest
dbhbin1995[2,]
lightperareacloudbin_all[2,]
alltree_light_95_filt2a <- alltree_light_95 %>%
  filter(dbh_corr >= 1.310926, dbh_corr <= 1.718526) %>%
  mutate(median_light_received= median(light_received)) %>%
  select(dbh_corr, median_light_received) 
alltree_light_95_filt2a[1:3,]
#median = 11.55611

#third_smallest
dbhbin1995[3,]
lightperareacloudbin_all[2,]
alltree_light_95_filt3a <- alltree_light_95 %>%
  filter(dbh_corr >= 1.718526, dbh_corr <= 2.25286) %>%
  mutate(median_light_received= median(light_received)) %>%
  select(dbh_corr, median_light_received) 
alltree_light_95_filt3a[1:3,]
# 22.36945
#combine
new_df2 <- data.frame(x = c(1.15546278812523, 1.514726, 1.985693), y = c(6.676022, 11.55611,  22.36945))

hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000,10000), 
                                             labels = c(1,10,100,1000,10000), limits=c(1,10000))

alpha_value <- 1
hexfill2 <- scale_fill_gradient(low = 'forestgreen', high = 'navy', trans = 'log', breaks = c(1, 10,100))
exl <- expression(atop('Intercepted Light', paste('per Individual (W)')))
exd <- 'Stem Diameter (cm)'

labels = trans_format("log10", math_format(10^.x))

unscaledlightbydbhcloudbin_fg[1,]
lightperareacloudbin_fg[1,]
ggplot() +
  theme_plant() +
  scale_x_log10(name = exd, limits = c(.9, 200)) +
  scale_y_log10(name = exl, breaks = c(1,100,10000, 1000000), limits = c(1,1000000), 
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received)) +
  hex_scale_log_colors +
  geom_pointrange(data = unscaledlightbydbhcloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  geom_point(data = new_df2, aes(x = x, y = y), 
             col = "black", size = 4, size = 4, shape = 21, fill = "red", stroke = 2) + 
  theme_plant() + theme(legend.position = "right", legend.text = element_text(size = 15), legend.title = element_text(size = 16))+
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))# +



