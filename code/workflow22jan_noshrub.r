# Workflow binning and plotting, with shrubs removed
# 22 Jan
# Updated 24 Jan to do all the ratios too.

library(dplyr)

load('C:/Users/Q/google_drive/ForestLight/data/rawdataobj_22jan.r')
group_names <- c('all','all_classified','fg1','fg2','fg3','fg4','fg5','unclassified')

get_rid_of_shrubs <- function(x) {
  x[x$grform %in% c('M','T','U') | is.na(x$grform), ]
}

alltreedat <- lapply(alltreedat, get_rid_of_shrubs)
fgdat <- lapply(fgdat, function(z) lapply(z, get_rid_of_shrubs))

source('code/allfunctions27july.r')

# Set number of bins       
numbins <- 20

# Make a version of alltreedat without the unclassified trees
alltreedat_classified <- lapply(alltreedat, function(x) subset(x, !is.na(fg)))

# Bin classified trees. (log binning of density)
allyeardbh_classified <- unlist(lapply(alltreedat_classified[2:6], '[', , 'dbh_corr'))
dbhbin_allclassified <- logbin(x = allyeardbh_classified, y = NULL, n = numbins)

# Bin all trees including unclassified
allyeardbh <- unlist(lapply(alltreedat[2:6], '[', , 'dbh_corr'))
dbhbin_all <- logbin(x = allyeardbh, y = NULL, n = numbins)

allyeardbh_fg <- lapply(fgdat, function(x) unlist(lapply(x[2:6], '[', , 'dbh_corr')))

# The dbhbin_ objects can be used as the edges argument in logbin_setedges()

# Bin density and individual production by year using the above generated bin edges.

dbhbin_all_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_all))
dbhbin_allclassified_byyear <- lapply(alltreedat_classified[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_allclassified))

dbhbin_fg_byyear <- list()

for (i in 1:6) {
  dbhbin_fg_byyear[[i]] <- lapply(fgdat[[i]][2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_allclassified))
}

bin_across_years <- function(binlist) {
  binvals <- do.call('cbind', lapply(binlist, '[', , 'bin_value'))
  binindivs <- do.call('cbind', lapply(binlist, '[', , 'bin_count'))
  data.frame(bin_midpoint = binlist[[1]]$bin_midpoint,
             bin_min = binlist[[1]]$bin_min,
             bin_max = binlist[[1]]$bin_max,
             bin_yvalue = apply(binvals, 1, median),
             bin_ymin = apply(binvals, 1, min),
             bin_ymax = apply(binvals, 1, max),
             mean_n_individuals = apply(binindivs, 1, mean))
}

dbhbin_all_5census <- bin_across_years(dbhbin_all_byyear)
dbhbin_allclassified_5census <- bin_across_years(dbhbin_allclassified_byyear)
dbhbin_fg_5census <- lapply(dbhbin_fg_byyear, bin_across_years)

# Combine into single df
densitybin_5census <- cbind(fg = rep(group_names, each = numbins), 
                            rbind(dbhbin_all_5census, dbhbin_allclassified_5census, do.call('rbind', dbhbin_fg_5census)))

# Individual production 
# Take the mean and 2.5, 50 (median), 97.5 quantiles within each bin.
# Do it across all years.

allyearprod <- unlist(lapply(alltreedat[2:6], '[', , 'production'))
allyearprod_classified <- unlist(lapply(alltreedat_classified[2:6], '[', , 'production'))
allyearprod_fg <- lapply(fgdat, function(x) unlist(lapply(x[2:6], '[', , 'production')))

fakebin_across_years <- function(dat_values, dat_classes, edges, mean = 'geometric', n_census = 5) {
  qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  # add some padding just in case
  mins <- edges$bin_min
  mins[1] <- 0
  maxes <- edges$bin_max
  maxes[length(maxes)] <- Inf
  
  binstats <- t(sapply(1:length(mins), function(x) {
    indivs <- dat_values[dat_classes >= mins[x] & dat_classes < maxes[x]]
    c(mean = ifelse(mean == 'geometric', exp(mean(log(indivs))), mean(indivs)), quantile(indivs, probs = qprobs))
  }))
  dimnames(binstats)[[2]] <- c('mean', 'q025', 'q25', 'median', 'q75', 'q975')
  data.frame(bin_midpoint = edges$bin_midpoint,
             bin_min = edges$bin_min,
             bin_max = edges$bin_max,
             mean_n_individuals = edges$bin_count / n_census,
             binstats)
}

prodbin_all_5census <- fakebin_across_years(dat_values = allyearprod, dat_classes = allyeardbh, edges = dbhbin_all)
prodbin_allclassified_5census <- fakebin_across_years(dat_values = allyearprod_classified, dat_classes = allyeardbh_classified, edges = dbhbin_allclassified)

prodbin_fg_5census <- list()

for (i in 1:6) {
  prodbin_fg_5census[[i]] <- fakebin_across_years(dat_values = allyearprod_fg[[i]], dat_classes = allyeardbh_fg[[i]], edges = dbhbin_allclassified)
}

indivproductionbin_5census <- cbind(fg = rep(group_names, each = numbins),
                                    rbind(prodbin_all_5census, prodbin_allclassified_5census, do.call('rbind', prodbin_fg_5census)))

# Total production
# Do the binning for each year separately, as for density, then find min, max, and median.

totalprodbin_alltree_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_all))
totalprodbin_allclassified_byyear <- lapply(alltreedat_classified[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_allclassified))

totalprodbin_fg_byyear <- list()

for (i in 1:6) {
  totalprodbin_fg_byyear[[i]] <- lapply(fgdat[[i]][2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_allclassified))
}

totalprodbin_all_5census <- bin_across_years(totalprodbin_alltree_byyear)
totalprodbin_allclassified_5census <- bin_across_years(totalprodbin_allclassified_byyear)
totalprodbin_fg_5census <- lapply(totalprodbin_fg_byyear, bin_across_years)

# Combine into single df
totalproductionbin_5census <- cbind(fg = rep(group_names, each = numbins), 
                                    rbind(totalprodbin_all_5census, totalprodbin_allclassified_5census, do.call('rbind', totalprodbin_fg_5census)))

# Total light received and crown area
# 1990 and 1995 only

## crown area
crownareabin_alltree_byyear <- lapply(alltreedat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownarea, edges = dbhbin_all))
crownareabin_allclassified_byyear <- lapply(alltreedat_classified[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownarea, edges = dbhbin_allclassified))

crownareabin_fg_byyear <- list()

for (i in 1:6) {
  crownareabin_fg_byyear[[i]] <- lapply(fgdat[[i]][2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownarea, edges = dbhbin_allclassified))
}

crownareabin_all_2census <- bin_across_years(crownareabin_alltree_byyear)
crownareabin_allclassified_2census <- bin_across_years(crownareabin_allclassified_byyear)
crownareabin_fg_2census <- lapply(crownareabin_fg_byyear, bin_across_years)

crownareabin_2census <- cbind(fg = rep(group_names, each = numbins), 
                              rbind(crownareabin_all_2census, crownareabin_allclassified_2census, do.call('rbind', crownareabin_fg_2census)))

## light received
lightreceivedbin_alltree_byyear <- lapply(alltreedat[2:3], function(z) with(subset(z, !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_all)))
lightreceivedbin_allclassified_byyear <- lapply(alltreedat_classified[2:3], function(z) with(subset(z, !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_allclassified)))

lightreceivedbin_fg_byyear <- list()

for (i in 1:6) {
  lightreceivedbin_fg_byyear[[i]] <- lapply(fgdat[[i]][2:3], function(z) with(subset(z, !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_allclassified)))
}

lightreceivedbin_all_2census <- bin_across_years(lightreceivedbin_alltree_byyear)
lightreceivedbin_allclassified_2census <- bin_across_years(lightreceivedbin_allclassified_byyear)
lightreceivedbin_fg_2census <- lapply(lightreceivedbin_fg_byyear, bin_across_years)

lightreceivedbin_2census <- cbind(fg = rep(group_names, each = numbins), 
                                  rbind(lightreceivedbin_all_2census, lightreceivedbin_allclassified_2census, do.call('rbind', lightreceivedbin_fg_2census)))

# Production and light received per meter squared of crown area.
# Divide production by crown area and bin (1990 and 1995)
# Divide light received by crown area and bin (1990 and 1995)
# These are "fake" bins because it's just an individual measurement

binprod <- function(dat, bindat, xvar) {
  dat <- do.call('rbind', dat)
  dat$prod_area <- dat$production/dat$crownarea
  dat$light_area <- dat$light_received/dat$crownarea
  dat <- subset(dat, !is.na(light_received))
  with(dat, fakebin_across_years(dat_values = prod_area, dat_classes = light_area, edges = bindat, n_census = 2))
}

# Bin the entire light received per crown area dataset for 1990 and 1995 into a single set of bin edges.
light_per_area_all <- unlist(lapply(alltreedat[2:3], '[', , 'light_received'))/unlist(lapply(alltreedat[2:3], '[', , 'crownarea'))
light_per_area_allclassified <- unlist(lapply(alltreedat_classified[2:3], '[', , 'light_received'))/unlist(lapply(alltreedat_classified[2:3], '[', , 'crownarea'))
light_per_area_fg <- lapply(fgdat, function(x) {
  unlist(lapply(x[2:3], '[', , 'light_received'))/unlist(lapply(x[2:3], '[', , 'crownarea'))
})

light_per_area_bins_all <- logbin(x = na.omit(light_per_area_all), n = numbins)
light_per_area_bins_allclassified <- logbin(x = na.omit(light_per_area_allclassified), n = numbins)
light_per_area_bins_fg <- lapply(light_per_area_fg, function(z) logbin(x = na.omit(z), n = numbins))

indivprodperareabin_alltree_2census <- binprod(dat = alltreedat[2:3], bindat = light_per_area_bins_all, xvar = 'production')
indivprodperareabin_allclassified_2census <- binprod(dat = alltreedat_classified[2:3], bindat = light_per_area_bins_allclassified, xvar = 'production')
indivprodperareabin_fg_2census <- list()
for (i in 1:6) {
  indivprodperareabin_fg_2census[[i]] <- binprod(dat = fgdat[[i]][2:3], bindat = light_per_area_bins_fg[[i]], xvar = 'production')
}


indivprodperareabin_2census <- cbind(fg = rep(group_names, each = numbins), 
                                     rbind(indivprodperareabin_alltree_2census, indivprodperareabin_allclassified_2census, do.call('rbind', indivprodperareabin_fg_2census)))

# Figure 5 ratio binning
# Previously we had shade:gap ratio of both production and density by light received bin, as well as average shade score for each light bin.
# This was repeated for diameter bins.
# Here do the same 6 binning types, but repeat twice, once for ratio fg2:fg4 (longlived pioneer to shortlived breeder)
# and once for ratio fg1:fg3 (fast:slow)

binscore <- function(dat, bindat, score_column, class_column) {
  dat <- do.call('rbind', dat)
  dat <- dat[!is.na(dat$light_received) & !is.na(dat[,score_column]), ]
  dat$light_area <- dat$light_received/dat$crownarea
  with(dat, fakebin_across_years(dat_values = dat[,score_column], dat_classes = dat[,class_column], edges = bindat, mean = 'arithmetic', n_census = 2))
}


totalprodbin_byyear_bydiam <- lapply(fgdat, function(w) lapply(w[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_allclassified)))

densitybin_byyear_bydiam <- lapply(fgdat, function(w) lapply(w[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_allclassified)))

totalprodbin_byyear_bylight <- lapply(fgdat, function(w) lapply(w[2:3], function(z) logbin_setedges(x = z$light_received/z$crownarea, y = z$production, edges = light_per_area_bins_allclassified)))

densitybin_byyear_bylight <- lapply(fgdat, function(w) lapply(w[2:3], function(z) logbin_setedges(x = z$light_received/z$crownarea, y = NULL, edges = light_per_area_bins_allclassified)))

# Breeder to pioneer by diameter
breeder_stats_bydiam <- list() # fg4 to fg2

for (i in 1:2) {
  breeder_stats_bydiam[[i]] <- data.frame(bin = 1:numbins,
                                          year = c(1990,1995)[i],
                                          breeder_production_ratio = totalprodbin_byyear_bydiam[[4]][[i]]$bin_value/totalprodbin_byyear_bydiam[[2]][[i]]$bin_value,
                                          breeder_density_ratio = densitybin_byyear_bydiam[[4]][[i]]$bin_value/densitybin_byyear_bydiam[[2]][[i]]$bin_value)  
}

breeder_stats_bydiam <- do.call('rbind', breeder_stats_bydiam)

breeder_stats_bydiam[breeder_stats_bydiam == Inf | is.na(breeder_stats_bydiam)] <- NA

breeder_stats_bydiam_2census <- breeder_stats_bydiam %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(breeder_production_ratio),
            density_ratio_mean = mean(breeder_density_ratio),
            production_ratio_min = min(breeder_production_ratio),
            production_ratio_max = max(breeder_production_ratio),
            density_ratio_min = min(breeder_density_ratio),
            density_ratio_max = max(breeder_density_ratio)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) %>%
  cbind(mean_n_individuals = densitybin_byyear_bydiam[[2]][[1]]$bin_count + densitybin_byyear_bydiam[[2]][[2]]$bin_count + densitybin_byyear_bydiam[[4]][[1]]$bin_count + densitybin_byyear_bydiam[[4]][[2]]$bin_count / 2)

breederscore_bin_bydiam_2census <- binscore(dat = alltreedat[2:3], bindat = dbhbin_allclassified, score_column = 'X2', class_column = 'dbh_corr')

# Breeder to pioneer by light received per unit crown area
breeder_stats_bylight <- list() # fg2 to fg4

for (i in 1:2) {
  breeder_stats_bylight[[i]] <- data.frame(bin = 1:numbins,
                                           year = c(1990,1995)[i],
                                           breeder_production_ratio = totalprodbin_byyear_bylight[[4]][[i]]$bin_value/totalprodbin_byyear_bylight[[2]][[i]]$bin_value,
                                           breeder_density_ratio = densitybin_byyear_bylight[[4]][[i]]$bin_value/densitybin_byyear_bylight[[2]][[i]]$bin_value)  
}

breeder_stats_bylight <- do.call('rbind', breeder_stats_bylight)

breeder_stats_bylight[breeder_stats_bylight == Inf | is.na(breeder_stats_bylight)] <- NA

breeder_stats_bylight_2census <- breeder_stats_bylight %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(breeder_production_ratio),
            density_ratio_mean = mean(breeder_density_ratio),
            production_ratio_min = min(breeder_production_ratio),
            production_ratio_max = max(breeder_production_ratio),
            density_ratio_min = min(breeder_density_ratio),
            density_ratio_max = max(breeder_density_ratio)) %>%
  cbind(densitybin_byyear_bylight[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) %>%
  cbind(mean_n_individuals = densitybin_byyear_bylight[[2]][[1]]$bin_count + densitybin_byyear_bylight[[2]][[2]]$bin_count + densitybin_byyear_bylight[[4]][[1]]$bin_count + densitybin_byyear_bylight[[4]][[2]]$bin_count / 2)

breederscore_bin_bylight_2census <- binscore(dat = alltreedat[2:3], bindat = light_per_area_bins_allclassified, score_column = 'X2', class_column = 'light_area')


# Fast to slow by diameter

fastslow_stats_bydiam <- list() # fg1 to fg3

for (i in 1:2) {
  fastslow_stats_bydiam[[i]] <- data.frame(bin = 1:numbins,
                                           year = c(1990,1995)[i],
                                           fastslow_production_ratio = totalprodbin_byyear_bydiam[[1]][[i]]$bin_value/totalprodbin_byyear_bydiam[[3]][[i]]$bin_value,
                                           fastslow_density_ratio = densitybin_byyear_bydiam[[1]][[i]]$bin_value/densitybin_byyear_bydiam[[3]][[i]]$bin_value)  
}

fastslow_stats_bydiam <- do.call('rbind', fastslow_stats_bydiam)

fastslow_stats_bydiam[fastslow_stats_bydiam == Inf | is.na(fastslow_stats_bydiam)] <- NA

fastslow_stats_bydiam_2census <- fastslow_stats_bydiam %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(fastslow_production_ratio),
            density_ratio_mean = mean(fastslow_density_ratio),
            production_ratio_min = min(fastslow_production_ratio),
            production_ratio_max = max(fastslow_production_ratio),
            density_ratio_min = min(fastslow_density_ratio),
            density_ratio_max = max(fastslow_density_ratio)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) %>%
  cbind(mean_n_individuals = densitybin_byyear_bydiam[[1]][[1]]$bin_count + densitybin_byyear_bydiam[[1]][[2]]$bin_count + densitybin_byyear_bydiam[[3]][[1]]$bin_count + densitybin_byyear_bydiam[[3]][[2]]$bin_count / 2)

fastslowscore_bin_bydiam_2census <- binscore(dat = alltreedat[2:3], bindat = dbhbin_allclassified, score_column = 'X1', class_column = 'dbh_corr')


# Performer to recruiter by light received per unit crown area
fastslow_stats_bylight <- list() 

for (i in 1:2) {
  fastslow_stats_bylight[[i]] <- data.frame(bin = 1:numbins,
                                            year = c(1990,1995)[i],
                                            fastslow_production_ratio = totalprodbin_byyear_bylight[[1]][[i]]$bin_value/totalprodbin_byyear_bylight[[3]][[i]]$bin_value,
                                            fastslow_density_ratio = densitybin_byyear_bylight[[1]][[i]]$bin_value/densitybin_byyear_bylight[[3]][[i]]$bin_value)  
}

fastslow_stats_bylight <- do.call('rbind', fastslow_stats_bylight)

fastslow_stats_bylight[fastslow_stats_bylight == Inf | is.na(fastslow_stats_bylight)] <- NA

fastslow_stats_bylight_2census <- fastslow_stats_bylight %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(fastslow_production_ratio),
            density_ratio_mean = mean(fastslow_density_ratio),
            production_ratio_min = min(fastslow_production_ratio),
            production_ratio_max = max(fastslow_production_ratio),
            density_ratio_min = min(fastslow_density_ratio),
            density_ratio_max = max(fastslow_density_ratio)) %>%
  cbind(densitybin_byyear_bylight[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) %>%
  cbind(mean_n_individuals = densitybin_byyear_bylight[[1]][[1]]$bin_count + densitybin_byyear_bylight[[3]][[2]]$bin_count + densitybin_byyear_bylight[[1]][[1]]$bin_count + densitybin_byyear_bylight[[3]][[2]]$bin_count / 2)

fastslowscore_bin_bylight_2census <- binscore(dat = alltreedat[2:3], bindat = light_per_area_bins_allclassified, score_column = 'X1', class_column = 'light_area')

# Export binned data

fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'breeder_stats_bydiam_2census', 'breederscore_bin_bydiam_2census', 'breeder_stats_bylight_2census', 'breederscore_bin_bylight_2census', 'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census', 'fastslow_stats_bylight_2census', 'fastslowscore_bin_bylight_2census')

for (i in file_names) {
  write.csv(get(i), file=file.path(fpdata, paste0('noshrub_', i,'.csv')), row.names = FALSE)
}

save(list = file_names, file = file.path(fpdata, 'noshrub_bin_object.RData'))

# Load data ---------------------------------------------------------------

# Loop through all the csv files and load them into R
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'breeder_stats_bydiam_2census', 'breederscore_bin_bydiam_2census', 'breeder_stats_bylight_2census', 'breederscore_bin_bylight_2census', 'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census', 'fastslow_stats_bylight_2census', 'fastslowscore_bin_bylight_2census')

for (i in file_names) {
  assign(i, read.csv(file.path(fpdata, paste0('noshrub_', i,'.csv')), stringsAsFactors = FALSE))
}

# Plot data ---------------------------------------------------------------

library(cowplot)
library(dplyr)

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.

guild_colors <- RColorBrewer::brewer.pal(5, 'Set1')
fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('fast','long-lived pioneer', 'slow', 'short-lived breeder', 'intermediate')

# Figure 3a
fig_3a <- indivproductionbin_5census %>%
  filter(fg %in% fg_names) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Figure 3b
fig_3b <- densitybin_5census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Figure 3c
fig_3c <- totalproductionbin_5census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Total production (kg ha'^-1,' y'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Figure 4a
fig_4a <- crownareabin_2census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Total crown area (m'^2, ' ha'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Figure 4b
fig_4b <- lightreceivedbin_2census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Total light received (W ha'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Figure 4c
fig_4c <- indivprodperareabin_2census %>%
  filter(fg %in% fg_names) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = expression(paste('Production per unit crown area (kg y'^-1, ' m'^-2,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

pdf('C:/Users/Q/google_drive/ForestLight/figs/figures_22jan2018/noshrub_figs3and4.pdf', height = 6, width = 7.5)
print(fig_3a + ggtitle('3A'))
print(fig_3b + ggtitle('3B'))
print(fig_3c + ggtitle('3C'))
print(fig_4a + ggtitle('4A'))
print(fig_4b + ggtitle('4B'))
print(fig_4c + ggtitle('4C'))
dev.off()

###
# Added 23 Jan: Ratio figures (fig 5)
# Note that the continuous figures, as they don't use groupings, use all the data.

# Breeder to pioneer by light
fig_blightprod <- breeder_stats_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Breeder to pioneer production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_blightdens <- breeder_stats_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Breeder to pioneer density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_blightscore <- breederscore_bin_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_continuous(name = 'Low = long-lived pioneer, high = short-lived breeder') +
  panel_border(colour = 'black')

# Breeder to pioneer by diameter
fig_bdiamprod <- breeder_stats_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Breeder to pioneer production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_bdiamdens <- breeder_stats_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Breeder to pioneer density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_bdiamscore <- breederscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_continuous(name = 'Low = long-lived pioneer, high = short-lived breeder') +
  panel_border(colour = 'black')


####

# Fast to slow by light
fig_flightprod <- fastslow_stats_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Fast to slow production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_flightdens <- fastslow_stats_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Fast to slow density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_flightscore <- fastslowscore_bin_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_continuous(name = 'Low = slow, high = fast') +
  panel_border(colour = 'black')

# Fast to slow by diameter
fig_fdiamprod <- fastslow_stats_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Fast to slow production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_fdiamdens <- fastslow_stats_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Fast to slow density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_fdiamscore <- fastslowscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_continuous(name = 'Low = slow, high = fast') +
  panel_border(colour = 'black')

pdf('C:/Users/Q/google_drive/ForestLight/figs/figures_22jan2018/noshrub_ratio_figures.pdf', height = 6, width = 6)
print(fig_blightdens)
print(fig_blightprod)
print(fig_blightscore)
print(fig_bdiamdens)
print(fig_bdiamprod)
print(fig_bdiamscore)
print(fig_flightdens)
print(fig_flightprod)
print(fig_flightscore)
print(fig_fdiamdens)
print(fig_fdiamprod)
print(fig_fdiamscore)
dev.off()