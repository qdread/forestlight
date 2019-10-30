# Script 3: Create log-bins for all variables.

# Script split off on 30 October 2019.
# Crown volume added on 21 March 2019.

# Procedure: use all 5 FGs (but not unclassified) together to get the bin edges
# Then apply those bin edges to each FG in isolation.
# That way each FG are given the same bin edges.

library(tidyverse)

load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')
group_names <- c('all','all_classified','fg1','fg2','fg3','fg4','fg5','unclassified')
source('code/allfunctions27july.r')

# Set number of bins       
numbins <- 20

# Make a version of alltreedat without the unclassified trees
alltreedat_classified <- map(alltreedat, ~ filter(., !is.na(fg)))

# Bin classified trees. (log binning of density)
allyeardbh_classified <- map(alltreedat_classified[-1], ~ pull(., dbh_corr)) %>% unlist
dbhbin_allclassified <- logbin(x = allyeardbh_classified, y = NULL, n = numbins)

# Bin all trees including unclassified
allyeardbh <- map(alltreedat[-1], ~ pull(., dbh_corr)) %>% unlist
dbhbin_all <- logbin(x = allyeardbh, y = NULL, n = numbins)

allyeardbh_fg <- fgdat %>% map(~ map(., ~ pull(., dbh_corr)) %>% unlist)

# The dbhbin_ objects can be used as the edges argument in logbin_setedges()

# Bin density and individual production by year using the above generated bin edges.

dbhbin_all_byyear <- alltreedat[-1] %>% map(~ logbin_setedges(x = .$dbh_corr, y = NULL, edges = dbhbin_all))
dbhbin_allclassified_byyear <- alltreedat_classified[-1] %>% map(~ logbin_setedges(x = .$dbh_corr, y = NULL, edges = dbhbin_allclassified))

dbhbin_fg_byyear <- fgdat %>%
  map(~ map(.[-1], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_allclassified)))

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
allyearprod <- map(alltreedat[-1], ~ pull(., production)) %>% unlist
allyearprod_classified <- map(alltreedat_classified[-1], ~ pull(., production)) %>% unlist
allyearprod_fg <- fgdat %>% map(~ map(., ~ pull(., production)) %>% unlist)

fakebin_across_years <- function(dat_values, dat_classes, edges, mean_type = 'geometric', n_census = 5) {
  qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  # add some padding just in case
  mins <- edges$bin_min
  mins[1] <- 0
  maxes <- edges$bin_max
  maxes[length(maxes)] <- Inf
  
  binstats <- t(sapply(1:length(mins), function(x) {
    indivs <- dat_values[dat_classes >= mins[x] & dat_classes < maxes[x]]
    if (mean_type == 'geometric') {
      mean_n <- exp(mean(log(indivs)))
      sd_n <- exp(sd(log(indivs)))
      ci_width <- qnorm(0.975) * sd(log(indivs)) / sqrt(length(indivs))
      ci_min <- exp(mean(log(indivs)) - ci_width)
      ci_max <- exp(mean(log(indivs)) + ci_width)
      
    } else {
      mean_n <- mean(indivs)
      sd_n <- sd(indivs)
      ci_width <- qnorm(0.975) * sd(indivs) / sqrt(length(indivs))
      ci_min <- mean_n - ci_width
      ci_max <- mean_n + ci_width
    }
    c(mean = mean_n, 
      sd = sd_n,
      quantile(indivs, probs = qprobs),
      ci_min = ci_min,
      ci_max = ci_max)
  }))
  dimnames(binstats)[[2]] <- c('mean', 'sd', 'q025', 'q25', 'median', 'q75', 'q975', 'ci_min', 'ci_max')
  data.frame(bin_midpoint = edges$bin_midpoint,
             bin_min = edges$bin_min,
             bin_max = edges$bin_max,
             mean_n_individuals = edges$bin_count / n_census,
             binstats)
}

prodbin_all_5census <- fakebin_across_years(dat_values = allyearprod, dat_classes = allyeardbh, edges = dbhbin_all)
prodbin_allclassified_5census <- fakebin_across_years(dat_values = allyearprod_classified, dat_classes = allyeardbh_classified, edges = dbhbin_allclassified)

prodbin_fg_5census <- map2(allyearprod_fg, allyeardbh_fg, ~ fakebin_across_years(dat_values = .x, dat_classes = .y, edges = dbhbin_allclassified))

indivproductionbin_5census <- cbind(fg = rep(group_names, each = numbins),
                                    rbind(prodbin_all_5census, prodbin_allclassified_5census, do.call('rbind', prodbin_fg_5census)))

# Total production
# Do the binning for each year separately, as for density, then find min, max, and median.

totalprodbin_alltree_byyear <- alltreedat[-1] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$production, edges = dbhbin_all))
totalprodbin_allclassified_byyear <- alltreedat_classified[-1] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$production, edges = dbhbin_allclassified))

totalprodbin_fg_byyear <- fgdat %>%
  map(~ map(.[-1], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_allclassified)))

totalprodbin_all_5census <- bin_across_years(totalprodbin_alltree_byyear)
totalprodbin_allclassified_5census <- bin_across_years(totalprodbin_allclassified_byyear)
totalprodbin_fg_5census <- lapply(totalprodbin_fg_byyear, bin_across_years)

# Combine into single df
totalproductionbin_5census <- cbind(fg = rep(group_names, each = numbins), 
                                    rbind(totalprodbin_all_5census, totalprodbin_allclassified_5census, do.call('rbind', totalprodbin_fg_5census)))

# Total light received and crown area
# 1990 and 1995 only

## crown area
crownareabin_alltree_byyear <- alltreedat[2:3] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$crownarea, edges = dbhbin_all))
crownareabin_allclassified_byyear <- alltreedat_classified[2:3] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$crownarea, edges = dbhbin_allclassified))

crownareabin_fg_byyear <- fgdat %>%
  map(~ map(.[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownarea, edges = dbhbin_allclassified)))

crownareabin_all_2census <- bin_across_years(crownareabin_alltree_byyear)
crownareabin_allclassified_2census <- bin_across_years(crownareabin_allclassified_byyear)
crownareabin_fg_2census <- lapply(crownareabin_fg_byyear, bin_across_years)

crownareabin_2census <- cbind(fg = rep(group_names, each = numbins), 
                              rbind(crownareabin_all_2census, crownareabin_allclassified_2census, do.call('rbind', crownareabin_fg_2census)))

## crown volume
crownvolumebin_alltree_byyear <- alltreedat[2:3] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$crownvolume, edges = dbhbin_all))
crownvolumebin_allclassified_byyear <- alltreedat_classified[2:3] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$crownvolume, edges = dbhbin_allclassified))

crownvolumebin_fg_byyear <- fgdat %>%
  map(~ map(.[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownvolume, edges = dbhbin_allclassified)))

crownvolumebin_all_2census <- bin_across_years(crownvolumebin_alltree_byyear)
crownvolumebin_allclassified_2census <- bin_across_years(crownvolumebin_allclassified_byyear)
crownvolumebin_fg_2census <- lapply(crownvolumebin_fg_byyear, bin_across_years)

crownvolumebin_2census <- cbind(fg = rep(group_names, each = numbins), 
                                rbind(crownvolumebin_all_2census, crownvolumebin_allclassified_2census, do.call('rbind', crownvolumebin_fg_2census)))

## light received
lightreceivedbin_alltree_byyear <- alltreedat[2:3] %>% map(~ with(filter(., !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_all)))
lightreceivedbin_allclassified_byyear <- alltreedat_classified[2:3] %>% map(~ with(filter(., !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_all)))

lightreceivedbin_fg_byyear <- fgdat %>%
  map(~ map(.[2:3], function(z) with(filter(z, !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_allclassified))))

lightreceivedbin_all_2census <- bin_across_years(lightreceivedbin_alltree_byyear)
lightreceivedbin_allclassified_2census <- bin_across_years(lightreceivedbin_allclassified_byyear)
lightreceivedbin_fg_2census <- lapply(lightreceivedbin_fg_byyear, bin_across_years)

lightreceivedbin_2census <- cbind(fg = rep(group_names, each = numbins), 
                                  rbind(lightreceivedbin_all_2census, lightreceivedbin_allclassified_2census, do.call('rbind', lightreceivedbin_fg_2census)))


## light received per unit crown area
lightperareabin_alltree_byyear <- alltreedat[2:3] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$light_received/.$crownarea, edges = dbhbin_all))
lightperareabin_allclassified_byyear <- alltreedat_classified[2:3] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$light_received/.$crownarea, edges = dbhbin_all))

lightperareabin_fg_byyear <- fgdat %>%
  map(~ map(.[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$light_received/z$crownarea, edges = dbhbin_allclassified)))

## light received per unit crown volume
lightpervolumebin_alltree_byyear <- alltreedat[2:3] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$light_received/.$crownvolume, edges = dbhbin_all))
lightpervolumebin_allclassified_byyear <- alltreedat_classified[2:3] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$light_received/.$crownvolume, edges = dbhbin_all))

lightpervolumebin_fg_byyear <- fgdat %>%
  map(~ map(.[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$light_received/z$crownvolume, edges = dbhbin_allclassified)))

## Combine the individual 1995 bins to data frames and then write them to R object.
crownareabins1995 <- rbind(data.frame(year = 1995, fg = 'all', crownareabin_allclassified_byyear[[2]]),
                           map2_dfr(crownareabin_fg_byyear, group_names[3:8], ~ data.frame(year = 1995, fg = .y, .x[[2]])))
crownvolumebins1995 <- rbind(data.frame(year = 1995, fg = 'all', crownvolumebin_allclassified_byyear[[2]]),
                             map2_dfr(crownvolumebin_fg_byyear, group_names[3:8], ~ data.frame(year = 1995, fg = .y, .x[[2]])))
lightreceivedbins1995 <- rbind(data.frame(year = 1995, fg = 'all', lightreceivedbin_allclassified_byyear[[2]]),
                               map2_dfr(lightreceivedbin_fg_byyear, group_names[3:8], ~ data.frame(year = 1995, fg = .y, .x[[2]])))
lightperareabins1995 <- rbind(data.frame(year = 1995, fg = 'all', lightperareabin_allclassified_byyear[[2]]),
                              map2_dfr(lightperareabin_fg_byyear, group_names[3:8], ~ data.frame(year = 1995, fg = .y, .x[[2]])))
lightpervolumebins1995 <- rbind(data.frame(year = 1995, fg = 'all', lightpervolumebin_allclassified_byyear[[2]]),
                                map2_dfr(lightpervolumebin_fg_byyear, group_names[3:8], ~ data.frame(year = 1995, fg = .y, .x[[2]])))

# Production and light received per meter squared of crown area.
# Divide production by crown area and bin (1990 and 1995)
# Divide light received by crown area and bin (1990 and 1995)
# These are "fake" bins because it's just an individual measurement

binprod <- function(dat, bindat) {
  dat <- do.call('rbind', dat)
  dat$prod_area <- dat$production/dat$crownarea
  dat$light_area <- dat$light_received/dat$crownarea
  dat <- subset(dat, !is.na(light_received))
  with(dat, fakebin_across_years(dat_values = prod_area, dat_classes = light_area, edges = bindat, n_census = 2))
}

# Bin the entire light received per crown area dataset for 1990 and 1995 into a single set of bin edges.
light_per_area_all <- alltreedat[2:3] %>% map(~ .$light_received/.$crownarea) %>% unlist
light_per_area_allclassified <- alltreedat_classified[2:3] %>% map(~ .$light_received/.$crownarea) %>% unlist

light_per_area_fg <- fgdat %>%
  map(~ map(.[2:3], function(z) z$light_received/z$crownarea) %>% unlist)

light_per_area_bins_all <- logbin(x = na.omit(light_per_area_all), n = numbins)
light_per_area_bins_allclassified <- logbin(x = na.omit(light_per_area_allclassified), n = numbins)
light_per_area_bins_fg <- map(light_per_area_fg, ~ logbin(x = na.omit(.), n = numbins))

indivprodperareabin_alltree_2census <- binprod(dat = alltreedat[2:3], bindat = light_per_area_bins_all)
indivprodperareabin_allclassified_2census <- binprod(dat = alltreedat_classified[2:3], bindat = light_per_area_bins_allclassified)
indivprodperareabin_fg_2census <- map2(fgdat, light_per_area_bins_fg, ~ binprod(dat = .x[2:3], bindat = .y))

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

totalprodbin_byyear_bydiam <- fgdat %>%
  map(~ map(.[-1], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_allclassified)))

densitybin_byyear_bydiam <- fgdat %>%
  map(~ map(.[-1], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_allclassified)))

totalprodbin_byyear_bylight <- fgdat %>%
  map(~ map(.[2:3], function(z) logbin_setedges(x = z$light_received/z$crownarea, y = z$production, edges = light_per_area_bins_allclassified)))

densitybin_byyear_bylight <- fgdat %>%
  map(~ map(.[2:3], function(z) logbin_setedges(x = z$light_received/z$crownarea, y = NULL, edges = light_per_area_bins_allclassified)))

# Breeder to pioneer by diameter (fg4 to fg2)
breeder_stats_bydiam <- tibble(fg_a_prod = totalprodbin_byyear_bydiam[[4]],
                               fg_b_prod = totalprodbin_byyear_bydiam[[2]],
                               fg_a_dens = densitybin_byyear_bydiam[[4]],
                               fg_b_dens = densitybin_byyear_bydiam[[2]],
                               year = c(1990,1995,2000,2005,2010)) %>%
  pmap(function(fg_a_prod, fg_b_prod, fg_a_dens, fg_b_dens, year) 
    data.frame(bin = 1:numbins,
               year = year,
               breeder_production_ratio = fg_a_prod$bin_value / fg_b_prod$bin_value,
               breeder_density_ratio = fg_a_dens$bin_value / fg_b_dens$bin_value,
               n_individuals = fg_a_dens$bin_count + fg_b_dens$bin_count)) %>%
  bind_rows %>%
  mutate_if(is.double, ~ if_else(is.finite(.x), .x, as.numeric(NA)))

breeder_stats_bydiam_2census <- breeder_stats_bydiam %>% 
  filter(year %in% c(1990, 1995)) %>%
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(breeder_production_ratio),
            density_ratio_mean = mean(breeder_density_ratio),
            production_ratio_min = min(breeder_production_ratio),
            production_ratio_max = max(breeder_production_ratio),
            density_ratio_min = min(breeder_density_ratio),
            density_ratio_max = max(breeder_density_ratio),
            mean_n_individuals = mean(n_individuals)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 

breederscore_bin_bydiam_2census <- binscore(dat = alltreedat[2:3], bindat = dbhbin_allclassified, score_column = 'X2', class_column = 'dbh_corr')

breeder_stats_bydiam_5census <- breeder_stats_bydiam %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(breeder_production_ratio),
            density_ratio_mean = mean(breeder_density_ratio),
            production_ratio_min = min(breeder_production_ratio),
            production_ratio_max = max(breeder_production_ratio),
            density_ratio_min = min(breeder_density_ratio),
            density_ratio_max = max(breeder_density_ratio),
            mean_n_individuals = mean(n_individuals)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 


# Breeder to pioneer by light received per unit crown area (fg4 to fg2)
breeder_stats_bylight <- tibble(fg_a_prod = totalprodbin_byyear_bylight[[4]],
                                fg_b_prod = totalprodbin_byyear_bylight[[2]],
                                fg_a_dens = densitybin_byyear_bylight[[4]],
                                fg_b_dens = densitybin_byyear_bylight[[2]],
                                year = c(1990,1995)) %>%
  pmap(function(fg_a_prod, fg_b_prod, fg_a_dens, fg_b_dens, year) 
    data.frame(bin = 1:numbins,
               year = year,
               breeder_production_ratio = fg_a_prod$bin_value / fg_b_prod$bin_value,
               breeder_density_ratio = fg_a_dens$bin_value / fg_b_dens$bin_value,
               n_individuals = fg_a_dens$bin_count + fg_b_dens$bin_count)) %>%
  bind_rows %>%
  mutate_if(is.double, ~ if_else(is.finite(.x), .x, as.numeric(NA)))

breeder_stats_bylight_2census <- breeder_stats_bylight %>% 
  filter(year %in% c(1990, 1995)) %>%
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(breeder_production_ratio),
            density_ratio_mean = mean(breeder_density_ratio),
            production_ratio_min = min(breeder_production_ratio),
            production_ratio_max = max(breeder_production_ratio),
            density_ratio_min = min(breeder_density_ratio),
            density_ratio_max = max(breeder_density_ratio),
            mean_n_individuals = mean(n_individuals)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 

breederscore_bin_bylight_2census <- binscore(dat = alltreedat[2:3], bindat = dbhbin_allclassified, score_column = 'X2', class_column = 'light_area')

# Fast to slow by diameter (fg1 to fg3)
fastslow_stats_bydiam <- tibble(fg_a_prod = totalprodbin_byyear_bydiam[[1]],
                                fg_b_prod = totalprodbin_byyear_bydiam[[3]],
                                fg_a_dens = densitybin_byyear_bydiam[[1]],
                                fg_b_dens = densitybin_byyear_bydiam[[3]],
                                year = c(1990,1995,2000,2005,2010)) %>%
  pmap(function(fg_a_prod, fg_b_prod, fg_a_dens, fg_b_dens, year) 
    data.frame(bin = 1:numbins,
               year = year,
               fastslow_production_ratio = fg_a_prod$bin_value / fg_b_prod$bin_value,
               fastslow_density_ratio = fg_a_dens$bin_value / fg_b_dens$bin_value,
               n_individuals = fg_a_dens$bin_count + fg_b_dens$bin_count)) %>%
  bind_rows %>%
  mutate_if(is.double, ~ if_else(is.finite(.x), .x, as.numeric(NA)))

fastslow_stats_bydiam_2census <- fastslow_stats_bydiam %>% 
  filter(year %in% c(1990, 1995)) %>%
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(fastslow_production_ratio),
            density_ratio_mean = mean(fastslow_density_ratio),
            production_ratio_min = min(fastslow_production_ratio),
            production_ratio_max = max(fastslow_production_ratio),
            density_ratio_min = min(fastslow_density_ratio),
            density_ratio_max = max(fastslow_density_ratio),
            mean_n_individuals = mean(n_individuals)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 

fastslowscore_bin_bydiam_2census <- binscore(dat = alltreedat[2:3], bindat = dbhbin_allclassified, score_column = 'X1', class_column = 'dbh_corr')

fastslow_stats_bydiam_5census <- fastslow_stats_bydiam %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(fastslow_production_ratio),
            density_ratio_mean = mean(fastslow_density_ratio),
            production_ratio_min = min(fastslow_production_ratio),
            production_ratio_max = max(fastslow_production_ratio),
            density_ratio_min = min(fastslow_density_ratio),
            density_ratio_max = max(fastslow_density_ratio),
            mean_n_individuals = mean(n_individuals)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 

# Fast to slow by light received per unit crown area (fg1 to fg3)
fastslow_stats_bylight <- tibble(fg_a_prod = totalprodbin_byyear_bylight[[1]],
                                 fg_b_prod = totalprodbin_byyear_bylight[[3]],
                                 fg_a_dens = densitybin_byyear_bylight[[1]],
                                 fg_b_dens = densitybin_byyear_bylight[[3]],
                                 year = c(1990,1995)) %>%
  pmap(function(fg_a_prod, fg_b_prod, fg_a_dens, fg_b_dens, year) 
    data.frame(bin = 1:numbins,
               year = year,
               fastslow_production_ratio = fg_a_prod$bin_value / fg_b_prod$bin_value,
               fastslow_density_ratio = fg_a_dens$bin_value / fg_b_dens$bin_value,
               n_individuals = fg_a_dens$bin_count + fg_b_dens$bin_count)) %>%
  bind_rows %>%
  mutate_if(is.double, ~ if_else(is.finite(.x), .x, as.numeric(NA)))

fastslow_stats_bylight_2census <- fastslow_stats_bylight %>% 
  filter(year %in% c(1990, 1995)) %>%
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(fastslow_production_ratio),
            density_ratio_mean = mean(fastslow_density_ratio),
            production_ratio_min = min(fastslow_production_ratio),
            production_ratio_max = max(fastslow_production_ratio),
            density_ratio_min = min(fastslow_density_ratio),
            density_ratio_max = max(fastslow_density_ratio),
            mean_n_individuals = mean(n_individuals)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 

fastslowscore_bin_bylight_2census <- binscore(dat = alltreedat[2:3], bindat = dbhbin_allclassified, score_column = 'X1', class_column = 'light_area')

# Export binned data

fpdata <- '~/google_drive/ForestLight/data/data_binned'
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'breeder_stats_bydiam_2census',  'breederscore_bin_bydiam_2census', 'breeder_stats_bylight_2census', 'breederscore_bin_bylight_2census', 'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census', 'fastslow_stats_bylight_2census', 'fastslowscore_bin_bylight_2census','fastslow_stats_bydiam_5census','breeder_stats_bydiam_5census')

for (i in file_names) {
  write.csv(get(i), file=file.path(fpdata, paste0(i,'.csv')), row.names = FALSE)
}

save(list = file_names, file = file.path(fpdata, 'bin_object.RData'))

save(crownareabins1995, crownvolumebins1995, lightreceivedbins1995, lightperareabins1995, lightpervolumebins1995,
     file = '~/google_drive/ForestLight/data/data_binned/area_and_volume_bins_1995.RData')
