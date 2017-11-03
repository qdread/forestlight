# Binning and error bars: all years combined ------------------------------

# Script copied and modified on 1 Nov 2017. Now use shade and gap (but not unclassified) together to get the bin edges
# Then apply those bin edges to the shade and gap in isolation.
# That way both shade and gap are given the same bin edges.

# Edited 3 Nov: add number of individuals as a column, also do geometric mean and add the central 50% quantiles

# Concatenate the 5 censuses (3-7) excluding 1985, and find log bin edges.
# Keep the size classes constant across density and production.

# All years together are used to create the 20 bin boundaries, then the binning is done separately for each year.
# Use the bounds from the dbh bin to create the bins for all other scalings.

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_04oct'
load(file.path(fpdata, 'rawdataobj.r'))

# Source code for log binning function
source(file.path(fpdata, 'allfunctions27july.r'))
library(dplyr)

# Set number of bins       
numbins <- 20

# Combine shade and gap (not unclassified) into single data frame.
shadegapdat <- lapply(alltreedat, function(x) subset(x, tol_wright %in% c('S','G')))

# Log bin density
allyeardbh_alltree <- unlist(lapply(alltreedat[2:6], '[', , 'dbh_corr'))
allyeardbh_shade <- unlist(lapply(shadedat[2:6], '[', , 'dbh_corr'))
allyeardbh_gap <- unlist(lapply(gapdat[2:6], '[', , 'dbh_corr'))
allyeardbh_unclassified <- unlist(lapply(unclassifieddat[2:6], '[', , 'dbh_corr'))

dbhbin_alltree <- logbin(x = allyeardbh_alltree, y = NULL, n = numbins)
dbhbin_shade <- logbin(x = allyeardbh_shade, y = NULL, n = numbins)
dbhbin_gap <- logbin(x = allyeardbh_gap, y = NULL, n = numbins)
dbhbin_unclassified <- logbin(x = allyeardbh_unclassified, y = NULL, n = numbins)

# Bin shade and gap.
allyeardbh_shadegap <- unlist(lapply(shadegapdat[2:6], '[', , 'dbh_corr'))
dbhbin_shadegap <- logbin(x = allyeardbh_shadegap, y = NULL, n = numbins)


# The dbhbin_ objects can be used as the edges argument in logbin_setedges()

# Bin density and individual production by year using the above generated bin edges.

dbhbin_alltree_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_alltree))
dbhbin_shade_byyear <- lapply(shadedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_shadegap))
dbhbin_gap_byyear <- lapply(gapdat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_shadegap))
dbhbin_unclassified_byyear <- lapply(unclassifieddat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_unclassified))

# Find median, minimum and maximum bin values across years.

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

dbhbin_alltree_5census <- bin_across_years(dbhbin_alltree_byyear)
dbhbin_shade_5census <- bin_across_years(dbhbin_shade_byyear)
dbhbin_gap_5census <- bin_across_years(dbhbin_gap_byyear)
dbhbin_unclassified_5census <- bin_across_years(dbhbin_unclassified_byyear)

# Combine into single df
densitybin_5census <- rbind(cbind(guild = 'all', dbhbin_alltree_5census),
                            cbind(guild = 'shade', dbhbin_shade_5census),
                            cbind(guild = 'gap', dbhbin_gap_5census),
                            cbind(guild = 'unclassified', dbhbin_unclassified_5census))

# Individual production 
# Take the mean and 2.5, 50 (median), 97.5 quantiles within each bin.
# Do it across all years.

allyearprod_alltree <- unlist(lapply(alltreedat[2:6], '[', , 'production'))
allyearprod_shade <- unlist(lapply(shadedat[2:6], '[', , 'production'))
allyearprod_gap <- unlist(lapply(gapdat[2:6], '[', , 'production'))
allyearprod_unclassified <- unlist(lapply(unclassifieddat[2:6], '[', , 'production'))

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

prodbin_alltree_5census <- fakebin_across_years(dat_values = allyearprod_alltree, dat_classes = allyeardbh_alltree, edges = dbhbin_alltree)
prodbin_shade_5census <- fakebin_across_years(dat_values = allyearprod_shade, dat_classes = allyeardbh_shade, edges = dbhbin_shadegap)
prodbin_gap_5census <- fakebin_across_years(dat_values = allyearprod_gap, dat_classes = allyeardbh_gap, edges = dbhbin_shadegap)
prodbin_unclassified_5census <- fakebin_across_years(dat_values = allyearprod_unclassified, dat_classes = allyeardbh_unclassified, edges = dbhbin_unclassified)

indivproductionbin_5census <- rbind(cbind(guild = 'all', prodbin_alltree_5census),
                                    cbind(guild = 'shade', prodbin_shade_5census),
                                    cbind(guild = 'gap', prodbin_gap_5census),
                                    cbind(guild = 'unclassified', prodbin_unclassified_5census))

# Total production
# Do the binning for each year separately, as for density, then find min, max, and median.
totalprodbin_alltree_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_alltree))
totalprodbin_shade_byyear <- lapply(shadedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_shadegap))
totalprodbin_gap_byyear <- lapply(gapdat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_shadegap))
totalprodbin_unclassified_byyear <- lapply(unclassifieddat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_unclassified))

totalprodbin_alltree_5census <- bin_across_years(totalprodbin_alltree_byyear)
totalprodbin_shade_5census <- bin_across_years(totalprodbin_shade_byyear)
totalprodbin_gap_5census <- bin_across_years(totalprodbin_gap_byyear)
totalprodbin_unclassified_5census <- bin_across_years(totalprodbin_unclassified_byyear)

totalproductionbin_5census <- rbind(cbind(guild = 'all', totalprodbin_alltree_5census),
                                    cbind(guild = 'shade', totalprodbin_shade_5census),
                                    cbind(guild = 'gap', totalprodbin_gap_5census),
                                    cbind(guild = 'unclassified', totalprodbin_unclassified_5census))

# Total light received and crown area
# 1990 and 1995 only

crownareabin_alltree_byyear <- lapply(alltreedat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownarea, edges = dbhbin_alltree))
crownareabin_shade_byyear <- lapply(shadedat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownarea, edges = dbhbin_shadegap))
crownareabin_gap_byyear <- lapply(gapdat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownarea, edges = dbhbin_shadegap))
crownareabin_unclassified_byyear <- lapply(unclassifieddat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownarea, edges = dbhbin_unclassified))

crownareabin_alltree_2census <- bin_across_years(crownareabin_alltree_byyear)
crownareabin_shade_2census <- bin_across_years(crownareabin_shade_byyear)
crownareabin_gap_2census <- bin_across_years(crownareabin_gap_byyear)
crownareabin_unclassified_2census <- bin_across_years(crownareabin_unclassified_byyear)

crownareabin_2census <- rbind(cbind(guild = 'all', crownareabin_alltree_2census),
                              cbind(guild = 'shade', crownareabin_shade_2census),
                              cbind(guild = 'gap', crownareabin_gap_2census),
                              cbind(guild = 'unclassified', crownareabin_unclassified_2census))

lightreceivedbin_alltree_byyear <- lapply(alltreedat[2:3], function(z) with(subset(z, !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_alltree)))
lightreceivedbin_shade_byyear <- lapply(shadedat[2:3], function(z) with(subset(z, !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_shadegap)))
lightreceivedbin_gap_byyear <- lapply(gapdat[2:3], function(z) with(subset(z, !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_shadegap)))
lightreceivedbin_unclassified_byyear <- lapply(unclassifieddat[2:3], function(z) with(subset(z, !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_alltree)))

lightreceivedbin_alltree_2census <- bin_across_years(lightreceivedbin_alltree_byyear)
lightreceivedbin_shade_2census <- bin_across_years(lightreceivedbin_shade_byyear)
lightreceivedbin_gap_2census <- bin_across_years(lightreceivedbin_gap_byyear)
lightreceivedbin_unclassified_2census <- bin_across_years(lightreceivedbin_unclassified_byyear)

lightreceivedbin_2census <- rbind(cbind(guild = 'all', lightreceivedbin_alltree_2census),
                                  cbind(guild = 'shade', lightreceivedbin_shade_2census),
                                  cbind(guild = 'gap', lightreceivedbin_gap_2census),
                                  cbind(guild = 'unclassified', lightreceivedbin_unclassified_2census))

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
light_per_area_shade <- unlist(lapply(shadedat[2:3], '[', , 'light_received'))/unlist(lapply(shadedat[2:3], '[', , 'crownarea'))
light_per_area_gap <- unlist(lapply(gapdat[2:3], '[', , 'light_received'))/unlist(lapply(gapdat[2:3], '[', , 'crownarea'))
light_per_area_unclassified <- unlist(lapply(unclassifieddat[2:3], '[', , 'light_received'))/unlist(lapply(unclassifieddat[2:3], '[', , 'crownarea'))
light_per_area_shadegap <- unlist(lapply(shadegapdat[2:3], '[', , 'light_received'))/unlist(lapply(shadegapdat[2:3], '[', , 'crownarea'))


light_per_area_bins_all <- logbin(x = na.omit(light_per_area_all), n = numbins)
light_per_area_bins_shade <- logbin(x = na.omit(light_per_area_shade), n = numbins)
light_per_area_bins_gap <- logbin(x = na.omit(light_per_area_gap), n = numbins)
light_per_area_bins_unclassified <- logbin(x = na.omit(light_per_area_unclassified), n = numbins)
light_per_area_bins_shadegap <- logbin(x = na.omit(light_per_area_shadegap), n = numbins)

indivprodperareabin_alltree_2census <- binprod(dat = alltreedat[2:3], bindat = light_per_area_bins_all, xvar = 'production')
indivprodperareabin_shade_2census <- binprod(dat = shadedat[2:3], bindat = light_per_area_bins_shadegap, xvar = 'production')
indivprodperareabin_gap_2census <- binprod(dat = gapdat[2:3], bindat = light_per_area_bins_shadegap, xvar = 'production')
indivprodperareabin_unclassified_2census <- binprod(dat = unclassifieddat[2:3], bindat = light_per_area_bins_unclassified, xvar = 'production')

indivprodperareabin_2census <- rbind(cbind(guild = 'all', indivprodperareabin_alltree_2census),
                                     cbind(guild = 'shade', indivprodperareabin_shade_2census),
                                     cbind(guild = 'gap', indivprodperareabin_gap_2census),
                                     cbind(guild = 'unclassified', indivprodperareabin_unclassified_2census))


# Figure 5: need to bin (a) the number of shade individuals over the number of gap individuals, (b) total shade production over total gap production, and (c) average shade tolerance in the bin
# The error bars should be based on different values for different years for a and b, and quantiles with all years put together for c.

# Have to rebin production and density because we have to have the same cutoffs for shade and gap so the numbers are comparable.
totalprodbin_shade_byyear_samebins <- lapply(shadedat[2:3], function(z) logbin_setedges(x = z$light_received/z$crownarea, y = z$production, edges = light_per_area_bins_shadegap))
totalprodbin_gap_byyear_samebins <- lapply(gapdat[2:3], function(z) logbin_setedges(x = z$light_received/z$crownarea, y = z$production, edges = light_per_area_bins_shadegap))

densitybin_shade_byyear_samebins <- lapply(shadedat[2:3], function(z) logbin_setedges(x = z$light_received/z$crownarea, y = NULL, edges = light_per_area_bins_shadegap))
densitybin_gap_byyear_samebins <- lapply(gapdat[2:3], function(z) logbin_setedges(x = z$light_received/z$crownarea, y = NULL, edges = light_per_area_bins_shadegap))

shadegap_stats <- list()

for (i in 1:2) {
  shadegap_stats[[i]] <- data.frame(bin = 1:numbins,
                                    year = c(1990,1995)[i],
                                    shade_gap_production_ratio = totalprodbin_shade_byyear_samebins[[i]]$bin_value/totalprodbin_gap_byyear_samebins[[i]]$bin_value,
                                    shade_gap_density_ratio = densitybin_shade_byyear_samebins[[i]]$bin_value/densitybin_gap_byyear_samebins[[i]]$bin_value)  
}

shadegap_stats <- do.call('rbind', shadegap_stats)

shadegap_stats[shadegap_stats == Inf | is.na(shadegap_stats)] <- NA

shadegap_stats_bin_2census <- shadegap_stats %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(shade_gap_production_ratio),
            density_ratio_mean = mean(shade_gap_density_ratio),
            production_ratio_min = min(shade_gap_production_ratio),
            production_ratio_max = max(shade_gap_production_ratio),
            density_ratio_min = min(shade_gap_density_ratio),
            density_ratio_max = max(shade_gap_density_ratio)) %>%
  cbind(light_per_area_bins_shadegap[,c('bin_midpoint', 'bin_min', 'bin_max')]) %>%
  cbind(mean_n_individuals = light_per_area_bins_shadegap$bin_count / 2)

# Shade tolerance score for each bin (continuous). Fake bin

binshadescore <- function(dat, bindat) {
  dat <- do.call('rbind', dat)
  dat <- subset(dat, !is.na(light_received) & !is.na(pca))
  dat$light_area <- dat$light_received/dat$crownarea
  with(dat, fakebin_across_years(dat_values = pca, dat_classes = light_area, edges = bindat, mean = 'arithmetic', n_census = 2))
}

shadescore_bin_2census <- binshadescore(dat = alltreedat[2:3], bindat = light_per_area_bins_shadegap)

#################################################################################
# Added 2 Nov.
# Same shade gap ratio and average shade score, but do it by diameter this time.
# Figure 5: need to bin (a) the number of shade individuals over the number of gap individuals, (b) total shade production over total gap production, and (c) average shade tolerance in the bin
# The error bars should be based on different values for different years for a and b, and quantiles with all years put together for c.

# Have to rebin production and density because we have to have the same cutoffs for shade and gap so the numbers are comparable.
totalprodbin_shade_byyear_bydiam <- lapply(shadedat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_shadegap))
totalprodbin_gap_byyear_bydiam <- lapply(gapdat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_shadegap))

densitybin_shade_byyear_bydiam <- lapply(shadedat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_shadegap))
densitybin_gap_byyear_bydiam <- lapply(gapdat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_shadegap))

shadegap_stats_bydiam <- list()

for (i in 1:2) {
  shadegap_stats_bydiam[[i]] <- data.frame(bin = 1:numbins,
                                    year = c(1990,1995)[i],
                                    shade_gap_production_ratio = totalprodbin_shade_byyear_bydiam[[i]]$bin_value/totalprodbin_gap_byyear_bydiam[[i]]$bin_value,
                                    shade_gap_density_ratio = densitybin_shade_byyear_bydiam[[i]]$bin_value/densitybin_gap_byyear_bydiam[[i]]$bin_value)  
}

shadegap_stats_bydiam <- do.call('rbind', shadegap_stats_bydiam)

shadegap_stats_bydiam[shadegap_stats_bydiam == Inf | is.na(shadegap_stats_bydiam)] <- NA

shadegap_stats_bin_2census_bydiam <- shadegap_stats_bydiam %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(shade_gap_production_ratio),
            density_ratio_mean = mean(shade_gap_density_ratio),
            production_ratio_min = min(shade_gap_production_ratio),
            production_ratio_max = max(shade_gap_production_ratio),
            density_ratio_min = min(shade_gap_density_ratio),
            density_ratio_max = max(shade_gap_density_ratio)) %>%
  cbind(dbhbin_shadegap[,c('bin_midpoint', 'bin_min', 'bin_max')]) %>%
  cbind(mean_n_individuals = dbhbin_shadegap$bin_count / 2)

# Shade tolerance score for each bin (continuous). Fake bin

binshadescore <- function(dat, bindat) {
  dat <- do.call('rbind', dat)
  dat <- subset(dat, !is.na(light_received) & !is.na(pca))
  with(dat, fakebin_across_years(dat_values = pca, dat_classes = dbh_corr, edges = bindat, mean = 'arithmetic', n_census = 2))
}

shadescore_bin_2census_bydiam <- binshadescore(dat = alltreedat[2:3], bindat = dbhbin_shadegap)



# Export binned data

fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_04oct'
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'shadegap_stats_bin_2census', 'shadescore_bin_2census', 'shadegap_stats_bin_2census_bydiam', 'shadescore_bin_2census_bydiam')

for (i in file_names) {
  write.csv(get(i), file=file.path(fpdata, paste0(i,'.csv')), row.names = FALSE)
}