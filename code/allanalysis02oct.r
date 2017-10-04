# Modified analysis 02 Oct 2017.

# Do all the function fits using Caroline's MLE code.
# Include 2 clusters rather than 3.
# Don't worry about quantiles.

load('C:/Users/Q/Dropbox/projects/forestlight/bcidata/bciqcrun.R')

# Modification 25 August: get rid of young and edge trees for all datasets. This will reduce the number of hectares but will be most correct.

for (i in 1:7) {
  assign(paste0('bci.full', i), subset(get(paste0('bci.full', i)), !young & !edge))
  assign(paste0('bci.stem', i), subset(get(paste0('bci.stem', i)), !young & !edge))
}

growth8590 <- subset(growth8590, !young & !edge)
growth9095 <- subset(growth9095, !young & !edge)

library(dplyr)

wright <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/wright2010.csv', stringsAsFactors = FALSE)
wright[wright == -99] <- NA # Unknown values were given a value of -99
wright$SPECIES.[109:110] <- c('simplex_var1', 'simplex_var2') # Correct duplicate named subspecies.
wright$Taxon <- with(wright, paste(GENUS., SPECIES.))

wright_df <- with(wright, data.frame(Taxon, mrt = MRT25SAP/100, rgr = RGR95SAP, logit_mrt = qlogis(MRT25SAP/100), log_rgr = log10(RGR95SAP), stringsAsFactors = FALSE))
wright_df <- subset(wright_df, !is.na(mrt))

wright_pca <- with(wright_df, prcomp(data.frame(qlogis(mrt), log10(rgr)), scale=TRUE, center=TRUE)) # 90% of variation on the growth-mortality single axis. Nice.
pca_scores <- wright_pca$x[,1]

set.seed(4930116)
k_bothtrans <- kmeans(x = wright_df[,c('logit_mrt','log_rgr')], centers=2, nstart=25)
wright_df$cluster <- k_bothtrans$cluster
wright_df$pca <- pca_scores

wright_df$tol_wright <- c('G','S')[wright_df$cluster]

# Correct wright_df entries that are not correct.
wright_df$Taxon[grep('Beilsc',wright_df$Taxon)] <- 'Beilschmiedia pendula'
wright_df$Taxon[grep('Cestrum',wright_df$Taxon)] <- 'Cestrum megalophyllum'
wright_df$Taxon[grep('phyllu arg',wright_df$Taxon)] <- 'Chrysophyllum argenteum'
wright_df$Taxon[grep('Coccol',wright_df$Taxon)] <- 'Coccoloba manzinellensis'
wright_df$Taxon[grep('Tabern',wright_df$Taxon)] <- 'Tabernaemontana arborea'
wright_df$Taxon[grep('var1',wright_df$Taxon)] <- 'Swartzia simplex_var.grandiflora'
wright_df$Taxon[grep('var2',wright_df$Taxon)] <- 'Swartzia simplex_var.ochnacea'
wright_df$Taxon[grep('colorado',wright_df$Taxon)] <- 'Eugenia coloradoensis'
wright_df$Taxon[grep('Croton',wright_df$Taxon)] <- 'Croton billbergianus'

bci_lookup <- read.delim('C:/Users/Q/Dropbox/projects/forestlight/bcidata/ViewTax.txt', stringsAsFactors = FALSE)

bci_lookup <- bci_lookup %>%
  mutate(Taxon = paste(Genus, SpeciesName)) 

taxmatch <- bci_lookup$Taxon %in% wright_df$Taxon

bci_lookup_wright <- left_join(bci_lookup, wright_df) %>%
  rename(sp = Mnemonic) %>%
  select(sp, mrt, rgr, pca, tol_wright)

# 2. Convert all census data frames to the correct units (dbh in cm, not mm; agb in kg, not Mg)
# Simultaneously, calculate biomass increments for each stem from the previous census to the current one. Combine everything into a single data frame if possible.
# Also get rid of the young (secondary) forest patches.

bci_production <- list()

for (i in 2:7) {
  dat_i <- get(paste0('bci.full',i))
  dat_iminus1 <- get(paste0('bci.full',i-1))
  bci_production[[i-1]] <- dat_i$agb_corr - dat_iminus1$agb_corr
  dat_i$production <- dat_i$agb_corr - dat_iminus1$agb_corr
  assign(paste0('bci.full', i), dat_i)
}

bci_production <- do.call(cbind, bci_production)
bci_production <- as.data.frame(bci_production)
names(bci_production) <- paste0('production', c('8085','8590','9095','9500','0005','0510'))

bcicensusdat <- list()

for (i in 2:7) {
  bcicensusdat[[i-1]] <- get(paste0('bci.full',i)) %>%
    filter(DFstatus == 'alive') %>%
    mutate(dbh_corr = dbh_corr/10,
           agb_corr = agb_corr * 1000,
           production = production * 1000) %>%
    filter(!young)
}


# 3. Join with shade tolerance (Wright) groups.

for (i in 1:6) {
  bcicensusdat[[i]] <- bcicensusdat[[i]] %>%
    mutate(Taxon = paste(Genus, Species)) %>%
    left_join(wright_df)
}

# 4. Join the appropriate years with Nadja's light data.
# bcicensusdat[[2]] is 1985-1990
# bcicensusdat[[3]] is 1990-1995

bcicensusdat[[2]] <- bcicensusdat[[2]] %>%
  mutate(tag = as.numeric(tag)) %>%
  left_join(growth8590 %>% select(tag, light, young, edge))

bcicensusdat[[3]] <- bcicensusdat[[3]] %>%
  mutate(tag = as.numeric(tag)) %>%
  left_join(growth9095 %>% select(tag, light, young, edge))

# 5. Run allometries to get total light received.

# Function to get a rough approximation of insolation by latitude.

insolation <- function(lat) {
  lat <- lat * pi/180 # to radians
  y <- sin(lat)
  0.25 * 1367 * (1 - 0.482 * (3*y^2 - 1)/2)
}

# Insolation at BCI, 9.2 degrees N
(insol_bci <- insolation(9.2))

# Function to get tree height and crown dimensions from dbh
# Use same parameters for all trees, taken from Bohlman and O'Brien

tp <- function(dbh) {
  h <- exp(.438 + .595 * log(dbh))    # Height
  cd <- exp(-.157 + .702 * log(dbh))  # Crown depth
  cr <- exp(-.438 + .658 * log(dbh))  # Crown radius
  cV <- exp(-.681 + 2.02 * log(dbh))  # Crown volume
  data.frame(h=h, cd=cd, cr=cr, cV=cV)
}

for (i in 2:3) {
  
  crowndim <- tp(bcicensusdat[[i]]$dbh_corr) 
  bcicensusdat[[i]]$crownarea <- pi * crowndim$cr^2
  bcicensusdat[[i]] <- transform(bcicensusdat[[i]], light_received = light * crownarea * insol_bci)
  
}

# Classification of light into 3 groups.
# This is percent of full irradiance, not the total light received.

for (i in 2:3) {
  
  light_groups <- cut(bcicensusdat[[i]]$light, breaks = 3)
  table(light_groups)
  light_groupcodes <- factor(light_groups, 
                             labels = c('Low light','Intermediate light','High light'))
  
  bcicensusdat[[i]]$light_group <- light_groupcodes
  
}

alltreedat <- list()
shadedat <- gapdat <- unclassifieddat <- list()

for (i in 1:6) {
  
  alltreedat[[i]] <- subset(bcicensusdat[[i]], 
                            !is.na(dbh_corr) & production > 0)
  shadedat[[i]] <- subset(alltreedat[[i]], tol_wright == 'S')
  gapdat[[i]] <- subset(alltreedat[[i]], tol_wright == 'G')
  unclassifieddat[[i]] <- subset(alltreedat[[i]], is.na(tol_wright))
}

# for 2 and 3, make additional data frames that have only trees with light measurements.
alltree_light_90 <- subset(alltreedat[[2]], !is.na(light))
alltree_light_95 <- subset(alltreedat[[3]], !is.na(light))
shade_light_90 <- subset(shadedat[[2]], !is.na(light))
shade_light_95 <- subset(shadedat[[3]], !is.na(light))
gap_light_90 <- subset(gapdat[[2]], !is.na(light))
gap_light_95 <- subset(gapdat[[3]], !is.na(light))
unclassified_light_90 <- subset(unclassifieddat[[2]], !is.na(light))
unclassified_light_95 <- subset(unclassifieddat[[3]], !is.na(light))


# Binning and error bars: all years combined ------------------------------

# Concatenate the 5 censuses (3-7) excluding 1985, and find log bin edges.
# Keep the size classes constant across density and production.

# All years together are used to create the 20 bin boundaries, then the binning is done separately for each year.
# Use the bounds from the dbh bin to create the bins for all other scalings.

source('code/allfunctions27july.r')
numbins <- 20

allyeardbh_alltree <- unlist(lapply(alltreedat[2:6], '[', , 'dbh_corr'))
allyeardbh_shade <- unlist(lapply(shadedat[2:6], '[', , 'dbh_corr'))
allyeardbh_gap <- unlist(lapply(gapdat[2:6], '[', , 'dbh_corr'))
allyeardbh_unclassified <- unlist(lapply(unclassifieddat[2:6], '[', , 'dbh_corr'))

dbhbin_alltree <- logbin(x = allyeardbh_alltree, y = NULL, n = numbins)
dbhbin_shade <- logbin(x = allyeardbh_shade, y = NULL, n = numbins)
dbhbin_gap <- logbin(x = allyeardbh_gap, y = NULL, n = numbins)
dbhbin_unclassified <- logbin(x = allyeardbh_unclassified, y = NULL, n = numbins)

# The dbhbin_ objects can be used as the edges argument in logbin_setedges()

# Bin density and individual production by year using the above generated bin edges.
# dbh
dbhbin_alltree_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_alltree))
dbhbin_shade_byyear <- lapply(shadedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_shade))
dbhbin_gap_byyear <- lapply(gapdat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_gap))
dbhbin_unclassified_byyear <- lapply(unclassifieddat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_unclassified))

# Find median, minimum and maximum bin values across years.

bin_across_years <- function(binlist) {
  binvals <- do.call('cbind', lapply(binlist, '[', , 'bin_value'))
  data.frame(bin_midpoint = binlist[[1]]$bin_midpoint,
             bin_min = binlist[[1]]$bin_min,
             bin_max = binlist[[1]]$bin_max,
             bin_yvalue = apply(binvals, 1, median),
             bin_ymin = apply(binvals, 1, min),
             bin_ymax = apply(binvals, 1, max))
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

# Individual production (stupid)
# Take the mean and 2.5, 25, 50, 75, 97.5 quantiles within each bin.
# Do it across all years.

allyearprod_alltree <- unlist(lapply(alltreedat[2:6], '[', , 'production'))
allyearprod_shade <- unlist(lapply(shadedat[2:6], '[', , 'production'))
allyearprod_gap <- unlist(lapply(gapdat[2:6], '[', , 'production'))
allyearprod_unclassified <- unlist(lapply(unclassifieddat[2:6], '[', , 'production'))

fakebin_across_years <- function(dat_values, dat_classes, edges) {
  qprobs <- c(0.025, 0.5, 0.975)
  # add some padding just in case
  mins <- edges$bin_min
  mins[1] <- 0
  maxes <- edges$bin_max
  maxes[length(maxes)] <- Inf
  
  binstats <- t(sapply(1:length(mins), function(x) {
    indivs <- dat_values[dat_classes >= mins[x] & dat_classes < maxes[x]]
    c(mean = mean(indivs), quantile(indivs, probs = qprobs))
  }))
  dimnames(binstats)[[2]] <- c('mean', 'q025', 'median', 'q975')
  data.frame(bin_midpoint = edges$bin_midpoint,
             bin_min = edges$bin_min,
             bin_max = edges$bin_max,
             binstats)
}

prodbin_alltree_5census <- fakebin_across_years(dat_values = allyearprod_alltree, dat_classes = allyeardbh_alltree, edges = dbhbin_alltree)
prodbin_shade_5census <- fakebin_across_years(dat_values = allyearprod_shade, dat_classes = allyeardbh_shade, edges = dbhbin_shade)
prodbin_gap_5census <- fakebin_across_years(dat_values = allyearprod_gap, dat_classes = allyeardbh_gap, edges = dbhbin_gap)
prodbin_unclassified_5census <- fakebin_across_years(dat_values = allyearprod_unclassified, dat_classes = allyeardbh_unclassified, edges = dbhbin_unclassified)

indivproductionbin_5census <- rbind(cbind(guild = 'all', prodbin_alltree_5census),
                                    cbind(guild = 'shade', prodbin_shade_5census),
                                    cbind(guild = 'gap', prodbin_gap_5census),
                                    cbind(guild = 'unclassified', prodbin_unclassified_5census))

# Total production
# Do the binning for each year separately, as for density, then find min, max, and median.
totalprodbin_alltree_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_alltree))
totalprodbin_shade_byyear <- lapply(shadedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_shade))
totalprodbin_gap_byyear <- lapply(gapdat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_gap))
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
crownareabin_shade_byyear <- lapply(shadedat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownarea, edges = dbhbin_shade))
crownareabin_gap_byyear <- lapply(gapdat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownarea, edges = dbhbin_gap))
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
lightreceivedbin_shade_byyear <- lapply(shadedat[2:3], function(z) with(subset(z, !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_alltree)))
lightreceivedbin_gap_byyear <- lapply(gapdat[2:3], function(z) with(subset(z, !is.na(light_received)), logbin_setedges(x = dbh_corr, y = light_received, edges = dbhbin_alltree)))
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
  with(dat, fakebin_across_years(dat_values = prod_area, dat_classes = light_area, edges = bindat))
}

# Bin the entire light received per crown area dataset for 1990 and 1995 into a single set of bin edges.
light_per_area_all <- unlist(lapply(alltreedat[2:3], '[', , 'light_received'))/unlist(lapply(alltreedat[2:3], '[', , 'crownarea'))
light_per_area_shade <- unlist(lapply(shadedat[2:3], '[', , 'light_received'))/unlist(lapply(shadedat[2:3], '[', , 'crownarea'))
light_per_area_gap <- unlist(lapply(gapdat[2:3], '[', , 'light_received'))/unlist(lapply(gapdat[2:3], '[', , 'crownarea'))
light_per_area_unclassified <- unlist(lapply(unclassifieddat[2:3], '[', , 'light_received'))/unlist(lapply(unclassifieddat[2:3], '[', , 'crownarea'))

light_per_area_bins_all <- logbin(x = na.omit(light_per_area_all), n = numbins)
light_per_area_bins_shade <- logbin(x = na.omit(light_per_area_shade), n = numbins)
light_per_area_bins_gap <- logbin(x = na.omit(light_per_area_gap), n = numbins)
light_per_area_bins_unclassified <- logbin(x = na.omit(light_per_area_unclassified), n = numbins)

indivprodperareabin_alltree_2census <- binprod(dat = alltreedat[2:3], bindat = light_per_area_bins_all, xvar = 'production')
indivprodperareabin_shade_2census <- binprod(dat = shadedat[2:3], bindat = light_per_area_bins_shade, xvar = 'production')
indivprodperareabin_gap_2census <- binprod(dat = gapdat[2:3], bindat = light_per_area_bins_gap, xvar = 'production')
indivprodperareabin_unclassified_2census <- binprod(dat = unclassifieddat[2:3], bindat = light_per_area_bins_unclassified, xvar = 'production')

indivprodperareabin_2census <- rbind(cbind(guild = 'all', indivprodperareabin_alltree_2census),
                              cbind(guild = 'shade', indivprodperareabin_shade_2census),
                              cbind(guild = 'gap', indivprodperareabin_gap_2census),
                              cbind(guild = 'unclassified', indivprodperareabin_unclassified_2census))


# Figure 5: need to bin (a) the number of shade individuals over the number of gap individuals, (b) total shade production over total gap production, and (c) average shade tolerance in the bin
# The error bars should be based on different values for different years for a and b, and quantiles with all years put together for c.

# Have to rebin production and density because we have to have the same cutoffs for shade and gap so the numbers are comparable.
totalprodbin_shade_byyear_samebins <- lapply(shadedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = light_per_area_bins_all))
totalprodbin_gap_byyear_samebins <- lapply(gapdat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = light_per_area_bins_all))

densitybin_shade_byyear_samebins <- lapply(shadedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = light_per_area_bins_all))
densitybin_gap_byyear_samebins <- lapply(gapdat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = light_per_area_bins_all))

shadegap_stats <- list()

for (i in 1:5) {
  shadegap_stats[[i]] <- data.frame(bin = 1:numbins,
                                    year = c(1990,1995,2000,2005,2010)[i],
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
  cbind(light_per_area_bins_all[,c(1,4,5)])
  
# Shade tolerance score for each bin (continuous). Fake bin

binshadescore <- function(dat, bindat) {
  dat <- do.call('rbind', dat)
  dat <- subset(dat, !is.na(light_received) & !is.na(pca))
  dat$light_area <- dat$light_received/dat$crownarea
  with(dat, fakebin_across_years(dat_values = pca, dat_classes = light_area, edges = bindat))
}

shadescore_bin_2census <- binshadescore(dat = alltreedat[2:3], bindat = light_per_area_bins_all)

# Export binned data

fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_04oct'
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'shadegap_stats_bin_2census', 'shadescore_bin_2census')

for (i in file_names) {
  write.csv(get(i), file=file.path(fpdata, paste0(i,'.csv')), row.names = FALSE)
}

# Scalings ----------------------------------------------------------------

# Required scalings:
# 1. Fig 3a. Individual production scaling.
# 2. Fig 3b. Density scaling.
# 3. Fig 3c. Total production scaling. (binned)
# 4. Fig 4a. Total crown area scaling (binned)
# 5. Fig 4b. Total light received scaling (binned)
# 6. Fig 4c. Production by light received
# 7. Fig 5. Shade to gap ratios, Shade to gap production ratios, Shade tolerance score by light received regression


# Individual production scaling -------------------------------------------

# Straight line fit and curved line fit.

# Problem: Is there even a point to fitting this? It looks like a modeled thing.

# Plot them.
library(cowplot)
ggplot(shadedat[[6]], aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point()
ggplot(shadedat[[6]], aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_hex()
ggplot(gapdat[[6]], aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point() +
  stat_smooth(method='lm', color = 'blue') + stat_smooth(method='auto', color = 'red')

nls_fit_gap <- nls(logprod ~ C * (exp(k * logdbh)), 
    data = transform(gapdat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)), 
    start = c(C=100, k=-2), 
    control = list(maxiter=500))

nls_fit_shade <- nls(logprod ~ C * (exp(k * logdbh)), 
    data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)), 
    start = c(C=100, k=-2), 
    control = list(maxiter=500))

linear_fit_gap <- lm(logprod ~ logdbh, 
                      data = transform(gapdat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))

library(mgcv)
gam_fit_gap <- gam(logprod ~ logdbh, 
                   data = transform(gapdat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)), k=3)

plot(gam_fit_gap, pages=1, residuals=T, all.terms=T)
summary(gam_fit_gap)

gam_fit_shade <- gam(logprod ~ s(logdbh, bs='cr', k = 3), 
                   data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))
linear_fit_shade <- lm(logprod ~ logdbh, 
                     data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))

predict(gam_fit_shade, newdata=)

ggplot(shadedat[[6]], aes(x=dbh_corr, y=production)) + 
  geom_point() +
  stat_function(geom='line', fun = function(x,beta) beta[1] + beta[2]*log10(x), args = list(beta=linear_fit_shade$coefficients), color = 'red', size=2) +
  geom_line(data=data.frame(cbind(shadedat[[6]], ypred=predict(gam_fit_shade))), aes(y=ypred), color='skyblue', size=2) +
  scale_x_log10() + scale_y_log10() 

ggplot(shadedat[[6]], aes(x=dbh_corr, y=production)) + 
  geom_point() +
  stat_smooth(method = 'lm', se = FALSE, color = 'red', size = 2) +
  stat_smooth(method = 'lm', formula = y ~ splines::bs(x,3), se = FALSE, color = 'goldenrod', size = 2) +
  stat_smooth(method = 'lm', formula = y ~ splines::bs(x,knots=c(1,2.7)), se = FALSE, color = 'green', size = 2) +
  stat_smooth(method = 'lm', formula = y ~ x + I(x^0.5), se = FALSE, color = 'blue', size = 2) +
  scale_x_log10() + scale_y_log10() 

spline_fit_shade <- lm(logprod ~ splines::bs(logdbh,3), 
                       data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))  
spline_fit2_shade <- lm(logprod ~ splines::bs(logdbh, knots = c(1, 2.7)), 
                       data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))  
sqrt_fit_shade <- lm(logprod ~ logdbh + I(logdbh^0.5),
                        data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))

AIC(sqrt_fit_shade)
AIC(spline_fit_shade)
AIC(spline_fit2_shade)
AIC(linear_fit_shade)

library(lme4)
mixed_sqrt_shade <- lmer(logprod ~ logdbh + I(logdbh^0.5) + (1|sp),
     data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))
mixed_lm_shade <- lmer(logprod ~ logdbh + (1|sp),
                         data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))
mixed_spline_shade <- lmer(logprod ~ splines::bs(logdbh, knots = c(1, 2.7)) + (1|sp),
                       data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))

AIC(mixed_spline_shade)
AIC(mixed_lm_shade)
AIC(mixed_sqrt_shade)
# Plot the fits hopefully with predictive intervals

nls_pars_gap <- nls_fit_gap$m$getPars()

ggplot(gapdat[[6]], aes(x=dbh_corr, y=production)) + 
  stat_function(geom='line', fun = function(x,beta) beta[1] + beta[2]*log10(x), args = list(beta=linear_fit_gap$coefficients), color = 'red', size=2) +
  stat_function(geom='line', fun = function(x,pars) (pars[1] * exp(pars[2] * log10(x))), args = list(pars=nls_pars_gap), color = 'skyblue', size=2) +
  scale_x_log10() + scale_y_log10() + geom_point()
  
  



# Test code below ---------------------------------------------------------



# Check whether there is some artifact in the tree production scalings
ggplot(shadedat[[6]] %>% filter(sp=='hybapr'), aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point()
ggplot(shadedat[[6]] %>% filter(sp=='caseac'), aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point()
ggplot(shadedat[[6]] %>% filter(sp=='tachve'), aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point()
ggplot(shadedat[[6]] %>% filter(sp %in% unique(shadedat[[6]]$sp[1:50])), aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point() + facet_wrap(~sp)

# These appear to be modeled. Fitting curves is maybe a bad idea?

# get basal area increments.
bci.full7$basalareaincrement <- pi*(bci.full7$dbh_corr/2)^2 - pi*(bci.full6$dbh_corr/2)^2
bci.full7$diameterincrement <- bci.full7$dbh_corr - bci.full6$dbh_corr

bci.full7 %>%
  filter(sp %in% unique(bci.full7$sp[1:100])) %>%
  ggplot(aes(x=dbh_corr, y=basalareaincrement)) +
  scale_x_log10() + scale_y_log10() +
  geom_point() +
  facet_wrap(~ sp)

bci.full7 %>%
  filter(sp %in% unique(bci.full7$sp[1:50])) %>%
  ggplot(aes(x=dbh_corr, y=diameterincrement)) +
  geom_point() +
  facet_wrap(~ sp)
