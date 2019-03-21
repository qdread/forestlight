# Full data processing pipeline
# June 2018 version: Use the fg5 classification.

# Last modified 25 June 2018 to use Nadja's newest FG classification

# Each bin dataframe should include the following
# 1 all individuals
# 2 all individuals classified in a FG
# 3-7 FGs 1 through 5
# 8 individuals not classified in a FG

group_names <- c('all','all_classified','fg1','fg2','fg3','fg4','fg5','unclassified')

load('~/google_drive/ForestLight/data/BCI_raw/bcidata/bciqcrun.R')

# Modification 25 August: get rid of young and edge trees for all datasets. This will reduce the number of hectares but will be most correct.

for (i in 1:7) {
  assign(paste0('bci.full', i), subset(get(paste0('bci.full', i)), !young & !edge))
  assign(paste0('bci.stem', i), subset(get(paste0('bci.stem', i)), !young & !edge))
}

growth8590 <- subset(growth8590, !young & !edge)
growth9095 <- subset(growth9095, !young & !edge)

library(dplyr)

# Load Nadja's data (new functional groups 25 June)
# fg5 is the new column (we originally used fg from the older df)
fgbci <- read.table('~/google_drive/ForestLight/data/Ruger/fgroups_dynamics_new.txt', stringsAsFactors = FALSE)

# Correct functional groups so that: 1 fast, 2 pioneer, 3 slow, 4 breeder, 5 intermediate
# Old 1,2,3,4,5 --> New 2,3,1,4,5
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5))

# Currently X1new is high for slow species and low for fast species
# Currently X2new is high for pioneer species and low for breeder species
# Correct these
fgbci$PC_slow_to_fast <- -fgbci$X1new
fgbci$PC_breeder_to_pioneer <- fgbci$X2new

# 2. Convert all census data frames to the correct units (dbh in cm, not mm; agb in kg, not Mg)
# Simultaneously, calculate biomass increments for each stem from the previous census to the current one. Combine everything into a single data frame if possible.
# Use harmonic mean to annualize the 5-year biomass increment into a 1-year biomass increment using the small end of the interval
# Also get rid of the young (secondary) forest patches.

### Function for annualizing biomass increment
# census interval should be 5 years and desired new interval should be 1 year
annual_increment <- function(agb_old, agb_new, census_interval = 5, new_interval = 1){
  rate <-  (agb_new / agb_old)^(1/census_interval) - 1
  agb_oldplus1year <- agb_old * (1 + rate)^new_interval
  return(agb_oldplus1year - agb_old)
}

bci_production <- list()

for (i in 2:7) {
  dat_i <- get(paste0('bci.full',i))
  dat_iminus1 <- get(paste0('bci.full',i-1))
  bci_production[[i-1]] <- annual_increment(agb_old = dat_iminus1$agb_corr, agb_new = dat_i$agb_corr)
  dat_i$production <- bci_production[[i-1]]
  assign(paste0('bci.full', i), dat_i)
}

bci_production <- do.call(cbind, bci_production)
bci_production <- as.data.frame(bci_production)
names(bci_production) <- paste0('production', c('8085','8590','9095','9500','0005','0510'))

# Amendment 11 Oct. : Add dbh increment, used to detect production outliers.

for (i in 2:7) {
  dat_i <- get(paste0('bci.full',i))
  dat_iminus1 <- get(paste0('bci.full',i-1))
  dbhlastcensus <- dat_iminus1$dbh_corr
  dbhlastcensus[is.na(dbhlastcensus)] <- 0
  dat_i$dbhinc <- dat_i$dbh_corr - dbhlastcensus
  dat_i$dbhlastcensus <- dbhlastcensus
  assign(paste0('bci.full', i), dat_i)
}




bcicensusdat <- list()

for (i in 2:7) {
  bcicensusdat[[i-1]] <- get(paste0('bci.full',i)) %>%
    filter(DFstatus == 'alive') %>%
    mutate(dbh_corr = dbh_corr/10,
           agb_corr = agb_corr * 1000,
           production = production * 1000,
           dbhinc = dbhinc/10,
           dbhlastcensus = dbhlastcensus/10) %>%
    filter(!young)
}

####
# Amendment 11 Oct.
# Remove production outliers.
# Criteria: 
# 1. If the tree was not in the previous census and its current dbh is >10, remove.
# 2. If the tree grew over 20 cm dbh between censuses, remove.

bcicensusdat <- lapply(bcicensusdat, function(x) {
  x$production[x$dbhinc > 20 | (x$dbhinc > 10 & x$dbhlastcensus == 0)] <- NA
  x
})


# 3. Join with shade tolerance (Rüger) groups.
# Here is where the fg5 is assigned as the fg classification we are using.

fgbci_less <- fgbci %>%
  select(sp, grform, fg5, PC_slow_to_fast, PC_breeder_to_pioneer) %>%
  rename(fg = fg5, X1 = PC_slow_to_fast, X2 = PC_breeder_to_pioneer) %>%
  mutate(sp = tolower(sp))

for (i in 1:6) {
  bcicensusdat[[i]] <- bcicensusdat[[i]] %>%
    left_join(fgbci_less)
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

######
# If it's needed to edit the allometries for different species groups, add it here.
######

for (i in 2:3) {
  
  crowndim <- tp(bcicensusdat[[i]]$dbh_corr) 
  bcicensusdat[[i]]$crownarea <- pi * crowndim$cr^2
  bcicensusdat[[i]]$crownvolume <- crowndim$cV
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
fgdat <- replicate(6, list())

for (i in 1:6) {
  
  alltreedat[[i]] <- subset(bcicensusdat[[i]], 
                            !is.na(dbh_corr) & production > 0)
}

for (i in 1:6) {
  for (j in 1:5) {
    fgdat[[j]][[i]] <- subset(alltreedat[[i]], fg == j)
  }
  fgdat[[6]][[i]] <- subset(alltreedat[[i]], is.na(fg))
}

# for 2 and 3, make additional data frames that have only trees with light measurements.
alltree_light_90 <- subset(alltreedat[[2]], !is.na(light))
alltree_light_95 <- subset(alltreedat[[3]], !is.na(light))

light_fg_90 <- lapply(fgdat, function(x) subset(x[[2]], !is.na(light)))
light_fg_95 <- lapply(fgdat, function(x) subset(x[[3]], !is.na(light)))


# Save raw data as an object
save(alltreedat, fgdat, alltree_light_90, alltree_light_95, light_fg_90, light_fg_95, file = '~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')

####################################################################

# Binning and error bars: all years combined ------------------------------

# Script copied and modified on 22 Jan 2018. 
# Crown volume added on 21 March 2019.
# Procedure: use all 5 FGs (but not unclassified) together to get the bin edges
# Then apply those bin edges to each FG in isolation.
# That way each FG are given the same bin edges.

library(dplyr)

load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')
group_names <- c('all','all_classified','fg1','fg2','fg3','fg4','fg5','unclassified')
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

## crown volume
crownvolumebin_alltree_byyear <- lapply(alltreedat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownvolume, edges = dbhbin_all))
crownvolumebin_allclassified_byyear <- lapply(alltreedat_classified[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownvolume, edges = dbhbin_allclassified))

crownvolumebin_fg_byyear <- list()

for (i in 1:6) {
  crownvolumebin_fg_byyear[[i]] <- lapply(fgdat[[i]][2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$crownvolume, edges = dbhbin_allclassified))
}

crownvolumebin_all_2census <- bin_across_years(crownvolumebin_alltree_byyear)
crownvolumebin_allclassified_2census <- bin_across_years(crownvolumebin_allclassified_byyear)
crownvolumebin_fg_2census <- lapply(crownvolumebin_fg_byyear, bin_across_years)

crownvolumebin_2census <- cbind(fg = rep(group_names, each = numbins), 
                              rbind(crownvolumebin_all_2census, crownvolumebin_allclassified_2census, do.call('rbind', crownvolumebin_fg_2census)))

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

## light received per unit crown area
lightperareabin_alltree_byyear <- lapply(alltreedat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$light_received/z$crownarea, edges = dbhbin_all))
lightperareabin_allclassified_byyear <- lapply(alltreedat_classified[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$light_received/z$crownarea, edges = dbhbin_allclassified))

lightperareabin_fg_byyear <- list()

for (i in 1:6) {
  lightperareabin_fg_byyear[[i]] <- lapply(fgdat[[i]][2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$light_received/z$crownarea, edges = dbhbin_allclassified))
}

## light received per unit crown volume
lightpervolumebin_alltree_byyear <- lapply(alltreedat[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$light_received/z$crownvolume, edges = dbhbin_all))
lightpervolumebin_allclassified_byyear <- lapply(alltreedat_classified[2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$light_received/z$crownvolume, edges = dbhbin_allclassified))

lightpervolumebin_fg_byyear <- list()

for (i in 1:6) {
  lightpervolumebin_fg_byyear[[i]] <- lapply(fgdat[[i]][2:3], function(z) logbin_setedges(x = z$dbh_corr, y = z$light_received/z$crownvolume, edges = dbhbin_allclassified))
}

## Combine the individual 1995 bins to data frames and then write them to R object.
library(purrr)
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

save(crownareabins1995, crownvolumebins1995, lightreceivedbins1995, lightperareabins1995, lightpervolumebins1995,
     file = '~/google_drive/ForestLight/data/area_and_volume_bins_1995.RData')

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

fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_june2018_alternativecluster'
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'breeder_stats_bydiam_2census', 'breederscore_bin_bydiam_2census', 'breeder_stats_bylight_2census', 'breederscore_bin_bylight_2census', 'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census', 'fastslow_stats_bylight_2census', 'fastslowscore_bin_bylight_2census')

for (i in file_names) {
  write.csv(get(i), file=file.path(fpdata, paste0(i,'.csv')), row.names = FALSE)
}

save(list = file_names, file = file.path(fpdata, 'bin_object.RData'))

# Load data ---------------------------------------------------------------

# Loop through all the csv files and load them into R
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_june2018_alternativecluster'
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'breeder_stats_bydiam_2census', 'breederscore_bin_bydiam_2census', 'breeder_stats_bylight_2census', 'breederscore_bin_bylight_2census', 'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census', 'fastslow_stats_bylight_2census', 'fastslowscore_bin_bylight_2census')

for (i in file_names) {
  assign(i, read.csv(file.path(fpdata, paste0(i,'.csv')), stringsAsFactors = FALSE))
}

# Plot data ---------------------------------------------------------------

library(cowplot)
library(dplyr)

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.
global_diam_xlimits <- c(1, 316) # Maximum and minimum for x axis if it's diameter.

guild_colors <- RColorBrewer::brewer.pal(5, 'Set1')
fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('fast','long-lived pioneer', 'slow', 'short-lived breeder', 'intermediate')

# Figure 3a
fig_3a <- indivproductionbin_5census %>%
  filter(fg %in% fg_names, !is.na(mean), mean > 0) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group') +
  panel_border(colour = 'black')

# Figure 3b
fig_3b <- densitybin_5census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 3c
fig_3c <- totalproductionbin_5census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total production (kg ha'^-1,' y'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 4a
fig_4a <- crownareabin_2census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total crown area (m'^2, ' ha'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 4b
fig_4b <- lightreceivedbin_2census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total light received (W ha'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 4c
fig_4c <- indivprodperareabin_2census %>%
  filter(fg %in% fg_names, !is.na(mean), mean > 0) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = expression(paste('Production per unit crown area (kg y'^-1, ' m'^-2,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group') +
  panel_border(colour = 'black')

pdf('C:/Users/Q/google_drive/ForestLight/figs/figures_june2018_newcluster/figs3and4.pdf', height = 6, width = 7.5)
  print(fig_3a + ggtitle('3A'))
  print(fig_3b + ggtitle('3B'))
  print(fig_3c + ggtitle('3C'))
  print(fig_4a + ggtitle('4A'))
  print(fig_4b + ggtitle('4B'))
  print(fig_4c + ggtitle('4C'))
dev.off()

# Plot functional groups
ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, color = factor(fg5))) +
  geom_point() +
  labs(x = 'X1 slow to fast', y = 'X2 breeders to pioneers') +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group')
ggsave('C:/Users/Q/google_drive/ForestLight/figs/figures_june2018_newcluster/fg5plot (newest groups).pdf')
###
# Added 23 Jan: Ratio figures (fig 5)
# Note that the continuous figures, as they don't use groupings, use all the data.

# Breeder to pioneer by light
fig_blightprod <- breeder_stats_bylight_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Breeder to pioneer production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_blightdens <- breeder_stats_bylight_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
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
  scale_y_continuous(name = 'Low = short-lived breeder, high = long-lived pioneer') +
  panel_border(colour = 'black')

# Breeder to pioneer by diameter
fig_bdiamprod <- breeder_stats_bydiam_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = 'Breeder to pioneer production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_bdiamdens <- breeder_stats_bydiam_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = 'Breeder to pioneer density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_bdiamscore <- breederscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_continuous(name = 'Low = short-lived breeder, high = long-lived pioneer') +
  panel_border(colour = 'black')


####

# Fast to slow by light
fig_flightprod <- fastslow_stats_bylight_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Fast to slow production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_flightdens <- fastslow_stats_bylight_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
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
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = 'Fast to slow production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_fdiamdens <- fastslow_stats_bydiam_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = 'Fast to slow density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_fdiamscore <- fastslowscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_continuous(name = 'Low = slow, high = fast') +
  panel_border(colour = 'black')

pdf('C:/Users/Q/google_drive/ForestLight/figs/figures_june2018_newcluster/ratio_figures.pdf', height = 6, width = 6)
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