# Full data processing pipeline
# Include new classification from Nadja
# Later, include option to split by growth form or remove shrubs from the dataset.

# Each bin dataframe should include the following
# 1 all individuals
# 2 all individuals classified in a FG
# 3-7 FGs 1 through 5
# 8 individuals not classified in a FG

group_names <- c('all','all_classified','fg1','fg2','fg3','fg4','fg5','unclassified')

load('C:/Users/Q/Dropbox/projects/forestlight/bcidata/bciqcrun.R')

# Modification 25 August: get rid of young and edge trees for all datasets. This will reduce the number of hectares but will be most correct.

for (i in 1:7) {
  assign(paste0('bci.full', i), subset(get(paste0('bci.full', i)), !young & !edge))
  assign(paste0('bci.stem', i), subset(get(paste0('bci.stem', i)), !young & !edge))
}

growth8590 <- subset(growth8590, !young & !edge)
growth9095 <- subset(growth9095, !young & !edge)

library(dplyr)

# Load Nadja's data
fgbci <- read.table('C:/Users/Q/google_drive/ForestLight/data/Ruger/functional_groups_BCI.txt')

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

fgbci_less <- fgbci %>%
  select(sp, grform, fg) %>%
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
save(alltreedat, fgdat, alltree_light_90, alltree_light_95, light_fg_90, light_fg_95, file = 'C:/Users/Q/google_drive/ForestLight/data/rawdataobj_22jan.r')

####################################################################

# Binning and error bars: all years combined ------------------------------

# Script copied and modified on 22 Jan 2018. 
# Procedure: use all 5 FGs (but not unclassified) together to get the bin edges
# Then apply those bin edges to each FG in isolation.
# That way each FG are given the same bin edges.

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

# For now, do not include the figure 5 ratio binning.

# Export binned data

fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census')

for (i in file_names) {
  write.csv(get(i), file=file.path(fpdata, paste0(i,'.csv')), row.names = FALSE)
}

save(list = file_names, file = file.path(fpdata, 'bin_object.RData'))

# Load data and plot
# 04 Oct 2017


# Load data ---------------------------------------------------------------

# Loop through all the csv files and load them into R
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census')

for (i in file_names) {
  assign(i, read.csv(file.path(fpdata, paste0(i,'.csv')), stringsAsFactors = FALSE))
}

# Plot data ---------------------------------------------------------------

library(cowplot)
library(dplyr)

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.

guild_colors <- RColorBrewer::brewer.pal(5, 'Set1')
fg_names <- paste('fg',1:5,sep='')

# Figure 3a
fig_3a <- indivproductionbin_5census %>%
  filter(fg %in% fg_names) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
  scale_color_manual(values = guild_colors) +
  panel_border(colour = 'black')

# Figure 3b
fig_3b <- densitybin_5census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  scale_color_manual(values = guild_colors) +
  panel_border(colour = 'black')

# Figure 3c
fig_3c <- totalproductionbin_5census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Total production (kg ha'^-1,' y'^-1,')'))) +
  scale_color_manual(values = guild_colors) +
  panel_border(colour = 'black')

# Figure 4a
fig_4a <- crownareabin_2census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Total crown area (m'^2, ' ha'^-1,')'))) +
  scale_color_manual(values = guild_colors) +
  panel_border(colour = 'black')

# Figure 4b
fig_4b <- lightreceivedbin_2census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Total light received (W ha'^-1,')'))) +
  scale_color_manual(values = guild_colors) +
  panel_border(colour = 'black')

# Figure 4c
fig_4c <- indivprodperareabin_2census %>%
  filter(fg %in% fg_names) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = expression(paste('Production per unit crown area (kg y'^-1, ' m'^-2,')'))) +
  scale_color_manual(values = guild_colors) +
  panel_border(colour = 'black')

pdf('C:/Users/Q/google_drive/ForestLight/figs/figures_22jan2018/figs3and4.pdf', height = 6, width = 6)
  print(fig_3a + ggtitle('3A'))
  print(fig_3b + ggtitle('3B'))
  print(fig_3c + ggtitle('3C'))
  print(fig_4a + ggtitle('4A'))
  print(fig_4b + ggtitle('4B'))
  print(fig_4c + ggtitle('4C'))
dev.off()

# Make a figure to show the FGs
ggplot(fgbci, aes(x = X1, y = X2, color = factor(fg))) +
  geom_point() +
  labs(x = 'X1 performance/recruitment', y = 'X2 shade/gap') +
  scale_color_manual(values = guild_colors)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/figures_22jan2018/funcgroups.pdf', height=6, width=6)
