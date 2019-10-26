# ALL DATA DUMPS FOR (FINAL) FORESTLIGHT MODEL FITTING
# This script creates lists and writes them to Rdump text files that can be read into CmdStan for model fitting.
# QDR / ForestLight / 25 Oct 2019

# Needed data dumps:
# DBH and biomass production for density-biomass growth scalings
# DBH and diameter growth rate for density-diameter growth scalings
# Light per crown area and growth per crown area for growth/area-light/area scalings

# --------- #
# Load data #
# --------- #

# file paths for input and output
fpdata <- '~/google_drive/ForestLight/data'
fpdump <- '~/Dropbox/projects/forestlight/stanrdump'

# load raw data
load(file.path(fpdata, 'rawdataobj_alternativecluster.r'))

library(tidyverse)
library(rstan)

# -------------------------------- #
# Define function to create rdumps #
# -------------------------------- #

create_rdump <- function(dat, xvar, yvar, file_name = NULL, subsample = NULL, random_seed = NULL) {
  require(rstan)
  
  # Subset of data if specified
  if (!is.null(subsample) && nrow(dat) > subsample) {
	if !is.null(random_seed) set.seed(random_seed)
    dat <- dat[sample(nrow(dat), subsample, replace = FALSE), ]
  }
  
  # Minimum and maximum dbh of BCI trees
  ll <- 1
  ul <- 316

  x <- dat[, xvar]
  y <- dat[, yvar]
  xdat <- list(N = length(x), x = x, y = y, x_min = min(x), x_max = max(x), LL = ll, UL = ul)
  
  if (!is.null(file_name)) {
    with(xdat, stan_rdump(names(xdat), file = file_name))
  } else {
    return(xdat)
  }
}

# -------------------------------- #
# Create rdumps and write to files #
# -------------------------------- #

# Do 1995 in all cases. For all cases, do all trees together and one for each FG
fgs <- c('fg1','fg2','fg3','fg4','fg5','unclassified')

### biomass growth ~ diameter scalings

create_rdump(alltreedat[[3]], 'dbh_corr', 'production', file_name = file.path(fpdump, 'dump_alltree_1995.r'))

imap(fgs, ~ create_rdump(fgdat[[.y]][[3]], 'dbh_corr' , 'production', file_name = file.path(fpdump, paste0('dump_', .x, '_1995.r'))))

### diameter growth ~ diameter scalings
create_rdump(alltreedat[[3]], 'dbh_corr', 'diam_growth_rate', file_name = file.path(fpdump, 'dump_diamgrowthscaling_alltree_1995.r'))

imap(fgs, ~ create_rdump(fgdat[[.y]][[3]], 'dbh_corr' , 'diam_growth_rate', file_name = file.path(fpdump, paste0('dump_diamgrowthscaling_', .x, '_1995.r'))))

### light per area ~ growth per area scalings

dat95 <- alltree_light_95 %>%
  select(dbh_corr, production, light_received, crownarea, crownvolume, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))
  
create_rdump(dat95, 'light_area', 'production_area', file_name = file.path(fpdump, 'dump_light_alltree_1995.r'))

dat95 %>%
	group_by(fg) %>%
	group_walk(~ create_rdump(., 'light_area', 'production_area', file_name = file.path(fpdump, paste0('dump_light_', .y, '_1995.r'))))
	
### individual incoming light ~ diameter scalings
create_rdump(dat95, 'light_received', 'dbh_corr', file_name = file.path(fpdump, 'dump_rawlightscaling_alltree_1995.r'))

dat95 %>%
	group_by(fg) %>%
	group_walk(~ create_rdump(., 'light_received', 'dbh_corr', file_name = file.path(fpdump, paste0('dump_rawlightscaling_', .y, '_1995.r'))))

### crown volume ~ diameter scalings
create_rdump(dat95, 'crownvolume', 'dbh_corr', file_name = file.path(fpdump, 'dump_volumescaling_alltree_1995.r'))

dat95 %>%
	group_by(fg) %>%
	group_walk(~ create_rdump(., 'crownvolume', 'dbh_corr', file_name = file.path(fpdump, paste0('dump_volumescaling_', .y, '_1995.r'))))

# ------------------------------------------------------------------------------ #
# Create CSVs of minima, maxima, number of individuals, and normalization totals #
# ------------------------------------------------------------------------------ #

# This info is used to make some corrections in the output after the models have been run, for plotting purposes.

years <- c(1985, 1990, 1995, 2000, 2005, 2010)

##### EDIT BELOW DIS SPOT ######

# Minima, maxima, and number of individuals 
valall <- map2_dfr(alltreedat, years,
				   ~ data.frame(fg = 'alltree', year = .y, xmin = min(.x$dbh_corr), n = nrow(.x)))

valfg <- imap_dfr(fgdat, function(dat, fg_name) {
	map2_dfr(dat, years,
			 ~ data.frame(fg = fg_name, year = .y, xmin = min(dat$dbh_corr), n = nrow(dat)))
})

min_n <- rbind(valall, valfg)
write_csv(min_n, file.path(fpdump, 'min_n.csv'))

# Minima, maxima, and number of individuals for 1995 only for trees with light measured
valall <- data.frame(fg = 'alltree', year = 1995, xmin = with(alltree_light_95, min(light_received/crownarea)), n = nrow(alltree_light_95))
valfg <- alltree_light_95 %>%
  group_by(fg) %>%
  summarize(xmin = min(light_received/crownarea),
            n = n())

min_n <- data.frame(fg = c('alltree', 'fg1','fg2','fg3','fg4','fg5','unclassified'),
                    year = 1995,
                    xmin = c(valall$xmin, valfg$xmin),
                    n = c(valall$n, valfg$n))

write_csv(min_n, file.path(fpdump, 'min_n_lighttrees.csv'))

# Production totals

alltreeprod <- map2_dfr(alltreedat, years, ~ data.frame(year = .y, fg = 'alltree', production = sum(.x$production)))
fgprod <- map2_dfr(alltreedat, years, function(x, y) {
  x <- x %>%
    mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg',fg))) %>%
    group_by(fg) %>%
    summarize(production = sum(production))
  data.frame(year = y, x)	
})

prodtotals <- rbind(alltreeprod, fgprod) %>%
  arrange(year, fg)

write_csv(prodtotals, file = file.path(fpdump, 'production_total.csv'))

# Total production for all trees that have light values
alltreeprod <- sum(alltree_light_95$production)
fgprod <- alltree_light_95 %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg',fg))) %>%
  group_by(fg) %>%
  summarize(production = sum(production))

prodtotals <- data.frame(year = 1995, 
                         rbind(c(fg = 'alltree', production = alltreeprod), fgprod))

write_csv(prodtotals, file = file.path(fpdump, 'production_total_lighttrees.csv'))

# Light received totals
alltreelightrec <- sum(dat95$light_received)
fglightrec <- dat95 %>%
  group_by(fg) %>%
  summarize(light_received = sum(light_received))

lightrectotals <- data.frame(year = 1995, 
                             rbind(c(fg = 'alltree', light_received = alltreelightrec), fglightrec))

write_csv(lightrectotals, file = file.path(fpdump, 'lightrec_total.csv'))

# Crown volume totals
alltreecrownvol <- sum(dat95$crownvolume)
fgcrownvol <- dat95 %>%
  group_by(fg) %>%
  summarize(crownvolume = sum(crownvolume))

crownvolumetotals <- data.frame(year = 1995, 
                             rbind(c(fg = 'alltree', crownvolume = alltreecrownvol), fgcrownvol))

write_csv(crownvolumetotals, file = file.path(fpdump, 'crownvol_total.csv'))
