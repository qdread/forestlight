# Dump forest light data so that we can run stan on the cluster.


# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_04oct'
load(file.path(fpdata, 'rawdataobj.r'))

# Each year and guild separate: indiv production (dbh and production)

prod_dump <- function(dat, fn) {
  require(rstan)
  x <- dat$dbh_corr
  y <- dat$production
  x_pred <- seq(min(x), max(x), length.out = 50)
  xdat <- list(N = length(x), N_pred = length(x_pred), x = x, y = y, x_pred = x_pred)
  with(xdat, stan_rdump(names(xdat), file = fn))
}

# Each year and guild separate: density (dbh only)
# Also add predicted values to this.
# I think this isn't necessary because prod_dump already contains all that info.

dens_dump <- function(dat, fn) {
  require(rstan)
  x <- dat$dbh_corr
  x_pred <- seq(min(x), max(x), length.out = 50)
  xdat <- list(N = length(x), N_pred = length(x_pred), min_x = min(x), x = x, x_pred = x_pred)
  with(xdat, stan_rdump(names(xdat), file = fn))
}

yrs <- c(1990, 1995, 2000, 2005, 2010)

for (i in 2:6) {
  prod_dump(dat = gapdat[[i]], fn = paste0('C:/Users/Q/Dropbox/projects/forestlight/stanrdump/dat_densprod_gap_', yrs[i-1], '.r'))
  prod_dump(dat = shadedat[[i]], fn = paste0('C:/Users/Q/Dropbox/projects/forestlight/stanrdump/dat_densprod_shade_', yrs[i-1], '.r'))
}

# All years for each guild, prod and dens.
allygap <- lapply(gapdat[2:6], function(x) x[,c('dbh_corr','production')])
allygap <- do.call('rbind', allygap)
allyshade <- lapply(shadedat[2:6], function(x) x[,c('dbh_corr','production')])
allyshade <- do.call('rbind', allyshade)

prod_dump(dat = allygap, fn = 'C:/Users/Q/Dropbox/projects/forestlight/stanrdump/dat_densprod_gap_allyrs.r')
prod_dump(dat = allyshade, fn = 'C:/Users/Q/Dropbox/projects/forestlight/stanrdump/dat_densprod_shade_allyrs.r')
