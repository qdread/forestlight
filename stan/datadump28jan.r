# Create Rdump for new grouped data to run with CMDSTAN.

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/'
load(file.path(fpdata, 'rawdataobj_22jan.r'))

prod_dump <- function(dat, to_file = FALSE, fn = NULL) {
  require(rstan)
  x <- dat$dbh_corr
  y <- dat$production
  xdat <- list(N = length(x), x = x, y = y, x_min = min(x))
  if (to_file) {
    with(xdat, stan_rdump(names(xdat), file = fn))
  } else {
    return(xdat)
  }
}

fpdump <- 'C:/Users/Q/Dropbox/projects/forestlight/stanrdump'
years <- c(1990, 1995, 2000, 2005, 2010)
fgs <- c('fg1','fg2','fg3','fg4','fg5','unclassified')

for (i in 2:6) {
  prod_dump(alltreedat[[i]], to_file = TRUE, fn = file.path(fpdump, paste0('dump_alltree_', years[i-1], '.r')))
  for (j in 1:6) {
    prod_dump(fgdat[[j]][[i]], to_file = TRUE, fn = file.path(fpdump, paste0('dump_', fgs[j], '_', years[i-1], '.r')))
  }
}

