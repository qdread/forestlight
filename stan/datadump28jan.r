# Create Rdump for new grouped data to run with CMDSTAN.

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'
load(file.path(fpdata, 'rawdataobj_22jan.r'))

prod_dump <- function(dat, to_file = FALSE, fn = NULL, subsample = NULL) {
  require(rstan)
  if (!is.null(subsample) && nrow(dat) > subsample) {
    dat <- dat[sample(nrow(dat), subsample, replace = FALSE), ]
  }
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

# Get minimum x values and numbers of individuals that are used for plotting.
valall <- do.call('rbind', lapply(alltreedat, function(x) data.frame(xmin = min(x$dbh_corr), n = nrow(x))))
valall <- data.frame(fg = 'alltree', year = c(1985, 1990, 1995, 2000, 2005, 2010), valall)
valfg <- do.call('rbind', do.call('rbind', lapply(fgdat, function(x) lapply(x, function(x) data.frame(xmin = min(x$dbh_corr), n = nrow(x))))))

valfg <- data.frame(fg = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
                    year = rep(c(1985, 1990, 1995, 2000, 2005, 2010), each = 6),
                    valfg)

min_n <- rbind(valall, valfg)
write.csv(min_n, 'C:/Users/Q/Dropbox/projects/forestlight/stanoutput/min_n.csv', row.names = FALSE)

# Added 3 Feb. 2018
# Take a subset of the data so that we can fit the Weibull in a reasonable amount of time.

n_sub <- 25000
set.seed(574)

for (i in 2:6) {
  prod_dump(alltreedat[[i]], to_file = TRUE, fn = file.path(fpdump, paste0('ssdump_alltree_', years[i-1], '.r')), subsample = n_sub)
  for (j in 1:6) {
    prod_dump(fgdat[[j]][[i]], to_file = TRUE, fn = file.path(fpdump, paste0('ssdump_', fgs[j], '_', years[i-1], '.r')), subsample = n_sub)
  }
}