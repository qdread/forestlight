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
# Edited 8 Feb to subsample an even smaller number

n_sub <- 5000
set.seed(574)

for (i in 2:6) {
  prod_dump(alltreedat[[i]], to_file = TRUE, fn = file.path(fpdump, paste0('ssdump_alltree_', years[i-1], '.r')), subsample = n_sub)
  for (j in 1:6) {
    prod_dump(fgdat[[j]][[i]], to_file = TRUE, fn = file.path(fpdump, paste0('ssdump_', fgs[j], '_', years[i-1], '.r')), subsample = n_sub)
  }
}

# Added 7 Feb. 2018
# Create initial values dumps for each of the parameters

fpdump <- 'C:/Users/Q/Dropbox/projects/forestlight/stanrdump'

pareto_inits <- list(alpha = 2)
weibull_inits <- list(shape = 2, scale = 1)
powerlaw_inits <- list(beta0 = 1, beta1 = 1, sigma = 0.1)
powerlawexp_inits <- list(a = 0.1, b = 0.1, c = 1)

init_ppow <- c(pareto_inits, powerlaw_inits)
init_wpow <- c(weibull_inits, powerlaw_inits)
init_pexp <- c(pareto_inits, powerlaw_inits, powerlawexp_inits)
init_wexp <- c(weibull_inits, powerlaw_inits, powerlawexp_inits)
with(init_ppow, stan_rdump(list = names(init_ppow), file = file.path(fpdump, 'init_ppow.R')))
with(init_wpow, stan_rdump(list = names(init_wpow), file = file.path(fpdump, 'init_wpow.R')))
with(init_pexp, stan_rdump(list = names(init_pexp), file = file.path(fpdump, 'init_pexp.R')))
with(init_wexp, stan_rdump(list = names(init_wexp), file = file.path(fpdump, 'init_wexp.R')))
