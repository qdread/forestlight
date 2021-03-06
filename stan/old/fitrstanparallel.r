# Fit stan models for all years with full datasets
# Use RStan on HPCC
# 16 Feb 2018
# Edited 27 March: add models that output log likelihoods with each iteration.

NC <- 3
NI <- 4000
NW <- 3500

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

z <- expand.grid(year = c(1990, 1995, 2000, 2005, 2010),
				 fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified', 'alltree'),
				 model_name = c('ppow', 'wpow', 'pexp', 'wexp'),
				 stringsAsFactors = FALSE)

# Load data (changing file path if necessary)
fpdata <- '~/forestlight/stanrdump'
load(file.path(fpdata, 'rawdataobj_22jan.r'))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = NC) 

prod_dump <- function(dat, to_file = FALSE, fn = NULL, subsample = NULL) {
  require(rstan)
  if (!is.null(subsample) && nrow(dat) > subsample) {
    dat <- dat[sample(nrow(dat), subsample, replace = FALSE), ]
  }
  x <- dat$dbh_corr
  y <- dat$production
  xdat <- list(N = length(x), x = x, y = y, x_min = min(x), LL = 1.1, UL = 316)
  if (to_file) {
    with(xdat, stan_rdump(names(xdat), file = fn))
  } else {
    return(xdat)
  }
}

year_idx <- which(z$year[task] == c(1985, 1990, 1995, 2000, 2005, 2010))

if (z$fg[task] == 'alltree') {
	dat <- prod_dump(alltreedat[[year_idx]])
} else {
	fg_idx <- which(z$fg[task] == c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified'))
	dat <- prod_dump(fgdat[[fg_idx]][[year_idx]])
}

model_file <- c('model_ppow_withlik.stan', 'model_wpow_withlik.stan', 'model_pexp_withlik.stan', 'model_wexp_withlik.stan')[which(z$model_name[task] == c('ppow', 'wpow', 'pexp', 'wexp'))]
stan_model_dens_prod <- stan_model(file = file.path('~/forestlight/stancode', model_file), model_name = z$model_name[task])

fit <- sampling(stan_model_dens_prod, data = dat, chains = NC, iter = NI, warmup = NW)

save(fit, file = file.path('~/forestlight/stanoutput', paste0('fit_', z$model_name[task], '_', z$fg[task], '_', z$year[task], '.r')))


