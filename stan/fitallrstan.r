# Fit stan models for 1995.
# Edit 06 Feb 2018: Use the new models and some test data perhaps.

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'
load(file.path(fpdata, 'rawdataobj_22jan.r'))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

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

n_sub <- 5000
set.seed(574)

data1995_alltree <- prod_dump(alltreedat[[3]], subsample = n_sub)
data1995_byfg <- lapply(fgdat, function(x) prod_dump(x[[3]], subsample = n_sub))

stanmodel_paretoxpower <- stan_model(file = 'stan/model_ppow.stan', model_name = 'paretoxpow')
stanmodel_paretoxexp <- stan_model(file = 'stan/model_pexp.stan', model_name = 'paretoxexp')
stanmodel_weibullxpower <- stan_model(file = 'stan/model_wpow.stan', model_name = 'weibullxpow')
stanmodel_weibullxexp <- stan_model(file = 'stan/model_wexp.stan', model_name = 'weibullxexp')

NC <- 3
NI <- 6000
NW <- 5000

fit_ppow_all <- sampling(stanmodel_paretoxpower, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_pexp_all <- sampling(stanmodel_paretoxexp, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_wpow_all <- sampling(stanmodel_weibullxpower, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_wexp_all <- sampling(stanmodel_weibullxexp, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)

fit_ppow_fg <- lapply(data1995_byfg, function(x) sampling(stanmodel_paretoxpower, data = x, chains = NC, iter = NI, warmup = NW))
fit_pexp_fg <- lapply(data1995_byfg, function(x) sampling(stanmodel_paretoxexp, data = x, chains = NC, iter = NI, warmup = NW))
fit_wpow_fg <- lapply(data1995_byfg, function(x) sampling(stanmodel_weibullxpower, data = x, chains = NC, iter = NI, warmup = NW))
fit_wexp_fg <- lapply(data1995_byfg, function(x) sampling(stanmodel_weibullxexp, data = x, chains = NC, iter = NI, warmup = NW))

save(list = grep('fit_', ls(), value = TRUE), file = 'C:/Users/Q/Dropbox/projects/forestlight/stanoutput/localsubsamplefit5000.RData')


# summary(fit_ppow)
# summary(fit_pexp)
# summary(fit_wpow)
# summary(fit_wexp)
# 
# library(bayesplot)
# 
# mcmc_trace(as.array(fit_ppow))
# mcmc_trace(as.array(fit_pexp))
# mcmc_trace(as.array(fit_wpow))
# mcmc_trace(as.array(fit_wexp))

# fit_pareto_byfg <- lapply(data1995_byfg, function(x) sampling(stanmodel_paretoxpower,
#                                                               data = x,
#                                                               chains = 3,
#                                                               iter = 10000,
#                                                               warmup = 5000))
# 
# fit_weibull_byfg <- lapply(data1995_byfg, function(x) sampling(stanmodel_weibullxexp,
#                                                                data = x,
#                                                                chains = 3,
#                                                                iter = 10000,
#                                                                warmup = 5000))
# 
# save(fit_pareto_byfg, fit_weibull_byfg, file = 'C:/Users/Q/Dropbox/projects/forestlight/stanoutput/paretoweibullfits.RData')