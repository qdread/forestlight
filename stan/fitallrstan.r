# Fit stan models for 1995.
# Edit 06 Feb 2018: Use the new models and some test data perhaps.

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'
load(file.path(fpdata, 'rawdataobj_22jan.r'))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

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

data1995_alltree <- prod_dump(alltreedat[[3]])
data1995_byfg <- lapply(fgdat, function(x) prod_dump(x[[3]]))

stanmodel_paretoxpower <- stan_model(file = 'stan/model_ppow.stan', model_name = 'paretoxpow')
stanmodel_paretoxexp <- stan_model(file = 'stan/model_pexp.stan', model_name = 'paretoxexp')
stanmodel_weibullxpower <- stan_model(file = 'stan/model_wpow.stan', model_name = 'weibullxpow')
stanmodel_weibullxexp <- stan_model(file = 'stan/model_wexp.stan', model_name = 'weibullxexp')

# Test each model.
fit_ppow <- sampling(stanmodel_paretoxpower, data = data1995_byfg[[1]], chains = 2, iter = 1500, warmup = 1000)
fit_pexp <- sampling(stanmodel_paretoxexp, data = data1995_byfg[[1]], chains = 2, iter = 1500, warmup = 1000)
fit_wpow <- sampling(stanmodel_weibullxpower, data = data1995_byfg[[1]], chains = 2, iter = 1500, warmup = 1000)
fit_wexp <- sampling(stanmodel_weibullxexp, data = data1995_byfg[[1]], chains = 2, iter = 1500, warmup = 1000)

summary(fit_ppow)
summary(fit_pexp)
summary(fit_wpow)
summary(fit_wexp)

library(bayesplot)

mcmc_trace(as.array(fit_ppow))
mcmc_trace(as.array(fit_pexp))
mcmc_trace(as.array(fit_wpow))
mcmc_trace(as.array(fit_wexp))

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