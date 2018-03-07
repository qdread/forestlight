# Fit stan models for 1995.
# Edit 06 Feb 2018: Use the new models and some test data perhaps.
# Edit 13 Feb 2018: Add the new specification of Weibulls

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
  xdat <- list(N = length(x), x = x, y = y, x_min = min(x), LL = 1.1, UL = 316) # Truncation limits from data
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

dbh_pred_100 <- logbin(x = allyeardbh, y = NULL, n = 100)$bin_mid

ci_ppow <- dens_prod_ci(fit_ppow_all, seq(1.1,315,length.out=100), 'pareto', 'powerlaw', min_n$xmin[3], min_n$n[3])
ci_pexp <- dens_prod_ci(fit_pexp_all, seq(1.1,315,length.out=100), 'pareto', 'powerlawexp', min_n$xmin[3], min_n$n[3])
ci_wpow <- dens_prod_ci(fit_wpow_all, seq(1.1,315,length.out=100), 'weibull', 'powerlaw', min_n$xmin[3], min_n$n[3])
ci_wexp <- dens_prod_ci(fit_wexp_all, seq(1.1,315,length.out=100), 'weibull', 'powerlawexp', min_n$xmin[3], min_n$n[3])

ci_df <- cbind(fg = 'alltree', dens_model = rep(c('pareto','weibull'), each=nrow(ci_ppow)*2), prod_model = rep(c('powerlaw','powerlawexp'), each=nrow(ci_ppow)), rbind(ci_ppow, ci_pexp, ci_wpow, ci_wexp))

save(fit_ppow_all, fit_pexp_all, fit_wpow_all, fit_wexp_all, ci_df, file = 'C:/Users/Q/Dropbox/projects/forestlight/stanoutput/localsubsamplefit5000_05Mar_alltree.RData')

fit_ppow_fg <- lapply(data1995_byfg, function(x) sampling(stanmodel_paretoxpower, data = x, chains = NC, iter = NI, warmup = NW))
fit_pexp_fg <- lapply(data1995_byfg, function(x) sampling(stanmodel_paretoxexp, data = x, chains = NC, iter = NI, warmup = NW))
fit_wpow_fg <- lapply(data1995_byfg, function(x) sampling(stanmodel_weibullxpower, data = x, chains = NC, iter = NI, warmup = NW))
fit_wexp_fg <- lapply(data1995_byfg, function(x) sampling(stanmodel_weibullxexp, data = x, chains = NC, iter = NI, warmup = NW))

save(list = grep('fit_', ls(), value = TRUE), file = 'C:/Users/Q/Dropbox/projects/forestlight/stanoutput/localsubsamplefit5000_13feb.RData')

###
library(rstan)
load('C:/Users/Q/Dropbox/projects/forestlight/stanoutput/localsubsamplefit5000_13feb.RData')

summary(fit_ppow_all)[[1]]
summary(fit_pexp_all)[[1]]
summary(fit_wpow_all)[[1]]
summary(fit_wexp_all)[[1]]


library(bayesplot)

mcmc_trace(as.array(fit_ppow_all))
mcmc_trace(as.array(fit_pexp_all))
mcmc_trace(as.array(fit_wpow_all))
mcmc_trace(as.array(fit_wexp_all))
