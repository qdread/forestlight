# Fit piecewise function to density
# 20 Feb 2018
# At first do nothing to constrain the two functions meeting at break point

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
  #xdat <- list(N = length(x), x = x, y = y, x_min = min(x))
  xdat <- list(N = length(x), x = x)
  if (to_file) {
    with(xdat, stan_rdump(names(xdat), file = fn))
  } else {
    return(xdat)
  }
}

n_sub <- 5000
set.seed(517)

data1995_alltree <- prod_dump(alltreedat[[3]], subsample = n_sub)

stanmodel_dens_we <- stan_model(file = 'stan/model_dens_wexp.stan', model_name = 'dens_we')
stanmodel_dens_pe <- stan_model(file = 'stan/model_dens_pexp.stan', model_name = 'dens_pe')
stanmodel_dens_pe2 <- stan_model(file = 'stan/model_dens_pexpmatch.stan', model_name = 'dens_pe_match')
stanmodel_dens_we2 <- stan_model(file = 'stan/model_dens_wexpmatch.stan', model_name = 'dens_we_match')

NC <- 2
NI <- 2000
NW <- 1500

fit_dwe <- sampling(stanmodel_dens_we, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_dpe <- sampling(stanmodel_dens_pe, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_dpe2 <- sampling(stanmodel_dens_pe2, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_dwe2 <- sampling(stanmodel_dens_we2, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)

summary(fit_dwe)
summary(fit_dpe)
summary(fit_dpe2)
summary(fit_dwe2)

library(bayesplot)
mcmc_trace(as.array(fit_dwe))
mcmc_trace(as.array(fit_dpe))
mcmc_trace(as.array(fit_dpe2))
mcmc_trace(as.array(fit_dwe2))

# Generate fitted values and uncertainty
ci_piecewise <- function(fit, func, dbh_pred, x_min = NULL, n_indiv = 1) {
  pdf_paretoexp <- function(x, xmin, alpha, L, beta) {
    ifelse(x <= L, (alpha * xmin^alpha) / (x ^ (alpha+1)), beta * exp(-beta * x))
  }
  pdf_weibullexp <- function(x, m, n, L, beta) {
    ifelse(x <= L, (m/n) * (x/n)^(m-1) * exp(-(x/n)^m), beta * exp(-beta * x))
  }
  
  pars <- as.data.frame(do.call('cbind', extract(fit)))
  
  if (func == 'weibullexp') {
    dens_pred <- with(pars, sapply(dbh_pred, pdf_weibullexp, m = m, n = n, L = L, beta = beta))
  }
  if (func == 'paretoexp') {
    dens_pred <- with(pars, sapply(dbh_pred, pdf_paretoexp, xmin = x_min, alpha = alpha, L = L, beta = beta))
  }
  
  dens_pred <- dens_pred * n_indiv
  
  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  dens_pred_quant <- apply(dens_pred, 2, quantile, probs = qprobs)
  
  out <- data.frame(dbh = dbh_pred,
                    t(dens_pred_quant))

  setNames(out, c('dbh', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
}

# Confidence interval plots.
# Load density bin midpoints to get x values at which to predict.
densitybin_5census <- read.csv('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/densitybin_5census.csv', stringsAsFactors = FALSE)
dbh_pred_bins <- densitybin_5census[densitybin_5census$fg == 'all', 'bin_midpoint']

# Minimum x values for Pareto.
min_n <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/stanoutput/min_n.csv', stringsAsFactors = FALSE)

n_i <- min_n$n[min_n$year == 1995 & min_n$fg == 'alltree']
xmin_i <- min_n$xmin[min_n$year == 1995 & min_n$fg == 'alltree']

ci_we_1995 <- ci_piecewise(fit = fit_dwe, func = 'weibullexp', dbh_pred = dbh_pred_bins, n_indiv = n_i)
ci_pe_1995 <- ci_piecewise(fit = fit_dpe2, func = 'paretoexp', dbh_pred = dbh_pred_bins, x_min = xmin_i, n_indiv = n_i)

par(mfrow=c(1,2))

with(ci_we_1995, plot(log10(dbh), log10(q975), type = 'l', lty = 3, ylab = 'Weibull piecewise'))
with(ci_we_1995, lines(log10(dbh), log10(q025), type = 'l', lty = 3))
with(ci_we_1995, lines(log10(dbh), log10(q50), type = 'l', lty = 1))

with(ci_pe_1995, plot(log10(dbh), log10(q975), type = 'l', lty = 3, ylab = 'Pareto piecewise'))
with(ci_pe_1995, lines(log10(dbh), log10(q025), type = 'l', lty = 3))
with(ci_pe_1995, lines(log10(dbh), log10(q50), type = 'l', lty = 1))
