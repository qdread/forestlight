# Fit Weibull with truncation at minimum diam


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
data1995_alltree_full <- prod_dump(alltreedat[[3]], subsample = NULL)


trunctwice_weib <- '
data {
  int<lower=0> N;
  vector<lower=0>[N] x;
}

transformed data {
  real<lower=0> x_min;
  real<lower=0> x_max;
  x_min = min(x);
  x_max = max(x);
}

parameters {
  // Weibull density
  real<lower=0> mu;
  real<lower=0> nu;
}

model {
  // Priors: Weibull density
  mu ~ lognormal(1, 1);
  nu ~ lognormal(-1, 1);
  
  // Likelihood: Weibull density
  for (i in 1:N) {
    x[i] ~ weibull(mu, nu) T[x_min,x_max];
  }
}
'

trunctwice_weib <- '
data {
int<lower=0> N;
real<lower=0> UL;
real<lower=0> LL;
vector<lower=0>[N] x;
}

parameters {
// Weibull density
real<lower=0> mu;
real<lower=0> nu;
}

model {
// Priors: Weibull density
mu ~ lognormal(1, 1);
nu ~ lognormal(1, 1);

// Likelihood: Weibull density
for (i in 1:N) {
x[i] ~ weibull(mu, nu) T[UL, LL];
}
}
'

stanmodel_trunc2w <- stan_model(model_code = trunctwice_weib, model_name = 'trunc2_weib')

NC <- 1
NI <- 2000
NW <- 1500

lowlim <- 1.1
upplim <- 316

fit_tw2 <- sampling(stanmodel_trunc2w, data = c(data1995_alltree, UL=lowlim, LL=upplim), chains = NC, iter = NI, warmup = NW)
fit_tw2 <- sampling(stanmodel_trunc2w, data = data1995_alltree_full, chains = NC, iter = NI, warmup = NW)
summary(fit_tw2)

library(bayesplot)
mcmc_trace(as.array(fit_tw2))

# Plot fitted values.
densitybin_5census <- read.csv('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/densitybin_5census.csv', stringsAsFactors = FALSE)
dbh_pred_bins <- densitybin_5census[densitybin_5census$fg == 'all', 'bin_midpoint']

w_fitted <- dweibull(dbh_pred_bins, 0.428, 0.657)


min_n <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/stanoutput/min_n.csv', stringsAsFactors = FALSE)

n_i <- min_n$n[min_n$year == 1995 & min_n$fg == 'alltree']
xmin_i <- min_n$xmin[min_n$year == 1995 & min_n$fg == 'alltree']

area_core <- 42.84

# Remove truncations
upper_lower <- pweibull(c(lowlim,upplim), 0.428, 0.657)

w2 <- w_fitted / sum(1-upper_lower)
w_fitted_n <- w2/sum(w2) * n_i/area_core
w_fitted_n <- (w_fitted / sum(w_fitted) / sum(1-upper_lower)) * n_i / area_core

global_diam_xlimits <- c(1, 316)

library(dplyr)
library(cowplot)

p_all_dens <- all_dens_1995 %>%
  filter(!is.na(bin_value), bin_value > 0) %>%
  mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
  mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
  ggplot() +
  geom_point(aes(x = bin_midpoint, y = bin_yvalue), color = 'black') +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  panel_border(colour = 'black')


p_all_dens + geom_line(data = data.frame(x = dbh_pred_bins, y = w_fitted_n), aes(x,y))
