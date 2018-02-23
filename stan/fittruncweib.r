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

# Define model with *truncated* weibull so that the log pdf is zero below the minimum diameter.
trunc_weib <- '
functions {
  real truncweib_log (real x, real m, real n, real xmin) {
    if (x < xmin) {
      return 0;
    } else {
      return log(m / n) + (m - 1) * log(x / n) - (x / n)^m;
    }
  }
}

data {
  int<lower=0> N;
  vector<lower=0>[N] x;
}

transformed data {
  real<lower=0> x_min;
  x_min = min(x);
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
    x[i] ~ truncweib(mu, nu, x_min);
  }
}
'

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
  nu ~ lognormal(1, 1);
  
  // Likelihood: Weibull density
  for (i in 1:N) {
    x[i] ~ weibull(mu, nu) T[x_min,x_max];
  }
}
'

trunctwice_weib_manual <- '
functions {
  real tweib_log (real x, real m, real n, real xmin, real xmax, real k) {
    if (x < xmin || x > xmax) {
      return negative_infinity();
    } else {
      return log(1/k) + (m-1)*log(x/n) - (x/n)^m; 
    }
  }
}

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

transformed parameters {
real K;
K = exp(-(x_min/nu)^mu) - exp(-(x_max/nu)^mu);
}

model {
// Priors: Weibull density
mu ~ lognormal(1, 1);
nu ~ lognormal(1, 1);

// Likelihood: Weibull density
for (i in 1:N) {
x[i] ~ tweib(mu, nu, x_min, x_max, K);
}
}
'

trunctwice_weib_manual <- '
functions {
real tweib_log (real x, real m, real n, real xmin, real xmax) {
if (x < xmin || x > xmax) {
return negative_infinity();
} else {
return (m-1)*log(x/n) - (x/n)^m; 
}
}
}

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
nu ~ lognormal(1, 1);

// Likelihood: Weibull density
for (i in 1:N) {
x[i] ~ tweib(mu, nu, x_min, x_max);
}
}
'

# Express the normalizing constant as a function of mu, nu, and xmin after doing the derivation longhand.

trunc_weib <- '
functions {
  real tweib_log (real x, real m, real n, real xmin) {
    if (x < xmin) {
      return negative_infinity();
    } else {
      return -m * log(m*n) + (xmin/n)^m + (m-1)*log(x/n) - (x/n)^m; // maybe add log(x/n) instead of log(x) in 3rd term
    }
  }
}

data {
int<lower=0> N;
vector<lower=0>[N] x;
}

transformed data {
real<lower=0> x_min;
x_min = min(x);
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
x[i] ~ tweib(mu, nu, x_min);
}
}
'

stanmodel_truncw <- stan_model(model_code = trunc_weib, model_name = 'trunc_weib')
stanmodel_trunc2w <- stan_model(model_code = trunctwice_weib, model_name = 'trunc2_weib')
stanmodel_trunc2manw <- stan_model(model_code = trunctwice_weib_manual, model_name = 'trunc2man_weib')

NC <- 1
NI <- 2000
NW <- 1500

fit_tw <- sampling(stanmodel_truncw, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_tw2 <- sampling(stanmodel_trunc2w, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_tw2 <- sampling(stanmodel_trunc2w, data = data1995_alltree_full, chains = NC, iter = NI, warmup = NW)
fit_tw2man <- sampling(stanmodel_trunc2manw, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
summary(fit_tw)
summary(fit_tw2)

library(bayesplot)
mcmc_trace(as.array(fit_tw))
mcmc_trace(as.array(fit_tw2))

# Plot fitted values.
densitybin_5census <- read.csv('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/densitybin_5census.csv', stringsAsFactors = FALSE)
dbh_pred_bins <- densitybin_5census[densitybin_5census$fg == 'all', 'bin_midpoint']

pdf_weibull <- function(x, m, n) (m/n) * (x/n)^(m-1) * exp(-(x/n)^m)
w_fitted <- pdf_weibull(dbh_pred_bins, 0.419, 0.585)

min_n <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/stanoutput/min_n.csv', stringsAsFactors = FALSE)

n_i <- min_n$n[min_n$year == 1995 & min_n$fg == 'alltree']
xmin_i <- min_n$xmin[min_n$year == 1995 & min_n$fg == 'alltree']

area_core <- 42.84

w_fitted_n <- w_fitted * n_i / area_core

global_diam_xlimits <- c(1, 316)

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
