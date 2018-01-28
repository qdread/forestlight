# Fit stan models for 1995.

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/'
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

# Define models
model_pareto_x_powerlaw <- 
  'data {
int<lower=0> N;
vector[N] x;
vector[N] y;
real<lower=0> x_min;
}
transformed data {
vector[N] logx;
vector[N] logy;
logx = log(x)/log(10);
logy = log(y)/log(10);
}
parameters {
real<lower=0, upper=5> alpha;
real a;
real b;
real<lower=0> sigma;
}
model {
// Priors
alpha ~ lognormal(1, 1) T[0, 5];
b ~ normal(2, 1);
a  ~ normal(0, 10);
sigma ~ exponential(0.01);

// Likelihood
x ~ pareto(x_min, alpha);
{
  vector[N] mu;
  for (i in 1:N) mu[i] = a + b * logx[i];
  logy ~ normal(mu, sigma);
}
}'

model_weibull_x_powerexp <- 
  'data {
int<lower=0> N;
vector[N] x;
vector[N] y;
}
transformed data {
vector[N] logx;
vector[N] logy;
logx = log(x)/log(10);
logy = log(y)/log(10);
}
parameters {
real<lower=0> shape;
real<lower=0> scale;
real a1;
real b1;
real a;
real b;
real c;
real<lower=0> sigma;
}
model {
// Priors
shape ~ lognormal(1, 1);
scale ~ lognormal(1, 1);
a ~ normal(1, 2);
b ~ normal(0, 10);
c ~ normal(1, 2);
b1 ~ normal(2, 1);
a1  ~ normal(0, 10);
sigma ~ exponential(0.01);

// Likelihood
x ~ weibull(shape, scale);
{
  vector[N] mu;
  for (i in 1:N) mu[i] = (a1 + b1 * logx[i]) * (a * logx[i] ^ b + c);
  logy ~ normal(mu, sigma);
}
}'

stanmodel_paretoxpower <- stan_model(model_code = model_pareto_x_powerlaw)
stanmodel_weibullxexp <- stan_model(model_code = model_weibull_x_powerexp)

fit_pareto_byfg <- lapply(data1995_byfg, function(x) sampling(stanmodel_paretoxpower,
                                                              data = x,
                                                              chains = 3,
                                                              iter = 10000,
                                                              warmup = 5000))

fit_weibull_byfg <- lapply(data1995_byfg, function(x) sampling(stanmodel_weibullxexp,
                                                               data = x,
                                                               chains = 3,
                                                               iter = 10000,
                                                               warmup = 5000))

save(fit_pareto_byfg, fit_weibull_byfg, file = 'C:/Users/Q/Dropbox/projects/forestlight/stanoutput/paretoweibullfits.RData')