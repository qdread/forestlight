// Fit Weibull to density.

data {
  int<lower=0> N;
  real x[N];
}
parameters {
  real<lower=0> shape;
  real<lower=0> scale;
  real<lower=0> xtrue[N];
  real<lower=0> sigma;
}
model {
  // Priors
  shape ~ lognormal(1, 1);
  scale ~ lognormal(1, 1);

  // Likelihood
  xtrue ~ weibull(shape, scale);
  x ~ normal(xtrue, sigma);
}