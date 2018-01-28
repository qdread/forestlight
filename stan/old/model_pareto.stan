// Powerlaw for density fit (Pareto)

data {
  int<lower=0> N;
  real x[N];
}
parameters {
  real<lower=0, upper=20> xmin;
  real<lower=0, upper=5> alpha;
  real<lower=xmin> xtrue[N];
  real<lower=0> sigma;
}
model {
  // Priors
  xmin ~ lognormal(1, 1) T[0,20];
  alpha ~ lognormal(1, 1) T[0, 5];

  // Likelihood
  xtrue ~ pareto(xmin, alpha);
  x ~ normal(xtrue, sigma);
}