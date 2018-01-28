// Hierarchical model including both density and production
// Density as Pareto and production as strict power law
// If the exponents add to zero, should get energy equivalence

data {
  int<lower=0> N;
  int<lower=0> N_pred;
  real x[N];
}
parameters {
  real<lower=0, upper=20> xmin;
  real<lower=0, upper=5> alpha;
  real<lower=xmin> xtrue[N];
  real<lower=0> sigma;
}
model {
	
  // DENSITY MODEL (Pareto)  
  
  // Priors
  xmin ~ lognormal(1, 1) T[0,20];
  alpha ~ lognormal(1, 1) T[0, 5];

  // Likelihood
  xtrue ~ pareto(xmin, alpha);
  x ~ normal(xtrue, sigma);
  
  // PRODUCTION MODEL (Power law)
  
  b ~ normal(2, 1);
  a  ~ normal(0, 10);
  sigma ~ exponential(0.01);
  {
	  vector[N] mu;
	  for (i in 1:N) mu[i] = a + b * logx[i];
	  logy ~ normal(mu, sigma);
  }
    
}

generated quantities {
	// Edit this to now predict total production instead of individual production.
	vector[N_pred] y_pred;
	for (i in 1:N_pred) y_pred[i] = 10^(a + b * logx_pred[i]);
}