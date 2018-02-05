// Model of density and production
// Density as Weibull
// Production as power law
// Created on 05 Feb

data {
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
	real a;
	real<lower=0> b;
	real<lower=0> sigma;
}
model {
	// Priors
	alpha ~ lognormal(1, 1) T[0, 5];
	b ~ normal(0, 2);
	a  ~ normal(0, 10);
	sigma ~ exponential(0.01);

	// Likelihood
	x ~ weibull(shape, scale);
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = a + b * logx[i]
	  logy ~ normal(mu, sigma);
	}
}