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
	real beta0;
	real<lower=0> beta1;
	real<lower=0> sigma;
}
model {
	// Priors
	shape ~ lognormal(1, 1);
	scale ~ lognormal(1, 1);
	beta0  ~ normal(0, 10);
	beta1 ~ normal(0, 2);
	sigma ~ exponential(0.01);

	// Likelihood
	x ~ weibull(shape, scale);
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = beta0 + beta1 * logx[i];
	  logy ~ normal(mu, sigma);
	}
}
