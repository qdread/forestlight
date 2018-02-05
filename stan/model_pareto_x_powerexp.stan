// Model of density and production
// Density as Pareto
// Production as power law times exponential
// Created on 05 Feb

data {
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
	real<upper=0> a1;
	real<lower=0> b1;
	real<lower=0> a;
	real<upper=0> b;
	real c;
	real<lower=0> sigma;
}
model {
	// Priors
	alpha ~ lognormal(1, 1) T[0, 5];
	a ~ normal(0, 5);
	b ~ normal(0, 1);
	c ~ normal(0, 10);
	b1 ~ normal(0, 2);
	a1  ~ normal(0, 10);
	sigma ~ exponential(0.01);

	// Likelihood
	x ~ pareto(x_min, alpha);
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = (a1 + b1 * logx[i]) * (a * logx[i] ^ b + c);
	  logy ~ normal(mu, sigma);
	}
}