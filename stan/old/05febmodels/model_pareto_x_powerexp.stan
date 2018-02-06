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
	real<upper=0> beta0;
	real<lower=0> beta1;
	real<lower=0> a;
	real<upper=0> b;
	real c;
	real<lower=0> sigma;
}
model {
	// Priors
	alpha ~ lognormal(1, 1) T[0, 5];
	a ~ normal(0, 1);
	b ~ normal(0, 5);
	c ~ normal(0, 1);
	beta0  ~ normal(0, 10);
	beta1 ~ normal(0, 2);
	sigma ~ exponential(0.01);

	// Likelihood
	x ~ pareto(x_min, alpha);
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = (beta0 + beta1 * logx[i]) * (a * x[i] ^ b + c);
	  logy ~ normal(mu, sigma);
	}
}
