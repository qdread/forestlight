// "Modules" for different models in Stan
// Edited 09 Feb 2018: weakly informative priors.

// Pareto density

data {
	int<lower=0> N;
	vector<lower=0>[N] x;
	vector<lower=0>[N] y;
	real<lower=0> x_min;
}
transformed data {
	vector[N] logx;
	vector[N] logy;
	logx = log(x);
	logy = log(y);
}
parameters {
	real<lower=0, upper=5> alpha;
}
model {
	// Priors
	alpha ~ lognormal(1, 1) T[0, 5];

	// Likelihood
	x ~ pareto(x_min, alpha);
}

// Weibull density

data {
	int<lower=0> N;
	vector<lower=0>[N] x;
	vector<lower=0>[N] y;
}
//same transformed data as pareto
parameters {
	real<lower=0> shape;
	real<lower=0> scale;
}
model {
	// Priors
	shape ~ lognormal(1, 1);
	scale ~ lognormal(1, 1);
	// Likelihood
	x ~ weibull(shape, scale);
}

// Power law production

parameters {
	real<lower=0> beta0;
	real<lower=0> beta1;
	real<lower=0> sigma;
}
model {
	// Priors
	beta0  ~ normal(5, 2);
	beta1 ~ normal(0.5, 1);	
	sigma ~ exponential(0.01);

	// Likelihood
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = -beta0 + beta1 * logx[i];
	  logy ~ normal(mu, sigma);
	}
}

// Power law times exponential production

parameters {
	real<lower=0> beta0;
	real<lower=0> beta1;
	real<lower=0> a;
	real<lower=0> b;
	real c;
	real<lower=0> sigma;
}
model {
	// Priors
	a ~ normal(5, 5);
	b ~ normal(0.5, 1);
	c ~ normal(5, 10);
	beta0 ~ normal(5, 2);
	beta1 ~ normal(0.5, 1);		
	sigma ~ exponential(0.01);

	// Likelihood
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = -beta0 + beta1 * logx[i] + log(-a * x[i] ^ -b + c);
	  logy ~ normal(mu, sigma);
	}
}