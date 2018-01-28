// Model of density and production
// Density as Weibull
// Production as power law times exponential

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
}