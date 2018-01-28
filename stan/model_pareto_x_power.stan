// Model of density and production
// Density as Pareto
// Production as power law

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
}