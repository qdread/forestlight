// Stan model 06 Feb
// Pareto density
// Power law times exponential production

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
	// Pareto density
	real<lower=0, upper=5> alpha;
	// Power law times exponential production
	real<lower=0> beta0;
	real<lower=0> beta1;
	real<lower=0> a;
	real<lower=0> b;
	real c;
	real<lower=0> sigma;
}

model {
	// Priors: Pareto density
	alpha ~ lognormal(1, 1) T[0, 5];
	// Priors: Power law times exponential production
	a ~ normal(5, 5);
	b ~ normal(0.5, 1);
	c ~ normal(5, 10);
	beta0 ~ normal(5, 2);
	beta1 ~ normal(0.5, 1);		
	sigma ~ exponential(0.01);
	
	// Likelihood: Pareto density
	x ~ pareto(x_min, alpha);
	// Likelihood: Power law times exponential production
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = -beta0 + beta1 * logx[i] + log(-a * x[i] ^ -b + c);
	  logy ~ normal(mu, sigma);
	}
}
