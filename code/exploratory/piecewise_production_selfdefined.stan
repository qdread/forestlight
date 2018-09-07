// Stan model for two part piecewise production
// With hinge function (self-defined with thanks to Gelman)

functions {
	real logistic_hinge(real x, real x0, real a, real b0, real b1, real delta) { 
	  real xdiff = x - x0;
	  return a + b0 * xdiff + (b1 - b0) * delta * log1p_exp(xdiff / delta);
	}
}

data {
	int<lower=0> N;
	vector<lower=0>[N] x;
	vector<lower=0>[N] y;
	real<lower=0> x_min;
	real<lower=0> x_max;
}

transformed data {
	vector[N] logx;
	vector[N] logy;
	logx = log(x);
	logy = log(y);
}

parameters {
	real<lower=0> x0;
	real a;
	real<lower=0> b0;
	real<lower=0> b1;
	real<lower=0> delta;
	real<lower=0> sigma;
}

model {
	// Priors: Hinged production
	b0 ~ lognormal(1, 1);
	b1 ~ lognormal(1, 1);
	delta ~ exponential(10);
	a ~ normal(0, 10);
	x0 ~ normal(0, 10);
	
	// Likelihood: Hinged production
	{
	  vector[N] mu;
	   
	  for (i in 1:N) {
		  mu[i] = logistic_hinge(logx[i], x0, a, b0, b1, delta);
	  }
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	// Log likelihood
	vector[N] log_lik_prod;
	
	for (i in 1:N) {
		log_lik_prod[i] = normal_lpdf(logy[i] | logistic_hinge(logx[i], x0, a, b0, b1, delta), sigma);
	}
	
}
