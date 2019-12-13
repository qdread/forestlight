functions {
	real logistic_hinge(real x, real x0, real beta0, real beta1_low, real beta1_high, real delta) { 
		real xdiff = x - log(x0);
		return log(beta0) + beta1_low * xdiff + (beta1_high - beta1_low) * delta * log1p_exp(xdiff / delta);
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
	// Pareto density
	real<lower=0, upper=5> alpha;
	
	// Hinged production
	real<lower=0> x0;
	real<lower=0> beta0;
	real<lower=0> beta1_low;
	real<lower=0> beta1_high;
	real<lower=0> delta;
	real<lower=0> sigma;
}

model {
	// Prior: Pareto density
	alpha ~ lognormal(1, 1) T[0, 5];
	
	// Priors: Hinged production
	beta0 ~ lognormal(1, 1);
	beta1_low ~ lognormal(1, 1);
	beta1_high ~ lognormal(1, 1);
	delta ~ exponential(10);
	x0 ~ lognormal(1, 1);
	
	sigma ~ exponential(0.1);
	
	// Likelihood: Pareto density
	x ~ pareto(x_min, alpha);
	
	// Likelihood: hinged production
	{
	  vector[N] mu;
	   
	  for (i in 1:N) {
		  mu[i] = logistic_hinge(logx[i], x0, beta0, beta1_low, beta1_high, delta);
	  }
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	vector[N] log_lik_dens; // Log-likelihood for getting info criteria later
	vector[N] log_lik_prod;
	
	for (i in 1:N) {
		log_lik_dens[i] = pareto_lpdf(x[i] | x_min, alpha);
		log_lik_prod[i] = normal_lpdf(logy[i] | logistic_hinge(logx[i], x0, beta0, beta1_low, beta1_high, delta), sigma);
	}
}
