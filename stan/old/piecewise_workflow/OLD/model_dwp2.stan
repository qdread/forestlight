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
	real<lower=0> LL; // Lower truncation limit
	real<lower=0> UL; // Upper truncation limit
}

transformed data {
	vector[N] logx;
	vector[N] logy;
	logx = log(x);
	logy = log(y);
}

parameters {
	// Weibull density
	real<lower=0> m;
	real<lower=0> n;
	
	// Hinged production
	real<lower=0> x0;
	real<lower=0> beta0;
	real<lower=0> beta1_low;
	real<lower=0> beta1_high;
	real<lower=0> delta;
	real<lower=0> sigma;
}

model {
	// Priors: Weibull density
	m ~ lognormal(1, 1);
	n ~ lognormal(1, 1);
	
	// Priors: Hinged production
	beta0 ~ lognormal(1, 1);
	beta1_low ~ lognormal(1, 1);
	beta1_high ~ lognormal(1, 1);
	delta ~ exponential(10);
	x0 ~ lognormal(1, 1);
	
	sigma ~ exponential(0.1);
	
	// Likelihood: Weibull density
	for (i in 1:N) x[i] ~ weibull(m, n) T[LL,UL];
	
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
	real k;
	
	k = -log_diff_exp(weibull_lcdf(UL | m, n), weibull_lcdf(LL | m, n));
	
	for (i in 1:N) {
		log_lik_dens[i] = weibull_lpdf(x[i] | m, n) + k;
		log_lik_prod[i] = normal_lpdf(logy[i] | logistic_hinge(logx[i], x0, beta0, beta1_low, beta1_high, delta), sigma);
	}
}
