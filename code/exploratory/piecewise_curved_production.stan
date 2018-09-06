// Stan model for two part piecewise production
// Lower portion is curved

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
	real<lower=x_min, upper=x_max> tau_p;
	// Power law production
	real<lower=0> beta1;
	real<lower=0> a;
	real<lower=0> b;
	real c;
	real<lower=0> sigma;
}

transformed parameters {
	vector[N] x2; // Indicator for whether x[i] > tau_p
	
	for (i in 1:N) {
		if (x[i] < tau_p) {
			x2[i] = 0;
		} else {
			x2[i] = 1;
		}
	}
}

model {
	// Priors: Power law production
	a ~ normal(5, 5);
	b ~ normal(0.5, 1);
	c ~ normal(5, 10);
	beta1 ~ normal(0.5, 1);		
	sigma ~ exponential(0.01);
	
	// Likelihood: Power law production
	{
	  vector[N] mu;
	   
	  for (i in 1:N) {
		  mu[i] = log(-a * x[i] ^ -b + c) + beta1 * (logx[i] - log(tau_p)) * x2[i];
	  }
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	// Log likelihood
	vector[N] log_lik_prod;
	
	for (i in 1:N) {
		log_lik_prod[i] = normal_lpdf(logy[i] | log(-a * x[i] ^ -b + c) + beta1 * (logx[i] - log(tau_p)) * x2[i], sigma);
	}
	
}
