// Stan model for two part piecewise production

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
	real<lower=x_min, upper=x_max> tau_p_low;
	real<lower=tau_p_low, upper=x_max> tau_p_high;
	// Power law production
	real<lower=0> beta0;
	real<lower=0> beta1_low;
	real<lower=0> beta1_mid;
	real<lower=0> beta1_high;
	real<lower=0> sigma;
}

transformed parameters {
	vector[N] x2; // Indicator for whether x[i] > tau_p_low
	vector[N] x3; // Indicator for whether x[i] > tau_p_high
	
	for (i in 1:N) {
		if (x[i] < tau_p_low && x[i] <= tau_p_high) {
			x2[i] = 0;
		} else {
			x2[i] = 1;
		}
		if (x[i] < tau_p_high) {
			x3[i] = 0;
		} else {
			x3[i] = 1;
		}
	}
}

model {
	// Priors: Power law production
	beta0 ~ normal(5, 2);
	beta1_low ~ normal(0.5, 1);
	beta1_mid ~ normal(0.5, 1);
	beta1_high ~ normal(0.5, 1);
	sigma ~ exponential(0.01);
	
	// Likelihood: Power law production
	{
	  vector[N] mu;
	   
	  for (i in 1:N) {
		  mu[i] = -beta0 + beta1_low * logx[i] + beta1_mid * (logx[i] - log(tau_p_low)) * x2[i] + beta1_high * (logx[i] - log(tau_p_high)) * x3[i];
	  }
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	// Log likelihood
	vector[N] log_lik_prod;
	
	for (i in 1:N) {
		log_lik_prod[i] = normal_lpdf(logy[i] | -beta0 + beta1_low * logx[i] + beta1_mid * (logx[i] - log(tau_p_low)) * x2[i] + beta1_high * (logx[i] - log(tau_p_high)) * x3[i], sigma);
	}
	
}
