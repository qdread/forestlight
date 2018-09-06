// Stan model for two part piecewise production
// With self defined function

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
	real<lower=0> beta0;
	real<lower=0> beta1;
	real<lower=0> a;
	real<lower=0> b;
	real c;
	real<lower=0> sigma;
}

model {
	// Priors: Power law production
	beta0 ~ normal(5, 2);
	beta1 ~ normal(0.5, 1);		
	beta1_high ~ normal(0.5, 1);
	sigma ~ exponential(0.01);
	
	// Likelihood: Power law production
	{
	  vector[N] mu;
	   
	  for (i in 1:N) {
		  if (x <= tau_p) {
			  mu[i] = ;
		  } else {
			  mu[i] = ;
		  }
	  }
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	// Log likelihood
	vector[N] log_lik_prod;
	
	for (i in 1:N) {
		log_lik_prod[i] = normal_lpdf(logy[i] | , sigma);
	}
	
}
