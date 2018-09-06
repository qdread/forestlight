// Stan model for square root production

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
	// Sqrt Power law production
	real beta0;
	real<lower=0> beta1;
	real<lower=0> sigma;
}

model {
	// Priors: Power law production
	beta0 ~ normal(0, 1);
	beta1 ~ normal(0.5, 1);
	sigma ~ exponential(0.01);
	
	// Likelihood: Power law production
	{
	  vector[N] mu;
	   
	  for (i in 1:N) {
		  mu[i] = sqrt( -beta0 + beta1 * logx[i] );
	  }
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	// Log likelihood
	vector[N] log_lik_prod;
	
	for (i in 1:N) {
		log_lik_prod[i] = normal_lpdf(logy[i] | sqrt( -beta0 + beta1 * logx[i] ), sigma);
	}
	
}
