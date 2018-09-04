// Stan model 06 Feb
// Pareto density
// Power law production

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
	real<lower=0> beta1_low;
	real<lower=0> beta1_high;
	real<lower=0> sigma;
}

model {
	// Priors: Power law production
	beta0 ~ normal(5, 2);
	beta1 ~ normal(0.5, 1);		
	sigma ~ exponential(0.01);
	
	// Likelihood: Power law production
	{
	  vector[N] mu;
	  real<lower=0> C_con;
	  
	  C_con = beta0 * tau_p ^ (beta1_high - beta1_low);
	  
	  for (i in 1:N) {
		  if (x[i] <= tau_p) mu[i] = ;
		  if (x[i] > tau_p) mu[i] = ;
	  }
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	// Log likelihood
	vector[N] log_lik_prod;
	
	for (i in 1:N) {
		log_lik_prod[i] = normal_lpdf(logy[i] | -beta0 + beta1 * logx[i], sigma);
	}
	
}
