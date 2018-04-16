// Stan model 06 Feb
// Pareto density
// Power law production

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
	// Power law production
	real<lower=0> beta0;
	real<lower=0> beta1;
	real<lower=0> sigma;
}

model {
	// Priors: Pareto density
	alpha ~ lognormal(1, 1) T[0, 5];
	// Priors: Power law production
	beta0 ~ normal(5, 2);
	beta1 ~ normal(0.5, 1);		
	sigma ~ exponential(0.01);
	
	// Likelihood: Pareto density
	x ~ pareto(x_min, alpha);
	// Likelihood: Power law production
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = -beta0 + beta1 * logx[i];
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	// Log likelihood and sum of the slopes
	// vector[N] log_lik_dens;
	// vector[N] log_lik_prod;
	real ee_slope;
	
	// for (i in 1:N) {
	// 	 log_lik_dens[i] = pareto_lpdf(x[i] | x_min, alpha);
	// 	 log_lik_prod[i] = normal_lpdf(logy[i] | -beta0 + beta1 * logx[i], sigma);
	// }
	
	ee_slope = beta1 - alpha - 1;
	
}
