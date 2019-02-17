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
	
	// Loglinear production
	real<lower=0> beta0; // Intercept
	real<lower=0> beta1; // Slope
	real<lower=0> sigma;
}

model {
	// Prior: Pareto density
	alpha ~ lognormal(1, 1) T[0, 5];
	
	// Priors: Loglinear production
	beta0 ~ lognormal(1, 1);
	beta1 ~ lognormal(1, 1);
	
	sigma ~ exponential(0.1);
	
	// Likelihood: Pareto density
	x ~ pareto(x_min, alpha);
	
	// Likelihood: Loglinear production
	{
	  vector[N] mu;
	   
	  for (i in 1:N) {
		  mu[i] = log(beta0) + beta1 * logx[i];
	  }
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	vector[N] log_lik_dens; // Log-likelihood for getting info criteria later
	vector[N] log_lik_prod;
	
	for (i in 1:N) {
		log_lik_dens[i] = pareto_lpdf(x[i] | x_min, alpha);
		log_lik_prod[i] = normal_lpdf(logy[i] | log(beta0) + beta1 * logx[i], sigma);
	}
}
