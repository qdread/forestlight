// Bertalanffy model with log link function
// Used to model production per unit area for individual trees as a function of incoming light per unit area
			
data {
	int<lower=0> N;
	vector<lower=0>[N] x;	// Light per area
	vector<lower=0>[N] y;	// Production per area
}

transformed data {
	vector[N] log_y;
	log_y = log(y);
}

parameters {
	real<lower=0> G;
	real<lower=0> b1;
	real<lower=0> k;
	real<lower=0> sigma;
}

model {
	// Priors
	G ~ lognormal(1, 1);
	b1 ~ lognormal(1, 1);
	k ~ lognormal(1, 1);
		
	// Likelihood
	{
		vector[N] mu;
		for (i in 1:N) mu[i] = G * (1 - b1 * exp(-k * x[i])) ^ 3;
		log_y ~ normal(log(mu), sigma);
	}
}

/* UNCOMMENT THIS IF YOU WANT TO CALCULATE LOG-LIKELIHOOD
generated quantities {
	vector[N] log_lik;
	for (i in 1:N) log_lik[i] = normal_lpdf(log_y[i] | log(G * (1 - b1 * exp(-k * x[i])) ^ 3;), sigma);
} 
*/
