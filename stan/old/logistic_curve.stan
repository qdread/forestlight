// Logistic model with log link function
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
	real<lower=min(y), upper=max(y)> L;
	real<lower=0> k;
	real<lower=min(x), upper=max(x)> x0;
	real<lower=0> sigma;
}

model {
	// Priors
	L ~ lognormal(1, 1);
	k ~ lognormal(1, 1);
	x0 ~ normal(200, 100);
	
	// Likelihood
	{
		vector[N] mu;
		mu = L ./ (1 + exp(-k*(x-x0)));
		log_y ~ normal(log(mu), sigma);
	}
}

/* UNCOMMENT THIS IF YOU WANT TO CALCULATE LOG-LIKELIHOOD
generated quantities {
	vector[N] log_lik;
	for (i in 1:N) log_lik[i] = normal_lpdf(log_y[i] | log(L / (1 + exp(-k*(x[i]-x0)))), sigma);
} 
*/
