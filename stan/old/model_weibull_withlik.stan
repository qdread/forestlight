// Stan model 06 Feb
// Weibull density
// Power law production
// Edited 03 Mar: Added truncations on Weibull

data {
	int<lower=0> N;
	vector<lower=0>[N] x;
	real<lower=0> UL; // Lower truncation limit
	real<lower=0> LL; // Upper truncation limit
	
	int<lower=0> N_pred;
	vector<lower=0>[N_pred] x_pred;
}

parameters {
	// Weibull density
	real<lower=0> m;
	real<lower=0> n;
}

model {
	// Priors: Weibull density
	m ~ lognormal(1, 1);
	n ~ lognormal(1, 1);
	
	// Likelihood: Weibull density
	for (i in 1:N) x[i] ~ weibull(m, n) T[LL,UL];
}

generated quantities {
	// Log likelihood
	vector[N_pred] log_lik_dens;
	real k;
	
	k = -log_diff_exp(weibull_lcdf(UL | m, n), weibull_lcdf(LL | m, n));
	
	for (i in 1:N_pred) {
		log_lik_dens[i] = weibull_lpdf(x_pred[i] | m, n) + k;
	}
	
}
