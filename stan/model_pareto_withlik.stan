// Stan model 06 Feb
// Pareto density only

data {
	int<lower=0> N;
	vector<lower=0>[N] x;
	real<lower=0> x_min;
	int<lower=0> N_pred;
	vector<lower=0>[N_pred] x_pred;
}

parameters {
	// Pareto density
	real<lower=0, upper=5> alpha;
}

model {
	// Priors: Pareto density
	alpha ~ lognormal(1, 1) T[0, 5];
	
	// Likelihood: Pareto density
	x ~ pareto(x_min, alpha);

}

generated quantities {
	// Log likelihood
	vector[N_pred] log_lik_dens;
	
	for (i in 1:N_pred) {
		log_lik_dens[i] = pareto_lpdf(x_pred[i] | x_min, alpha);
	}
	
}
