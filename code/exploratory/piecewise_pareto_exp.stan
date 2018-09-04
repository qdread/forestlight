functions {
	// Definition of PDF of the piecewise powerlaw (Pareto) and exponential distribution
	// Essentially a power law with positive slope below the cutoff value, and a power law with negative slope above it
	real paretoexp_lpdf(real x, real alpha, real lambda, real tau, real x_min) {
		real prob;
		real lprob;
		
		real C_con; // Continuity constant to make sure the two pieces meet
		real C_norm; // Normalization constant to make sure the pdf integrates to 1
		
		C_con = ( tau ^ (-1 - alpha) ) / ( lambda * exp(-lambda * tau) );
		C_norm = ((1/alpha) * (x_min ^ -alpha - tau ^ -alpha) + C_con * exp(-lambda * tau)) ^ -1;
		
		if (x < tau) prob = C_norm * ( x ^ (alpha - 1) );
		if (x >= tau) prob = C_con * C_norm * lambda * exp(-lambda * x);

		lprob = log(prob);

		return lprob;
	}
}

data {
	int<lower=0> N;
	vector<lower=0>[N] x;
    real<lower=0> x_min;
	real<lower=0> x_max;
}

parameters {
	real<lower=x_min, upper=x_max> tau;
	real<lower=0, upper=5> alpha;
	real<lower=0, upper=5> lambda;
}

model {
	// Prior
	// No prior set for tau (uniform on its interval)
	lambda ~ exponential(1) T[0, 5];					
	alpha ~ lognormal(1, 1) T[0, 5];
	// Likelihood
	for (i in 1:N) {
		target += paretoexp_lpdf(x[i] | alpha, lambda, tau, x_min);
	}
}

generated quantities {
	vector[N] log_lik; // Log-likelihood for getting info criteria later
	
	for (i in 1:N) {
		log_lik[i] = paretoexp_lpdf(x[i] | alpha, lambda, tau, x_min);
	}
}
