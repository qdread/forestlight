functions {
	// Definition of two part density function
	real logtriangular_lpdf(real x, real alpha_low, real alpha_high, real tau, real x_min) {
		real prob;
		real lprob;
		
		real C_con; // Continuity constant to make sure the two pieces meet
		real C_norm; // Normalization constant to make sure the pdf integrates to 1
		
		C_con = tau ^ -(alpha_high + alpha_low);
		C_norm = ( (C_con / alpha_low) * (tau ^ alpha_low - x_min ^ alpha_low) + ( tau ^ (-alpha_high) ) / alpha_high ) ^ -1;
		
		if (x < tau) prob = C_con * C_norm * ( x ^ (alpha_low - 1) );
		if (x >= tau) prob = C_norm * ( x ^ - (alpha_high + 1) );

		lprob = log(prob);

		return lprob;
	}
}

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
	// Two part density
	real<lower=x_min, upper=x_max> tau;
	real<lower=0, upper=10> alpha_low;
	real<lower=0, upper=5> alpha_high;
	
	// Loglinear production
	real<lower=0> beta0; // Intercept
	real<lower=0> beta1; // Slope
	real<lower=0> sigma;

}

model {
	// Prior: two part density
	// No prior set for tau (uniform on its interval)
	alpha_low ~ lognormal(1, 1) T[0, 10];	
	alpha_high ~ lognormal(1, 1) T[0, 5];
	
	// Priors: Loglinear production
	beta0 ~ lognormal(1, 1);
	beta1 ~ lognormal(1, 1);
	
	sigma ~ exponential(0.1);
	
	// Likelihood: two part density
	for (i in 1:N) {
		x[i] ~ logtriangular(alpha_low, alpha_high, tau, x_min);
	}
	
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
		log_lik_dens[i] = logtriangular_lpdf(x[i] | alpha_low, alpha_high, tau, x_min);
		log_lik_prod[i] = normal_lpdf(logy[i] | log(beta0) + beta1 * logx[i], sigma);
	}
}
