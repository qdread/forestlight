// Logistic regression to model mortality as a function of light received per unit crown area
// Single model fit with random effect for FG (random slope, random intercept model)
// QDR 07 Oct 2019
			
data {
	int<lower=0> N;				// Number of trees
	int<lower=0> M;				// Number of FGs
	vector<lower=0>[N] x;		// Light per area
	int<lower=0,upper=1> y[N];	// Mortality 1990-1995
	int<lower=1,upper=M> fg[N];	// Mapping to functional groups 1-M
}

transformed data {
	vector[N] log10_x;
	log10_x = log10(x);
}

parameters {
	real alpha;							// Mean intercept
	real alpha_fg[M];					// Deviation of each FG from mean intercept
	real<lower=0,upper=10> sigma_alpha;	// Variation in intercepts
	real beta;							// Mean slope
	real beta_fg[M];					// Deviation of each FG from mean slope
	real<lower=0,upper=10> sigma_beta;	// Variation in slopes
	
}

model {
	// No priors assigned at the moment for alpha and beta (flat priors)
	alpha_fg ~ normal(0, sigma_alpha);
	beta_fg ~ normal(0, sigma_beta);
		
	// Likelihood
	for (n in 1:N) {
		y[n] ~ bernoulli(inv_logit((alpha + alpha_fg[fg[n]]) + (beta + beta_fg[fg[n]]) * log10_x));
	}
}
