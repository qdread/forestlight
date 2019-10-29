// Logistic regression to model mortality of each FG separately as a function of light received per unit crown area
// QDR 07 Oct 2019
			
data {
	int<lower=0> N;
	vector<lower=0>[N] x;		// Light per area
	int<lower=0,upper=1> y[N];	// Mortality 1990-1995
}

transformed data {
	vector[N] log10_x;
	log10_x = log10(x);
}

parameters {
	real alpha;
	real beta;
}

model {
	// No priors assigned at the moment (both parameters just get flat priors)
		
	// Likelihood
	y ~ bernoulli_logit(alpha + beta * log10_x);
}
