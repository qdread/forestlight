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

generated quantities {
	//real<lower=0> max_slope;
	//real x_max;
	//real y_max;
	real log_slope;
	
	/* max_slope = (4.0/9.0) * k * G;
	x_max = (-1 / k) * log(1/(3*b1));
	y_max = G * (1 - b1 * exp(-k * x_max))^3;
	log_slope = max_slope * x_max / y_max; */
	
	// edit April 26: redo log slope
	{
		vector[N] all_slopes;
		vector[N] y_fit;
		for (i in 1:N) {
			y_fit[i] = G * (1 - b1 * exp(-k * x[i])) ^ 3;
			all_slopes[i] = 3*b1*k*G*((1-b1*exp(-k*x[i]))^2)*exp(-k*x[i])*(x[i]/y_fit[i]);
		}
		log_slope = max(all_slopes);
	}
	
	
	// Uncomment the lines below if you want to output the log-likelihood.
	// vector[N] log_lik;
	// for (i in 1:N) log_lik[i] = normal_lpdf(log_y[i] | log(G * (1 - b1 * exp(-k * x[i])) ^ 3;), sigma);
} 
