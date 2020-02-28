data {
	int<lower=0> N;
	vector<lower=0>[N] x;
	vector<lower=0>[N] y;
}

parameters {
	// Pareto density
	real<lower=0> x_min;
	real<lower=0, upper=5> alpha;
}

model {
	// Prior: Pareto density
	x_min ~ lognormal(1, 1);
	alpha ~ lognormal(1, 1);
	// x_min uniform on its interval.
	
	// Likelihood: Pareto density
	x ~ pareto(x_min, alpha);
}

