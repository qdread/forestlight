// Stan model 20 Feb
// Pareto piecewise to exponential density.
// Version 2. Add constant to make functions match at breakpoint

functions {
	real myexp_log (real x, real Cx, real beta) {
		return log(Cx) - beta * x;
	}
}

data {
	int<lower=0> N;
	vector<lower=0>[N] x;
}

transformed data {
	real<lower=0> min_x;
	real<lower=0> max_x;
	min_x = min(x);
	max_x = max(x);
}

parameters {
	// Pareto density
	real<lower=0, upper=5> alpha;
	real<lower=min_x> L; // breakpoint
	real<lower=0> beta;
}

transformed parameters {
	real Cx;
	Cx = ((alpha - min_x^alpha) / L^(alpha + 1)) / exp(-L * beta);
}

model {
	
	// Priors: Pareto to exponential density
	alpha ~ lognormal(1, 1) T[0, 5];
	L ~ uniform(min_x, max_x);
	beta ~ lognormal(1, 1);
		
	// Likelihood: Pareto to my customized exponential density
	for (i in 1:N) {
	  if (x[i] <= L) {
		  x[i] ~ pareto(min_x, alpha); 
	  } else {
		  x[i] ~ myexp(Cx, beta);
	  }
	  
	}

}
