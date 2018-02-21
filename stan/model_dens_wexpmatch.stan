// Stan model 20 Feb
// Weibull piecewise to exponential density.
// Version 2. Add constant to make the functions match at breakpoint.

functions {
  // Wrote this by taking the log of right hand side of equation 4b in Muller-Landau 2006
  // Must include the normalizing constant (same as for other Weibull)
  real myweib_log (real x, real m, real n) {
    return log(m / n) + (m - 1) * log(x / n) - (m * x / n);
  }
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
	// Weibull density
	real<lower=0> m;
	real<lower=0> n;
	real<lower=min_x> L; // breakpoint
	real<lower=0> beta;
}

transformed parameters {
	real Cx;
	Cx = ((m/n) * (L/n)^(m-1) * exp(-(L/n)^m)) / exp(-L * beta);
}

model {
	
	// Priors: Weibull to exponential density
	m ~ lognormal(1, 1);
    n ~ lognormal(1, 1);
	L ~ uniform(min_x, max_x);
	beta ~ lognormal(1, 1);
		
	// Likelihood: Weibull to exponential density
	for (i in 1:N) {
	  if (x[i] <= L) {
		  x[i] ~ myweib(m, n); 
	  } else {
		  x[i] ~ myexp(Cx, beta);
	  }
	  
	}

}
