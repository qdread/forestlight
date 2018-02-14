// Stan model 13 Feb
// Weibull density (manually specified)
// Power law production

functions {
  // Wrote this by taking the log of right hand side of equation 4b in Muller-Landau 2006
  // Must include the normalizing constant (same as for other Weibull)
  real myweib_log (real x, real m, real n) {
    return log(m / n) + (m - 1) * log(x / n) - (m * x / n);
  }
}

data {
	int<lower=0> N;
	vector<lower=0>[N] x;
	vector<lower=0>[N] y;
}

transformed data {
	vector[N] logx;
	vector[N] logy;
	logx = log(x);
	logy = log(y);
}

parameters {
	// Weibull density
	real<lower=0> m;
	real<lower=0> n;
	// Power law production
	real<lower=0> beta0;
	real<lower=0> beta1;
	real<lower=0> sigma;
}

model {
	// Priors: Weibull density
	m ~ lognormal(1, 1);
    n ~ lognormal(1, 1);
	// Priors: Power law production
	beta0 ~ normal(5, 2);
	beta1 ~ normal(0.5, 1);		
	sigma ~ exponential(0.01);
	
	// Likelihood: Weibull density
	for(i in 1:N) {
	  x[i] ~ myweib(m, n);  
	}
	// Likelihood: Power law production
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = -beta0 + beta1 * logx[i];
	  logy ~ normal(mu, sigma);
	}
}
