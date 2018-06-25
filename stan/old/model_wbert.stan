// Stan model 06 Feb
// Weibull density (truncated)
// Power law times Bertalanffy production

data {
	int<lower=0> N;
	vector<lower=0>[N] x;
	vector<lower=0>[N] y;
	real<lower=0> UL; // Lower truncation limit
	real<lower=0> LL; // Upper truncation limit
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
	// Power law times Bertalanffy production
	real<lower=0> beta0;
	real<lower=0> beta1;
	real<lower=0> beta2;
	real<lower=0> sigma;
}

model {
	// Priors: Weibull density
	m ~ lognormal(1, 1);
	n ~ lognormal(1, 1);
	// Priors: Power law times exponential production
	beta0 ~ normal(5, 2);
	beta1 ~ normal(0.5, 1);	
	beta2 ~ normal(0.5, 1);
	sigma ~ exponential(0.01);
	
	// Likelihood: Weibull density
	for (i in 1:N) x[i] ~ weibull(m, n) T[LL,UL];
	// Likelihood: Power law times exponential production
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = -beta0 + beta1 * logx[i] + log(1 - exp(-beta2 * x[i]));
	  logy ~ normal(mu, sigma);
	}
}
