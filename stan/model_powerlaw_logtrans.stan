// Power law for individual production fit.

data {
	int<lower=0> N;
	int<lower=0> N_pred;
	vector[N] x;
	vector[N] y;
	vector[N_pred] x_pred;
}

transformed data {
	vector[N] logx;
	vector[N] logy;
	vector[N_pred] logx_pred;
	logx = log(x)/log(10);
	logy = log(y)/log(10);
	logx_pred = log(x_pred)/log(10);
}

parameters {
	real a;
	real b;
	real<lower=0> sigma;
}

model {
	b ~ normal(2, 1);
	a  ~ normal(0, 10);
	sigma ~ exponential(0.01);
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = a + b * logx[i];
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	vector[N_pred] y_pred;
	for (i in 1:N_pred) y_pred[i] = 10^(a + b * logx_pred[i]);
}