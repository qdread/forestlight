// Power law times exponential fit for individual production

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
	real a1;
	real b1;
	real a;
	real b;
	real c;
	real<lower=0> sigma;
}

model {
	a ~ normal(1, 2);
	b ~ normal(0, 10);
	c ~ normal(1, 2);
	b1 ~ normal(2, 1);
	a1  ~ normal(0, 10);
	sigma ~ exponential(0.01);
	{
	  vector[N] mu;
	  for (i in 1:N) mu[i] = (a1 + b1 * logx[i]) * (a * logx[i] ^ b + c);
	  logy ~ normal(mu, sigma);
	}
}

generated quantities {
	vector[N_pred] y_pred;
	for (i in 1:N_pred) y_pred[i] = 10^((a1 + b1 * logx_pred[i]) * (a * logx_pred[i] ^ b + c));
}