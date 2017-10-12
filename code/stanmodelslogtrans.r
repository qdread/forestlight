# Stan models with log transformation
# Basic power law
# Power law plus exponential

model_powerlaw_logtrans <- '
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
	logx = log10(x);
	logy = log10(y);
	logx_pred = log10(x_pred);
  }
  
  parameters {
    real<lower=0> a;
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
'

# Power law times exponential
model_code_powerexp <- '
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
	logx = log10(x);
	logy = log10(y);
	logx_pred = log10(x_pred);
  }
  parameters {
    real<lower=0> a1;
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
      for (i in 1:N) mu[i] = (a + b * logx[i]) * (a * logx[i] ^ b + c);
      logy ~ normal(mu, sigma);
    }
  }
  generated quantities {
    vector[N_pred] y_pred;
    for (i in 1:N_pred) y_pred[i] = 10^((a + b * logx_pred[i]) * (a * logx_pred[i] ^ b + c));
  }
'
