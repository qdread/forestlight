# Bayesian credible intervals on individual production trends.
# Power law, and power law times exponential
# Add predicted values too.

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_04oct'
load(file.path(fpdata, 'rawdataobj.r'))


# Basic power law.
# https://stats.stackexchange.com/questions/199910/fitting-power-function-to-data
model_code_powerlaw <- '
  data {
    int<lower=0> N;
    int<lower=0> N_pred;
    vector[N] x;
    vector[N] y;
    vector[N_pred] x_pred;
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
      for (i in 1:N) mu[i] = a * x[i] ^ b;
      y ~ lognormal(mu, sigma);
    }
  }
  generated quantities {
    vector[N_pred] y_pred;
    for (i in 1:N_pred) y_pred[i] = a * x_pred[i] ^ b;
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
      for (i in 1:N) mu[i] = (a1 * x[i] ^ b1) ^ (a * x[i] ^ b + c);
      y ~ lognormal(mu, sigma);
    }
  }
  generated quantities {
    vector[N_pred] y_pred;
    for (i in 1:N_pred) y_pred[i] = (a1 * x_pred[i] ^ b1) ^ (a * x_pred[i] ^ b + c);
  }
'

estimate.model <- function(x, y, stanmodel) {
  N <- length(x)
  data <- list(N = N, x = x, y = y)
  fit <- sampling(stanmodel, data = data, chains = 3, iter = 2000, warmup = 1000)
  return(fit)
}

estimate_model_pred <- function(x, y, x_pred, stanmodel) {
  N <- length(x)
  N_pred <- length(x_pred)
  data <- list(N = N, N_pred = N_pred, x = x, y = y, x_pred = x_pred)
  fit <- sampling(stanmodel, data = data, chains = 3, iter = 2000, warmup = 1000)
  return(fit)
}

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

powerlaw_model <- stan_model(model_code = model_code_powerlaw)

gap_test_fit_powerlaw <- with(gapdat[[6]], estimate.model(x = dbh_corr, y = production, stanmodel = powerlaw_model))
shade_test_fit_powerlaw <- with(shadedat[[6]], estimate.model(x = dbh_corr, y = production, stanmodel = powerlaw_model))

gap_test_fit_powerlaw_pred <- with(gapdat[[6]], estimate_model_pred(x = dbh_corr, y = production, x_pred = seq(min(dbh_corr),max(dbh_corr),length.out = 10), stanmodel = powerlaw_model))

# Plot power law with confidence interval for gap trees 2010.
gap_summary <- summary(gap_test_fit_powerlaw_pred)


gap_sum_dat <- data.frame(dbh_corr = with(gapdat[[6]], seq(min(dbh_corr),max(dbh_corr),length.out = 10)),
                          pred_median = gap_summary$summary[4:13,'50%'],
                          pred_025 = gap_summary$summary[4:13,'2.5%'],
                          pred_975 = gap_summary$summary[4:13,'97.5%'])

library(cowplot)

ggplot(gap_sum_dat, aes(x=dbh_corr, y=pred_median, ymin=pred_025, ymax=pred_975)) +
  geom_ribbon(alpha = 0.5) +
  geom_line(color = 'red')+ 
  scale_x_log10() + scale_y_log10()

#################

powerlawexp_model <- stan_model(model_code = model_code_powerexp)

gap_test_fit_powerlawexp_pred <- with(gapdat[[6]], estimate_model_pred(x = dbh_corr, y = production, x_pred = seq(min(dbh_corr),max(dbh_corr),length.out = 10), stanmodel = powerlawexp_model))

gap_summary_exp <- summary(gap_test_fit_powerlawexp_pred)


gap_exp_sum_dat <- data.frame(dbh_corr = with(gapdat[[6]], seq(min(dbh_corr),max(dbh_corr),length.out = 10)),
                          pred_median = gap_summary_exp$summary[7:16,'50%'],
                          pred_025 = gap_summary_exp$summary[7:16,'2.5%'],
                          pred_975 = gap_summary_exp$summary[7:16,'97.5%'])

gap_all_sum_dat <- rbind(cbind(model = 'powerlaw', gap_sum_dat),
                         cbind(model = 'powerlaw_exp', gap_exp_sum_dat))

ggplot(gap_all_sum_dat, aes(x=dbh_corr, y=pred_median, ymin=pred_025, ymax=pred_975, group = model, color = model, fill = model)) +
  geom_ribbon(alpha = 0.5) +
  geom_line() + 
  scale_x_log10() + scale_y_log10()

##################
# With transformed data.

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
'

# Power law times exponential
model_powerlawexp_logtrans <- '
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
  for (i in 1:N) mu[i] = (a + b * logx[i]) * (a * logx[i] ^ b + c);
  logy ~ normal(mu, sigma);
}
}
generated quantities {
vector[N_pred] y_pred;
for (i in 1:N_pred) y_pred[i] = 10^((a + b * logx_pred[i]) * (a * logx_pred[i] ^ b + c));
}
'

powerlaw_trans_model <- stan_model(model_code = model_powerlaw_logtrans)
powerlawexp_trans_model <- stan_model(model_code = model_powerlawexp_logtrans)

gap_fit_powerlaw_trans <- with(gapdat[[6]], estimate_model_pred(x = dbh_corr, y = production, x_pred = seq(min(dbh_corr),max(dbh_corr),length.out = 10), stanmodel = powerlaw_trans_model))
gap_fit_powerlawexp_trans <- with(gapdat[[6]], estimate_model_pred(x = dbh_corr, y = production, x_pred = seq(min(dbh_corr),max(dbh_corr),length.out = 50), stanmodel = powerlawexp_trans_model))

gap_plaw_trans_summary <- summary(gap_fit_powerlaw_trans)
gap_exp_trans_summary <- summary(gap_fit_powerlawexp_trans)

# Plot.
yrows <- grep('y_pred', dimnames(gap_plaw_trans_summary$summary)[[1]])
gap_plaw_trans_sum_dat <- data.frame(dbh_corr = with(gapdat[[6]], seq(min(dbh_corr),max(dbh_corr),length.out = 10)),
                          pred_median = gap_plaw_trans_summary$summary[yrows,'50%'],
                          pred_025 = gap_plaw_trans_summary$summary[yrows,'2.5%'],
                          pred_975 = gap_plaw_trans_summary$summary[yrows,'97.5%'])
yrows <- grep('y_pred', dimnames(gap_exp_trans_summary$summary)[[1]])
gap_exp_trans_sum_dat <- data.frame(dbh_corr = with(gapdat[[6]], seq(min(dbh_corr),max(dbh_corr),length.out = 50)),
                              pred_median = gap_exp_trans_summary$summary[yrows,'50%'],
                              pred_025 = gap_exp_trans_summary$summary[yrows,'2.5%'],
                              pred_975 = gap_exp_trans_summary$summary[yrows,'97.5%'])

gap_all_trans_sum_dat <- rbind(cbind(model = 'powerlaw', gap_plaw_trans_sum_dat),
                         cbind(model = 'powerlaw_exp', gap_exp_trans_sum_dat))

ggplot(gap_all_trans_sum_dat) +
  geom_hex(aes(x=dbh_corr, y=production), data=gapdat[[6]]) +
  scale_fill_gradient(low = 'gray95', high='gray10') +
  scale_color_manual(values = c('indianred','dodgerblue')) +
  geom_line(aes(x=dbh_corr, y=pred_median, group = model, color = model), size=1) + 
  geom_line(aes(x=dbh_corr, y=pred_025, group = model, color = model), linetype = 'dotted', size=0.75) +
  geom_line(aes(x=dbh_corr, y=pred_975, group = model, color = model), linetype = 'dotted', size=0.75) +
  scale_x_log10() + scale_y_log10()

ggsave('C:/Users/Q/google_drive/ForestLight/figs/new_cutoff_plots/gap_indiv_production_predictioninterval_2010.png', height=6, width=7, dpi=400)

ggplot(gap_all_trans_sum_dat) +
  scale_color_manual(values = c('indianred','dodgerblue')) +
  scale_fill_manual(values = c('indianred','dodgerblue')) +
  geom_ribbon(aes(x=dbh_corr, ymin=pred_025, ymax=pred_975, group = model, fill = model), alpha=0.25) + 
  geom_line(aes(x=dbh_corr, y=pred_median, group = model, color = model)) + 
  scale_x_log10() + scale_y_log10()

ggsave('C:/Users/Q/google_drive/ForestLight/figs/new_cutoff_plots/gap_indiv_production_predictioninterval_2010nodata.png', height=6, width=7, dpi=400)
