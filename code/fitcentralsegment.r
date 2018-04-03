# Straight portion of power law.
# Use the bins


fp <- 'C:/Users/Q/google_drive/ForestLight/data/data_20mar2018' ## CHANGE PATH AS NEEDED

# Read all the csvs in directory.
for (i in dir(fp, pattern = '.csv')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}

library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

dat <- obs_totalprod %>%
  filter(year == 1995, fg == 'all')

linear_model <- '
data {
  int<lower=0> N;
  vector<lower=0>[N] x;
  vector<lower=0>[N] y;
}
parameters {
  real beta0;
  real beta1;
  real<lower=0> sigma;
}
model {
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  {
    vector[N] mu;
    for (i in 1:N) mu[i] = beta0 + beta1 * x[i];
    y ~ normal(mu, sigma);
  }
}
generated quantities {
  vector[N] ypred;
  for (i in 1:N) ypred[i] = (10^beta0) * (10^x[i])^beta1;
}
'

dat_linear <- with(dat, list(N = nrow(dat),
                             x = log10(bin_midpoint),
                             y = log10(bin_value)))

mod_linear <- stan_model(model_code = linear_model)

fit_linear <- sampling(mod_linear, data = dat_linear, chains = 3, iter = 6000, warmup = 5000, seed = 30303)
sum_linear <- summary(fit_linear)

sum_linear[[1]][c('beta0', 'beta1'), ]

dat_linear_straightportion <- with(dat[5:13,], list(N = length(5:13),
                             x = log10(bin_midpoint),
                             y = log10(bin_value)))
fit_linear_straight <- sampling(mod_linear, data = dat_linear_straightportion, chains = 3, iter = 6000, warmup = 5000, seed = 30303)
sum_linear_straight <- summary(fit_linear_straight)

sum_linear_straight[[1]][c('beta0', 'beta1'), ]

# Plot
library(ggplot2)
area_core <- 42.84

# Extract quantiles of fitted values
pred1 <- extract(fit_linear)$ypred
pred2 <- extract(fit_linear_straight)$ypred

library(purrr)

getpred <- function(x) {
x %>%
  as.data.frame %>%
  map(~ quantile(., probs = c(0.025, 0.5, 0.975))) %>%
  do.call(rbind, .) %>%
  as.data.frame %>%
  setNames(nm = c('q025', 'q50', 'q975'))
}

qpred1 <- data.frame(fit_to = 'all', dbh = dat$bin_midpoint, getpred(pred1))
qpred2 <- data.frame(fit_to = 'central section', dbh = dat$bin_midpoint[5:13], getpred(pred2))
qpred <- rbind(qpred1, qpred2)


ggplot(dat, aes(x = bin_midpoint)) + 
  geom_ribbon(data = qpred, aes(x = dbh, ymin = q025/area_core, ymax=q975/area_core, fill = fit_to), alpha = 0.2) +
  geom_line(data = qpred, aes(x = dbh, y = q50/area_core, color = fit_to)) +
  geom_point(aes(y = bin_value/area_core)) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  scale_x_log10(name = 'Diameter') +
  scale_y_log10(name = 'Biomass/year/hectare')
  
