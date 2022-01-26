# New mortality fit and plot
# QDR 25 Jan 2022

### FIT MODEL (remotely)

library(cmdstanr)
library(tidyverse)

mort <- read_csv('obs_mortalityindividuals.csv')
mort_mod <- cmdstan_model('mortreg_fg_v3.stan')


mort_data <- mort %>%
  filter(!fg %in% 'unclassified') %>%
  mutate(died = alive == 0) %>%
  select(fg, died, dbh)

mort_data_dump <- with(mort_data, list(N = nrow(mort_data), M = 5, x = dbh, y = as.numeric(died), fg = as.numeric(factor(fg))))

mort_fit <- mort_mod$sample(
  data = mort_data_dump,
  seed = 919,
  chains = 3,
  parallel_chains = 3,
  iter_warmup = 5000,
  iter_sampling = 1000,
  adapt_delta = 0.9,
  max_treedepth = 20
)

mort_fit$save_object(file = 'mort_fit.rds')


### EXTRACT FITTED VALUES (locally)

library(cmdstanr)
library(tidyverse)
library(forestscaling)

mort_fit <- readRDS('~/temp/forestlight/mort_fit.rds')
summ <- mort_fit$summary()

##### Get fitted values and their credible intervals
param_draws <- mort_fit$draws(variables = c('intercept', 'slope'), format = 'draws_df') %>%
  select(.draw, starts_with('intercept'), starts_with('slope')) %>%
  pivot_longer(-.draw, names_to = 'variable') %>%
  separate(variable, into = c('variable', 'fg')) %>%
  pivot_wider(id_cols = c(.draw, fg), names_from = variable)

dbh_pred <- seq(log10(1), log10(300), length.out = 101)
qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)

# Apply function to get fitted values
# Get quantiles in same pipe
logistic_fitted <- function(x, intercept, slope, ...) plogis(intercept + slope * x)

fitted_quant <- param_draws %>%
  group_by(fg, .draw) %>%
  group_modify(~ data.frame(dbh = 10^dbh_pred, y = logistic_fitted(x = dbh_pred, intercept = .x$intercept, slope = .x$slope))) %>%
  group_by(fg, dbh) %>%
  group_modify(~ as.data.frame(t(quantile(.x$y, probs = qprobs))) %>% 
                 setNames(c('q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975')))

write_csv(fitted_quant, '~/temp/forestlight/mortalitydbh_ci_by_fg.csv')
