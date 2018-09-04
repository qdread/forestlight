# Test stan model with small subset of the data

source('~/Dropbox/projects/forestlight/stanrdump_final/ssdump_alltree_1995.r') # Load dumped data

library(rstan)
library(bayesplot)
library(loo)
library(purrr)
library(dplyr)
library(lubridate)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)

pwmod <- stan_model('code/exploratory/piecewise_pareto_exp.stan')

n_iter <- 500
n_warm <- 400
n_chain <- 2

dummyinit <- function() list(alpha = 1, lambda = 1, tau = 50)

standata <- list(x = x, N = N, x_min = x_min, x_max = max(x))
random_seed <- 10101
pwfit <- sampling(pwmod, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = random_seed, pars = 'log_lik', include = FALSE, init = dummyinit)

pwsumm <- summary(pwfit)
pwsumm$summary
mcmc_trace(as.array(pwfit), pars = c('alpha','lambda','tau'))
