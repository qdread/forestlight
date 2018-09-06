# Test stan model with small subset of the data

source('~/Dropbox/projects/forestlight/stanrdump_final/ssdump_alltree_1995.r') # Load dumped data

code_dir <- 'code/exploratory' # Local
code_dir <- '~/forestlight/testcode' # Remote


source('~/forestlight/stanrdump/ssdump_alltree_1995.r')

library(rstan)
library(bayesplot)
library(loo)
library(purrr)
library(dplyr)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)

# Make an even smaller subset of the data.
n_sub <- 2000
set.seed(222)
idx <- sample(length(x), n_sub, replace = F)
x <- x[idx]
y <- y[idx]
N <- n_sub

pwmod <- stan_model(file.path(code_dir, 'piecewise_pareto_exp.stan'))

n_iter <- 5000
n_warm <- 4000
n_chain <- 3

dummyinit <- function() list(alpha = 1, lambda = 1, tau = 50)

standata <- list(x = x, N = N, x_min = x_min, x_max = max(x))
random_seed <- 10101
pwfit <- sampling(pwmod, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = random_seed, pars = 'log_lik', include = FALSE, init = dummyinit)

pwsumm <- summary(pwfit)
pwsumm$summary
mcmc_trace(as.array(pwfit), pars = c('alpha','lambda','tau'))

# Compare with just Pareto

paretomod <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/model_h1.stan')
paretomod <- stan_model(file.path(code_dir,'model_h1.stan'))
paretofit <- sampling(paretomod, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = random_seed, pars = 'log_lik', include = FALSE)

summary(paretofit)$summary

# Other models
mod2 <- stan_model(file.path(code_dir,'model_h2.stan'))
mod3 <- stan_model(file.path(code_dir,'model_h3.stan'))

fit2 <- sampling(mod2, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = random_seed, pars = 'log_lik', include = FALSE)
fit3 <- sampling(mod3, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = 555, pars = 'log_lik', include = FALSE)

# Production model
prodmod <- stan_model(file.path(code_dir,'piecewise_production.stan'))
standataprod <- list(x = x, y = y, N = N, x_min = x_min, x_max = max(x))
prodfit <- sampling(prodmod, data = standataprod, chains = n_chain, iter = n_iter, warmup = n_warm, seed = 333, pars = 'log_lik_prod', include = FALSE)

prod_fun <- function(x, beta0, beta1_low, beta1_high, tau_p) exp(-beta0 + beta1_low * log(x) + beta1_high * (log(x) - log(tau_p)) * (x >= tau_p))

prod_fun1p <- function(x, beta0, beta1_low, beta1_high, tau_p) exp(-beta0 + beta1_low * log(x) )

x_pred <- exp(seq(log(min(x)), log(max(x)), length.out=101))
y_pred <- prod_fun(x_pred, 5.68, 2.36, 0.37, 165.9)

plot(x,y,log='xy', ylim=c(1e-3,1e3))
lines(x_pred,y_pred,col='red')


prodmod3part <- stan_model(file.path(code_dir,'piecewise_production3part.stan'))
prodfit3part <- sampling(prodmod3part, data = standataprod, chains = n_chain, iter = n_iter, warmup = n_warm, seed = 333, pars = 'log_lik_prod', include = FALSE)

prodsumm3part <- summary(prodfit3part)
prodsumm3part$summary[1:7,]

prodmodsq <- stan_model(file.path(code_dir,'productionsqrt.stan'))
prodfitsq <- sampling(prodmodsq, data = standataprod, chains = n_chain, iter = n_iter, warmup = n_warm, seed = 444, pars = 'log_lik_prod', include = FALSE)

summary(prodfitsq)
mcmc_trace(as.array(prodfitsq))

prodmodcurve <- stan_model(file.path(code_dir,'piecewise_curved_production.stan'))
prodfitcurve <- sampling(prodmodcurve, data = standataprod, chains = 2, iter = 1500, warmup = 1000, seed = 676, pars = c('log_lik_prod', 'x2'), include = FALSE, control = list(max_treedepth = 20, adapt_delta = 0.8))

summary(prodfitcurve)$summary
