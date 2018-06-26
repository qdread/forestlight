# Process light fits

# Load fits and get CIs ---------------------------------------------------

fp <- '~/forestlight/stanoutput'

library(purrr)
library(tidyr)
library(dplyr)
library(rstan)

fgnames <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified', 'alltree')

mod_df <- expand.grid(fg=fgnames, year = c(1990, 1995), stringsAsFactors = FALSE)

# Define function to get credible intervals around the parameters and to get the interval around the fitted values

extract_light_ci <- function(fg, year) {
  getpar <- function(fit) as.data.frame(do.call(cbind, extract(fit)))
  fn_vonbert <- function(x, A, b, k, ...) A * (1 - b * exp(-k * x)) ^ 3
  light_pred <- exp(seq(log(1.1), log(412), length.out = 101))
  
  qprob <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  qname <- c('q025', 'q25', 'q50', 'q75', 'q975')
  
  filenames <- paste0('fit_light_', fg, '_', year, '_', 1:3, '.csv')
  fit <- read_stan_csv(file.path(fp, filenames))
  pars <- getpar(fit)
  
  summ_fit <- summary(fit)
  
  get_pars <- c('A','b','k', 'x_max', 'y_max', 'log_slope')
  
  param_cis <- cbind(data.frame(fg = fg, year = year), summ_fit[[1]][get_pars, ])
  
  param_cis <- cbind(parameter = dimnames(param_cis)[[1]], param_cis)
  names(param_cis)[7:11] <- c('q025', 'q25', 'q50', 'q75', 'q975')
  
  pred <- do.call(rbind, pmap(pars, fn_vonbert, x = light_pred))
  pred_quant <- pred %>%
    as.data.frame %>%
    map(~ quantile(., probs = qprob)) %>%
    do.call(rbind, .) %>%
    as.data.frame %>%
    setNames(nm = qname)
  pred_quant <- data.frame(fg = fg, year = year, light_area = light_pred, pred_quant)
  
  list(pars = param_cis, preds = pred_quant)
}

all_output <- pmap(mod_df, extract_light_ci)

# Combine the parameter confidence intervals and the predicted values into two single data frames.
all_pars <- do.call(rbind, map(all_output, 'pars'))
all_preds <- do.call(rbind, map(all_output, 'preds'))

write.csv(all_pars, file = '~/forestlight/lightbyarea_paramci_by_fg.csv', row.names = FALSE)
write.csv(all_preds, file = '~/forestlight/lightbyarea_predci_by_fg.csv', row.names = FALSE)

# Bayesian R2 for each fit ---------------------------------------------------

bayesian_rsquared_light <- function(fg, year) {
  require(rstan)
  require(purrr)
  
  # 1. Load CSVs with model fit as stanfit object
  fp <- '~/forestlight/stanoutput'
  files <- paste0('fit_light_', fg, '_', year, '_', 1:3, '.csv')
  fit <- read_stan_csv(file.path(fp, files))
  
  # 2. Load data
  fpdump <- '~/forestlight/stanrdump'
  dumpfile <- paste0('dump_light_', fg, '_', year, '.r')
  source(file.path(fpdump, dumpfile)) # Creates variables x and y.
  
  # 3. Extract parameter estimates.
  pars_to_get <- c('A', 'b', 'k') 
  pars <- extract(fit, pars_to_get)
  pars <- as.data.frame(do.call(cbind, pars))
  
  # 4. Plug in light (x) to get posterior estimates of linear predictor of production
  fn_vonbert <- function(x, A, b, k, ...) A * (1 - b * exp(-k * x)) ^ 3
  
  # Take the log of the fitted values
  prod_fitted <- log(do.call(rbind, pmap(pars, fn_vonbert, x = x)))
  
  # 5. Get residuals by subtracting log y from linear predictor
  resids <- -1 * sweep(prod_fitted, 2, log(y))
  
  # 6. Calculate variances and ratio
  pred_var <- apply(prod_fitted, 1, var)
  resid_var <- apply(resids, 1, var)
  r2s <- pred_var / (pred_var + resid_var)
  
  # Quantiles of rsq
  r2_quant <- quantile(r2s, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
  setNames(r2_quant, c('q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
}

fgnames <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified', 'alltree')

mod_df <- expand.grid(fg=fgnames, year = c(1990, 1995), stringsAsFactors = FALSE)

library(purrr)
r2s <- pmap(mod_df, bayesian_rsquared_light)
r2s <- cbind(mod_df, do.call(rbind, r2s))
write.csv(r2s, file = '~/forestlight/r2_light_by_fg.csv', row.names = FALSE)