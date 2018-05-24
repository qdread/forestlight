# Calculate Bayesian R2 of the production models, manually.

# Workflow
# -------------------------------------------------------------------------
# 1. Load the model fit
# 2. Load the data by sourcing the stan rdump for that model 
# 3. Extract the parameter estimates for each draw from the model fit
# 4. Plug the dbh (x) from the loaded data in to get the linear predictors for production
# 5. Get the residuals by subtracting the observed production (y) from predicted production
# 6. Calculate predicted variance / (predicted + residual variance) to get R2.

bayesian_rsquared_production <- function(dens_model, prod_model, fg, year) {
  require(rstan)
  require(purrr)

  # 1. Load CSVs with model fit as stanfit object
  fp <- '~/forestlight/stanoutput'
  files <- paste0('fit_', dens_model, 'x', prod_model, '_', fg, '_', year, '_', 1:3, '.csv')
  if (fg == 'alltree') files <- paste0('ss', files) # Use the 25K subset for all trees.
  fit <- read_stan_csv(file.path(fp, files))
  
  # 2. Load data
  fpdump <- '~/forestlight/stanrdump'
  dumpfile <- paste0('dump_', fg, '_', year, '.r')
  if (fg == 'alltree') dumpfile <- paste0('ss', dumpfile) # Use the 25K subset for all trees.
  source(file.path(fpdump, dumpfile)) # Creates variables x and y.
  
  # 3. Extract parameter estimates.
  pars_to_get <- c('beta0', 'beta1') 
  if (prod_model == 'exp') pars_to_get <- c(pars_to_get, 'a', 'b', 'c')
  
  pars <- extract(fit, pars_to_get)
  pars <- as.data.frame(do.call(cbind, pars))
  
  # 4. Plug in dbh (x) to get posterior estimates of linear predictor of production
  powerlaw_exp_log <- function(x, a, b, c, beta0, beta1) exp(-beta0) * x^beta1 * (-a * x ^ -b + c)
  powerlaw_log <- function(x, beta0, beta1) exp(-beta0) * x^beta1
  
  pmap(pars[,c('beta0','beta1')], powerlaw_log, x = x)
  
  # Take the log of the fitted values
  if (prod_model == 'power') {
    prod_fitted <- log(do.call(rbind, pmap(pars, powerlaw_log, x = x)))
  } else {
    prod_fitted <- log(do.call(rbind, pmap(pars, powerlaw_exp_log, x = x)))
  }

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

mod_df <- expand.grid(dens_model = c('pareto', 'weibull'),
                      prod_model = c('power', 'exp'),
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = seq(1990, 2010, 5), 
                      stringsAsFactors = FALSE)
					  
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
r2 <- bayesian_rsquared_production(mod_df$dens_model[task], mod_df$prod_model[task], mod_df$fg[task], mod_df$year[task])
save(r2, file = paste0('~/forestlight/stanoutput/fitinfo/r2_', task, '.r'))
			  
			  
###
# Combine output

r2list <- lapply(1:nrow(mod_df), function(i) {
	load(paste0('~/forestlight/stanoutput/fitinfo/r2_', i, '.r'))
	r2
})

r2df <- cbind(mod_df, do.call(rbind, r2list))
