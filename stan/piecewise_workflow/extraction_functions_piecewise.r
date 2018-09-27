# Source code for all functions to extract information from piecewise model fits
# Part of the final stan workflow
# QDR / Forestlight / 10 Sep 2018

# Information to get:
# Parameter values and credible intervals
# Fitted values and credible intervals (and prediction intervals)
# Fitted log slopes and credible intervals
# Information criteria
# "Bayesian" R-squared

# 0. Define all density and production functions.
# -----------------------------------------------

# Density model 1 part
pdf_pareto <- function(x, xmin, alpha) (alpha * xmin^alpha) / (x ^ (alpha+1))

# Density model 2 part
pdf_2part <- function(x, xmin, alpha_low, alpha_high, tau) {
	C_con <- tau ^ -(alpha_high + alpha_low)
	C_norm <- ( (C_con / alpha_low) * (tau ^ alpha_low - xmin ^ alpha_low) + ( tau ^ (-alpha_high) ) / alpha_high ) ^ -1
	
	prob <- case_when(
		x < tau ~ C_con * C_norm * ( x ^ (alpha_low - 1) ),
		x >= tau ~ C_norm * ( x ^ - (alpha_high + 1) )
	)
	return(prob)
}

# Density model 3 part
pdf_3part <- function(x, xmin, alpha_low, alpha_mid, alpha_high, tau_low, tau_high) {
	C_con_low <- tau_low ^ -(alpha_mid + alpha_low)
	C_con_high <- tau_high ^ (alpha_high - alpha_mid)
	C_norm <- ( (C_con_low / alpha_low) * (tau_low ^ alpha_low - xmin ^ alpha_low) + (1 / alpha_mid) * (tau_low ^ -alpha_mid - tau_high ^ -alpha_mid) + (C_con_high / alpha_high) * (tau_high ^ -alpha_high) ) ^ -1
	
	prob <- case_when(
		x < tau_low ~ C_con_low * C_norm * ( x ^ (alpha_low - 1) ),
		x >= tau_low & x <= tau_high ~ C_norm * ( x ^ - (alpha_mid + 1) ),
		x > tau_high ~ C_con_high * C_norm * ( x ^ - (alpha_high + 1) )
	)
	return(prob)
}

# Production model 1 part
powerlaw_log <- function(x, beta0, beta1) beta0 * x^beta1

# Production model 2 parts (hinged)
powerlaw_hinge_log <- function(x, x0, beta0, beta1_low, beta1_high, delta) {
	xdiff <- log(x) - log(x0)
	exp( log(beta0) + beta1_low * xdiff + (beta1_high - beta1_low) * delta * log(1 + exp(xdiff / delta)) )
}

# 1. Function for parameter values and their credible intervals
# -------------------------------------------------------------

param_values <- function(fit, pars) {
   summ_fit <- summary(fit)
   param_cis <- summ_fit[[1]][pars, ]

   param_cis <- data.frame(parameter = dimnames(param_cis)[[1]], param_cis)
   names(param_cis)[5:9] <- c('q025', 'q25', 'q50', 'q75', 'q975')

   return(param_cis)
}

# 2. Function to get fitted values and credible/prediction intervals
# ------------------------------------------------------------------

fitted_predicted_values <- function(fit, dbh_pred, dens_form, prod_form, total_prod, x_min = NULL, n_indiv = 1, ll = 1.1, ul = 316, pars_to_get, delete_samples = NULL) {
  require(purrr)
  require(pracma)
    
  pars <- as.data.frame(do.call('cbind', extract(fit, c(pars_to_get, 'sigma'))))
  
  if (!is.null(delete_samples)) {
    pars <- pars[-delete_samples,]
  }
  
  if (dens_form == 1) {
    dens_fitted <- sapply(dbh_pred, pdf_pareto, xmin = x_min, alpha = pars[,'alpha']) 
  }
  if (dens_form == 2) {
    dens_fitted <- sapply(dbh_pred, pdf_2part, xmin = x_min, alpha_low = pars[,'alpha_low'], alpha_high = pars[,'alpha_high'], tau = pars[,'tau'])
  }
  if (dens_form == 3) {
    dens_fitted <- sapply(dbh_pred, pdf_3part, xmin = x_min, alpha_low = pars[,'alpha_low'], alpha_mid = pars[,'alpha_mid'], alpha_high = pars[,'alpha_high'], tau_low = pars[,'tau_low'], tau_high = pars[,'tau_high'])
  }
  if (prod_form == 1) {
    prod_fitted <- sapply(dbh_pred, powerlaw_log, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])
  }
  if (prod_form == 2) {
    prod_fitted <- sapply(dbh_pred, powerlaw_hinge_log, beta0 = pars[,'beta0'], beta1_low = pars[,'beta1_low'], beta1_high = pars[,'beta1_high'], x0 = pars[,'x0'], delta = pars[,'delta'])
  }
  
  dens_fitted <- dens_fitted * n_indiv
  totalprod_fitted <- dens_fitted * prod_fitted
  
  
  # Integrate fitted total production and multiply total production fitted values by the 
  # ratio of total observed production and integral of fitted production
  # (Use trapezoidal integration)
  totalprod_fitted <- t(sapply(1:nrow(totalprod_fitted), function(i) {
    fitted_integral <- trapz(x = dbh_pred, y = totalprod_fitted[i,])
    totalprod_fitted[i,] * total_prod / fitted_integral
  }))
  
  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  dens_fitted_quant <- apply(dens_fitted, 2, quantile, probs = qprobs, na.rm = TRUE)
  prod_fitted_quant <- apply(prod_fitted, 2, quantile, probs = qprobs, na.rm = TRUE)
  totalprod_fitted_quant <- apply(totalprod_fitted, 2, quantile, probs = qprobs, na.rm = TRUE)
  
  # Use high and low normal quantiles to get prediction intervals
  prediction_quantile <- function(p, fitted_values) {
	pred_raw <- apply(fitted_values, 2, function(mus) exp(qnorm(p, mean = log(mus), sd = pars[,'sigma'])))
	apply(pred_raw, 2, median, na.rm = TRUE)
  }
  prod_pred_quant <- map_dfc(qprobs, prediction_quantile, fitted_values = prod_fitted)
  totalprod_pred_quant <- map_dfc(qprobs, prediction_quantile, fitted_values = totalprod_fitted)
  
  out <- data.frame(dbh = dbh_pred,
                    variable = rep(c('density','production','total_production','production_fitted','total_production_fitted'), each = length(dbh_pred)),
                    rbind(t(dens_fitted_quant), as.matrix(prod_pred_quant), as.matrix(totalprod_pred_quant), t(prod_fitted_quant), t(totalprod_fitted_quant)))
  
  setNames(out, c('dbh', 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
    
}

# 3. Function to get the fitted log slope values at many points
# -------------------------------------------------------------

fitted_slope_ci <- function(fit, dbh_pred, dens_form, prod_form, total_prod, x_min = NULL, n_indiv = 1, ll = 1.1, ul = 316, pars_to_get, delete_samples = NULL) {
  require(purrr)
  require(pracma)
    
  pars <- as.data.frame(do.call('cbind', extract(fit, pars = pars_to_get)))
  
  if (!is.null(delete_samples)) {
    pars <- pars[-delete_samples,]
  }
  
  if (dens_form == 1) {
    dens_fitted <- sapply(dbh_pred, pdf_pareto, xmin = x_min, alpha = pars[,'alpha']) 
  }
  if (dens_form == 2) {
    dens_fitted <- sapply(dbh_pred, pdf_2part, xmin = x_min, alpha_low = pars[,'alpha_low'], alpha_high = pars[,'alpha_high'], tau = pars[,'tau'])
  }
  if (dens_form == 3) {
    dens_fitted <- sapply(dbh_pred, pdf_3part, xmin = x_min, alpha_low = pars[,'alpha_low'], alpha_mid = pars[,'alpha_mid'], alpha_high = pars[,'alpha_high'], tau_low = pars[,'tau_low'], tau_high = pars[,'tau_high'])
  }
  if (prod_form == 1) {
    prod_fitted <- sapply(dbh_pred, powerlaw_log, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])
  }
  if (prod_form == 2) {
    prod_fitted <- sapply(dbh_pred, powerlaw_hinge_log, beta0 = pars[,'beta0'], beta1_low = pars[,'beta1_low'], beta1_high = pars[,'beta1_high'], x0 = pars[,'x0'], delta = pars[,'delta'])
  }
  
  dens_fitted <- dens_fitted * n_indiv
  totalprod_fitted <- dens_fitted * prod_fitted
  
  
  # Integrate fitted total production and multiply total production fitted values by the 
  # ratio of total observed production and integral of fitted production
  # (Use trapezoidal integration)
  totalprod_fitted <- t(sapply(1:nrow(totalprod_fitted), function(i) {
    fitted_integral <- trapz(x = dbh_pred, y = totalprod_fitted[i,])
    totalprod_fitted[i,] * total_prod / fitted_integral
  }))
  
  # Here, get the fitted log slope for all the predicted values. (x is dbh_pred and y is each column of fitted value)
  # Uses formula x/y dy/dx for log slope.
  log_slope <- function(x, y) (x[-1]/y[-1]) * diff(y)/diff(x)
  dens_fittedslopes <- map_dfr(as.data.frame(t(dens_fitted)), ~ log_slope(dbh_pred, .))
  prod_fittedslopes <- map_dfr(as.data.frame(t(prod_fitted)), ~ log_slope(dbh_pred, .))
  totalprod_fittedslopes <- map_dfr(as.data.frame(t(totalprod_fitted)), ~ log_slope(dbh_pred, .))
    
  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  dens_fittedslopes_quant <- apply(dens_fittedslopes, 1, quantile, probs = qprobs, na.rm = TRUE)
  prod_fittedslopes_quant <- apply(prod_fittedslopes, 1, quantile, probs = qprobs, na.rm = TRUE)
  totalprod_fittedslopes_quant <- apply(totalprod_fittedslopes, 1, quantile, probs = qprobs, na.rm = TRUE)
  
  out <- data.frame(dbh = dbh_pred[-1],
                    variable = rep(c('density','production','total_production'), each = length(dbh_pred)-1),
                    rbind(t(dens_fittedslopes_quant), t(prod_fittedslopes_quant), t(totalprod_fittedslopes_quant)))
  
  setNames(out, c('dbh', 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
    
}

# 4. Function to get the Bayesian R-squared and its quantiles
# -----------------------------------------------------------

bayesian_rsquared_production <- function(fit, x, y, prod_model) {
  # 3. Extract parameter estimates.
  production_par <- list(c('beta0', 'beta1'),
						 c('beta0', 'beta1_low', 'beta1_high', 'x0', 'delta'))
  
  pars_to_get <- production_par[[prod_model]] 
  
  pars <- extract(fit, pars_to_get)
  pars <- as.data.frame(do.call(cbind, pars))
  
  # 4. Plug in dbh (x) to get posterior estimates of linear predictor of production
    
  # Take the log of the fitted values
  if (prod_model == 1) {
    prod_fitted <- log(do.call(rbind, pmap(pars, powerlaw_log, x = x)))
  } else {
    prod_fitted <- log(do.call(rbind, pmap(pars, powerlaw_hinge_log, x = x)))
  }

  # 5. Get residuals by subtracting log y from linear predictor
  resids <- -1 * sweep(prod_fitted, 2, log(y))
  
  # 6. Calculate variances and ratio
  pred_var <- apply(prod_fitted, 1, var)
  resid_var <- apply(resids, 1, var)
  r2s <- pred_var / (pred_var + resid_var)
  
  # Quantiles of rsq
  r2_quant <- quantile(r2s, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), na.rm = TRUE)
  setNames(r2_quant, c('q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
}

# 5. Master function that calls all the above functions to extract all relevant fit information
# ------------------------------------------------------------------------------------------

extract_all_fit <- function(dens_model, prod_model, fg, year, xmin, n, total_production, use_subset = FALSE) {
  require(rstan)
  require(loo)

  # Load CSVs as stanfit object
  print('Loading stan fit . . .')
  fp <- '~/forestlight/stanoutput'
  files <- paste0('fit_d', dens_model, 'p', prod_model, '_', fg, '_', year, '_', 1:3, '.csv')
  #if (use_subset) files <- paste0('ss', files) # Use the 25K subset if needed.
  fit <- read_stan_csv(file.path(fp,files))
  
  # Get credible intervals of parameters
  print('Calculating credible intervals of parameters . . .')
  density_par <- list(c('alpha'),
					  c('alpha_low', 'alpha_high', 'tau'),
					  c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high'))
  production_par <- list(c('beta0', 'beta1'),
						 c('beta0', 'beta1_low', 'beta1_high', 'x0', 'delta'))
  
  get_pars <- c(density_par[[dens_model]], production_par[[prod_model]])
  
  param_cis <- cbind(data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg), param_values(fit, get_pars))

  # Get fitted values with credible intervals and prediction intervals
  print('Calculating fitted values and prediction intervals . . .')
  pred_interval <- fitted_predicted_values(fit, dbh_pred, dens_form = dens_model, prod_form = prod_model, total_prod = total_production, x_min = xmin, n_indiv = n, pars_to_get = get_pars)
  pred_interval <- data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg, pred_interval)
  
  # Get fitted log slopes
  print('Calculating log slopes . . .')
  fitted_slopes <- fitted_slope_ci(fit, dbh_pred, dens_form = dens_model, prod_form = prod_model, total_prod = total_production, x_min = xmin, n_indiv = n)
  fitted_slopes <- data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg, fitted_slopes)

  # Get WAIC and LOOIC
  ll_dens <- extract_log_lik(fit, 'log_lik_dens')
  ll_prod <- extract_log_lik(fit, 'log_lik_prod')
  print('Calculating WAIC . . .')
  waic_dens <- waic(ll_dens)
  waic_prod <- waic(ll_prod)
  print('Calculating LOOIC . . .')
  loo_dens <- loo(ll_dens)
  loo_prod <- loo(ll_prod)
  
  # Calculate R-squared
  print('Calculating Bayesian R-squared . . .')
  # Load the dump file from the model so that the R2 can be calculated
  fpdump <- '~/forestlight/stanrdump'
  dumpfile <- paste0('dump_', fg, '_', year, '.r')

  source(file.path(fpdump, dumpfile)) # Creates variables x and y.
  
  r2s <- bayesian_rsquared_production(fit, x, y, prod_model)
  
  list(waic_dens = waic_dens, 
       waic_prod = waic_prod, 
       loo_dens = loo_dens, 
       loo_prod = loo_prod, 
       param_cis = param_cis, 
       pred_interval = pred_interval,
	   fitted_slopes = fitted_slopes,
	   r2s = r2s)
}
