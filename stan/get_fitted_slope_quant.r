# Script to get credible interval of fitted slopes for density, individual production, and total production
# Work directly with samples from MCMC
# QDR/Forestlight/17 May 2018

fitted_slope_ci <- function(fit, dbh_pred, dens_form, prod_form, total_prod, x_min = NULL, n_indiv = 1, ll = 1.1, ul = 316, delete_samples = NULL) {
  require(purrr)
  require(pracma)
  pdf_pareto <- function(x, xmin, alpha) (alpha * xmin^alpha) / (x ^ (alpha+1))
  powerlaw_exp_log <- function(x, a, b, c, beta0, beta1) exp(-beta0) * x^beta1 * (-a * x ^ -b + c)
  powerlaw_log <- function(x, beta0, beta1) exp(-beta0) * x^beta1
  powerlaw_exp_log_pred <- function(x, a, b, c, beta0, beta1, sigma) exp(rnorm(length(beta0), mean = -beta0 + beta1 * log(x) + log(-a * x ^ -b + c), sd = sigma))
  powerlaw_log_pred <- function(x, beta0, beta1, sigma) exp(rnorm(length(beta0), mean = -beta0 + beta1 * log(x), sd = sigma))
  
  pars <- as.data.frame(do.call('cbind', extract(fit)))
  
  if (!is.null(delete_samples)) {
    pars <- pars[-delete_samples,]
  }
  
  if (dens_form == 'pareto') {
    dens_pred <- sapply(dbh_pred, pdf_pareto, xmin = x_min, alpha = pars[,'alpha']) 
  }
  if (dens_form == 'weibull') {
    # Must manually rescale and remove upper and lower truncations
    trunc_pts <- pmap(pars, function(m, n, ...) pweibull(q = c(ll,ul), shape = m, scale = n))
    dens_pred <- sapply(dbh_pred, dweibull, shape = pars[,'m'], scale = pars[,'n']) 
    dens_pred <- t(sapply(1:nrow(dens_pred), function(i) {
      x <- dens_pred[i,]/diff(trunc_pts[[i]])
      x
      }))
  }
  if (prod_form == 'powerlaw') {
    prod_fitted <- sapply(dbh_pred, powerlaw_log, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])
  }
  if (prod_form == 'powerlawexp') {
    prod_fitted <- sapply(dbh_pred, powerlaw_exp_log, a = pars[,'a'], b = pars[,'b'], c = pars[,'c'], beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])

  }
  
  dens_pred <- dens_pred * n_indiv
  totalprod_fitted <- dens_pred * prod_fitted
  
  
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
  dens_fittedslopes <- map_dfr(as.data.frame(t(dens_pred)), ~ log_slope(dbh_pred, .))
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
