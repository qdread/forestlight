# Function to get confidence intervals from fitted stan object
# Edit 29 Jan: Also functions to read csv files and functions to make stan plots
# Edit 02 Apr: Multiply fitted total production by this ratio: total observed production / integral of fitted total production
# Edit 04 May: both fitted value credible interval and prediction credible interval

diag_plots <- function(fit) {
  require(bayesplot)
  # Trace plot
  parnames <- names(fit)[-length(names(fit))]
  print(mcmc_trace(as.array(fit), pars = parnames))
  print(stan_diag(fit, 'sample'))
  print(stan_diag(fit, 'stepsize'))

}

dens_prod_ci <- function(fit, dbh_pred, dens_form, prod_form, total_prod, x_min = NULL, n_indiv = 1, ll = 1.1, ul = 316, delete_samples = NULL) {
  require(purrr)
  require(pracma)
  pdf_pareto <- function(x, xmin, alpha) (alpha * xmin^alpha) / (x ^ (alpha+1))
  powerlaw_exp_log <- function(x, a, b, c, beta0, beta1) exp(-beta0) * x^beta1 * (-a * x ^ -b + c)
  powerlaw_log <- function(x, beta0, beta1) exp(-beta0) * x^beta1
  powerlaw_bert_log <- function(x, beta0, beta1, beta2) exp(-beta0) * x^beta1 * (1 - exp(-beta2 * x))
  powerlaw_exp_log_pred <- function(x, a, b, c, beta0, beta1, sigma) exp(rnorm(length(beta0), mean = -beta0 + beta1 * log(x) + log(-a * x ^ -b + c), sd = sigma))
  powerlaw_log_pred <- function(x, beta0, beta1, sigma) exp(rnorm(length(beta0), mean = -beta0 + beta1 * log(x), sd = sigma))
  powerlaw_bert_log_pred <- function(x, beta0, beta1, beta2, sigma) exp(-beta0) * x^beta1 * (1 - exp(-beta2 * x)) + rnorm(length(beta0), mean = 0, sd = exp(sigma))
  
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
	prod_pred <- sapply(dbh_pred, powerlaw_log_pred, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'], sigma = pars[,'sigma'])
  }
  if (prod_form == 'powerlawexp') {
    prod_fitted <- sapply(dbh_pred, powerlaw_exp_log, a = pars[,'a'], b = pars[,'b'], c = pars[,'c'], beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])
	prod_pred <- sapply(dbh_pred, powerlaw_exp_log_pred, a = pars[,'a'], b = pars[,'b'], c = pars[,'c'], beta0 = pars[,'beta0'], beta1 = pars[,'beta1'], sigma = pars[,'sigma'])

  }
  if (prod_form == 'bertalanffy') {
    prod_fitted <- sapply(dbh_pred, powerlaw_bert_log, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'], beta2 = pars[,'beta2'])
    prod_pred <- sapply(dbh_pred, powerlaw_bert_log_pred, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'], beta2 = pars[,'beta2'], sigma = pars[,'sigma'])
  }
  
  dens_pred <- dens_pred * n_indiv
  totalprod_pred <- dens_pred * prod_pred
  totalprod_fitted <- dens_pred * prod_fitted
  
  
  # Integrate fitted total production and multiply total production fitted values by the 
  # ratio of total observed production and integral of fitted production
  # (Use trapezoidal integration)
  totalprod_pred <- t(sapply(1:nrow(totalprod_pred), function(i) {
    fitted_integral <- trapz(x = dbh_pred, y = totalprod_pred[i,])
    totalprod_pred[i,] * total_prod / fitted_integral
  }))
  totalprod_fitted <- t(sapply(1:nrow(totalprod_fitted), function(i) {
    fitted_integral <- trapz(x = dbh_pred, y = totalprod_fitted[i,])
    totalprod_fitted[i,] * total_prod / fitted_integral
  }))
  
  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  dens_pred_quant <- apply(dens_pred, 2, quantile, probs = qprobs)
  prod_pred_quant <- apply(prod_pred, 2, quantile, probs = qprobs)
  totalprod_pred_quant <- apply(totalprod_pred, 2, quantile, probs = qprobs)
  prod_fitted_quant <- apply(prod_fitted, 2, quantile, probs = qprobs)
  totalprod_fitted_quant <- apply(totalprod_fitted, 2, quantile, probs = qprobs)
  
  out <- data.frame(dbh = dbh_pred,
                    variable = rep(c('density','production','total_production','production_fitted','total_production_fitted'), each = length(dbh_pred)),
                    rbind(t(dens_pred_quant), t(prod_pred_quant), t(totalprod_pred_quant), t(prod_fitted_quant), t(totalprod_fitted_quant)))
  
  setNames(out, c('dbh', 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
    
}
