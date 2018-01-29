# Function to get confidence intervals from fitted stan object
# Edit 29 Jan: Also functions to read csv files and functions to make stan plots

diag_plots <- function(fit) {
  require(bayesplot)
  # Trace plot
  parnames <- names(fit)[-length(names(fit))]
  print(mcmc_trace(as.array(fit), pars = parnames))
  print(stan_diag(fit, 'sample'))
  print(stan_diag(fit, 'stepsize'))

}

dens_prod_ci <- function(fit, dbh_pred, dens_form, prod_form, x_min = NULL, n_indiv = 1) {
  pdf_pareto <- function(x, xmin, alpha) (alpha * xmin^alpha) / (x ^ (alpha+1))
  pdf_weibull <- function(x, shape, scale) (shape/scale) * (x/scale)^(shape-1) * exp(-(x/scale)^shape)
  powerlaw_exp_log <- function(x, a, b, c, a1, b1) 10^(a1 + b1 * log10(x)) * (a * log10(x) ^ b + c)
  powerlaw_log <- function(x, a, b) 10^(a + b * log10(x)) 
  
  pars <- do.call('cbind', extract(fit))

  if (dens_form == 'pareto') {
    dens_pred <- sapply(dbh_pred, pdf_pareto, xmin = x_min, alpha = pars[,'alpha'])
  }
  if (dens_form == 'weibull') {
    dens_pred <- sapply(dbh_pred, pdf_weibull, shape = pars[,'shape'], scale = pars[,'scale'])
  }
  if (prod_form == 'powerlaw') {
    prod_pred <- sapply(dbh_pred, powerlaw_log, a = pars[,'a'], b = pars[,'b'])
  }
  if (prod_form == 'powerlawexp') {
    prod_pred <- sapply(dbh_pred, powerlaw_exp_log, a = pars[,'a'], b = pars[,'b'], c = pars[,'c'], a1 = pars[,'a1'], b1 = pars[,'b1'])
  }
  
  dens_pred <- dens_pred * n_indiv
  totalprod_pred <- dens_pred * prod_pred
  
  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  dens_pred_quant <- apply(dens_pred, 2, quantile, probs = qprobs)
  prod_pred_quant <- apply(prod_pred, 2, quantile, probs = qprobs)
  totalprod_pred_quant <- apply(totalprod_pred, 2, quantile, probs = qprobs)
  
  out <- data.frame(dbh = dbh_pred,
                    variable = rep(c('density','production','total_production'), each = length(dbh_pred)),
                    rbind(t(dens_pred_quant), t(prod_pred_quant), t(totalprod_pred_quant)))
  
  setNames(out, c('dbh', 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
    
}
