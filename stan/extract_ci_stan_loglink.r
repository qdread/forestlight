# Function to get confidence intervals from fitted stan object
# Edit 29 Jan: Also functions to read csv files and functions to make stan plots
# Edit 30 Mar: Add the powerlaw and powerlaw * exponential without the log transform.

dens_prod_ci <- function(fit, dbh_pred, dens_form, prod_form, x_min = NULL, n_indiv = 1, delete_samples = NULL) {
  require(purrr)
  pdf_pareto <- function(x, xmin, alpha) (alpha * xmin^alpha) / (x ^ (alpha+1))
  #pdf_weibull <- function(x, m, n) (m/n) * (x/n)^(m-1) * exp(-(x/n)^m)
  powerlaw_exp_log <- function(x, a, b, c, beta0, beta1) beta0 * x^beta1 * (-a * x ^ -b + c)
  powerlaw_log <- function(x, beta0, beta1) beta0 * x^beta1
  powerlaw_bert_log <- function(x, beta0, beta1, beta2) beta0 * x^beta1 * (1 - exp(-beta2 * x))
  
  get_pars <- if (dens_form == 'pareto') 'alpha' else c('m', 'n')
  get_pars <- c(get_pars, c('beta0', 'beta1'))
  if (prod_form == 'powerlawexp') get_pars <- c(get_pars, c('a', 'b', 'c'))
  
  pars <- as.data.frame(do.call('cbind', extract(fit, pars = get_pars)))
  
  if (!is.null(delete_samples)) {
    pars <- pars[-delete_samples,]
  }
  
  if (dens_form == 'pareto') {
    dens_pred <- sapply(dbh_pred, pdf_pareto, xmin = x_min, alpha = pars[,'alpha'])
  }
  if (dens_form == 'weibull') {
    # Must manually rescale and remove upper and lower truncations
    ll <- 1.1
    ul <- 316
    trunc_pts <- pmap(pars, function(m, n, ...) pweibull(q = c(ll,ul), shape = m, scale = n))
    dens_pred <- sapply(dbh_pred, dweibull, shape = pars[,'m'], scale = pars[,'n']) 
    dens_pred <- t(sapply(1:nrow(dens_pred), function(i) {
      x <- dens_pred[i,]/diff(trunc_pts[[i]])
      x
      }))
  }
  if (prod_form == 'powerlaw') {
    prod_pred <- sapply(dbh_pred, powerlaw_log, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])
  }
  if (prod_form == 'powerlawexp') {
    prod_pred <- sapply(dbh_pred, powerlaw_exp_log, a = pars[,'a'], b = pars[,'b'], c = pars[,'c'], beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])
  }
  if (prod_form == 'bertalanffy') {
    prod_pred <- sapply(dbh_pred, powerlaw_bert_log, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'], beta2 = pars[,'beta2'])
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
