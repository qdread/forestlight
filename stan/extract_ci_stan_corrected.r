dens_prod_ci <- function(fit, dbh_pred, dens_form, prod_form, x_min = NULL, n_indiv = 1, delete_samples = NULL) {
  require(purrr)
  pdf_pareto <- function(x, xmin, alpha) (alpha * xmin^alpha) / (x ^ (alpha+1))
  cdf_pareto <- function(x, xmin, alpha) 1 - (xmin/x)^alpha
  powerlaw_exp_log <- function(x, a, b, c, beta0, beta1) exp(-beta0) * x^beta1 * (-a * x ^ -b + c)
  powerlaw_log <- function(x, beta0, beta1) exp(-beta0) * x^beta1
  powerlaw_bert_log <- function(x, beta0, beta1, beta2) exp(-beta0) * x^beta1 * (1 - exp(-beta2 * x))
  
  powerlaw_raw <- function(x, beta0, beta1) -beta0 + beta1*log(x)
  powerlaw_exp_raw <- function(x, a, b, c, beta0, beta1) -beta0 + beta1*log(x) + log(-a * x ^ -b + c)
  
  pars <- as.data.frame(do.call('cbind', extract(fit)))
  
  if (!is.null(delete_samples)) {
    pars <- pars[-delete_samples,]
  }
  
  if (dens_form == 'pareto') {
    dens_pred <- sapply(dbh_pred, pdf_pareto, xmin = x_min, alpha = pars[,'alpha'])
    dens_cdf <- sapply(dbh_pred, cdf_pareto, xmin = x_min, alpha = pars[,'alpha'])
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
    dens_cdf <- sapply(dbh_pred, pweibull, shape = pars[,'m'], scale = pars[,'n'])
    dens_cdf <- t(sapply(1:nrow(dens_cdf), function(i) {
      x <- (dens_cdf[i,] - trunc_pts[[i]][1])/diff(trunc_pts[[i]])
      x
    }))
    
  }
  if (prod_form == 'powerlaw') {
    prod_pred <- sapply(dbh_pred, powerlaw_log, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])
    prod_raw_pred <- sapply(dbh_pred, powerlaw_raw, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])
  }
  if (prod_form == 'powerlawexp') {
    prod_pred <- sapply(dbh_pred, powerlaw_exp_log, a = pars[,'a'], b = pars[,'b'], c = pars[,'c'], beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])
    prod_raw_pred <- sapply(dbh_pred, powerlaw_exp_raw, a = pars[,'a'], b = pars[,'b'], c = pars[,'c'], beta0 = pars[,'beta0'], beta1 = pars[,'beta1'])
  }
  if (prod_form == 'bertalanffy') {
    prod_pred <- sapply(dbh_pred, powerlaw_bert_log, beta0 = pars[,'beta0'], beta1 = pars[,'beta1'], beta2 = pars[,'beta2'])
  }
  
  # New way of predicting total production:
  # 1. Evaluate CDF at the dbh prediction points (do truncation if necessary). See above.
  # 2. Get predicted number of individuals in each bin using CDF
  # 3. Get predicted production at midpoint of bin
  # 4. Multiply the result of 2 and 3 together
  
  dbh_binwidth <- diff(dbh_pred)
  n_ind_bin <- t(apply(dens_cdf, 1, function(x) diff(x) * n_indiv)) # Step 2
  mids <- function(a) a[-length(a)] + diff(a)/2
  

  totalprod_pred <- n_ind_bin * t(apply(prod_pred, 1, mids))
  totalprod_pred <- t(sapply(1:nrow(totalprod_pred), function(i) totalprod_pred[i,] / dbh_binwidth))
  
  
  dens_pred <- dens_pred * n_indiv
  #totalprod_pred <- dens_pred * prod_pred
  
  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  dens_pred_quant <- apply(dens_pred, 2, quantile, probs = qprobs)
  prod_pred_quant <- apply(prod_pred, 2, quantile, probs = qprobs)
  totalprod_pred_quant <- apply(totalprod_pred, 2, quantile, probs = qprobs)
  
  out <- data.frame(dbh = c(dbh_pred, dbh_pred, mids(dbh_pred)),
                    variable = rep(c('density','production','total_production'), c(length(dbh_pred), length(dbh_pred), length(dbh_pred) - 1)),
                    rbind(t(dens_pred_quant), t(prod_pred_quant), t(totalprod_pred_quant)))

  
  # out <- data.frame(dbh = dbh_pred,
  #                   variable = rep(c('density','production','total_production'), each = length(dbh_pred)),
  #                   rbind(t(dens_pred_quant), t(prod_pred_quant), t(totalprod_pred_quant)))
  setNames(out, c('dbh', 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
}
