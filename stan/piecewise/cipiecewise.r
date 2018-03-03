ci_piecewise <- function(fit, func, dbh_pred, x_min = NULL, n_indiv = 1) {
  pdf_paretoexp <- function(x, xmin, alpha, L, beta, Cx) {
    ifelse(x <= L, (alpha * xmin^alpha) / (x ^ (alpha+1)), Cx * exp(-beta * x))
  }
  pdf_weibullexp <- function(x, m, n, L, beta) {
    ifelse(x <= L, (m/n) * (x/n)^(m-1) * exp(-(x/n)^m), Cx * exp(-beta * x))
  }
  
  pars <- as.data.frame(do.call('cbind', extract(fit)))
  
  if (func == 'weibullexp') {
    dens_pred <- with(pars, sapply(dbh_pred, pdf_weibullexp, m = m, n = n, L = L, beta = beta, Cx = Cx))
  }
  if (func == 'paretoexp') {
    dens_pred <- with(pars, sapply(dbh_pred, pdf_paretoexp, xmin = x_min, alpha = alpha, L = L, beta = beta, Cx = Cx))
  }
  
  dens_pred <- dens_pred * n_indiv
  
  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  dens_pred_quant <- apply(dens_pred, 2, quantile, probs = qprobs)
  
  out <- data.frame(dbh = dbh_pred,
                    t(dens_pred_quant))
  
  setNames(out, c('dbh', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
}
