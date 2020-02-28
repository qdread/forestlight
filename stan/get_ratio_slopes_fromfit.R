get_ratio_slopes_fromfit <- function(fit, fg_top, fg_bottom, year, xmin, n, use_subset = FALSE, n_chains = 3, scaling_var = 'dbh', fp = '~/forestlight/stanoutput', densityfitprefix = 'fit_density', productionfitprefix = 'fit_production', scalingtype = 'production') {
  # fit must be a named list with names density_top, production_top, density_bottom, production_bottom
  require(rstan)
  require(Brobdingnag)
  require(purrr)
  
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975) 
  # Uses formula x/y dy/dx for log slope.
  log_slope <- function(x, y) (x[-1]/y[-1]) * diff(y)/diff(x)  
  
  # names of parameters
  density_par <- list('1' = c('alpha'),
                      '2' = c('alpha_low', 'alpha_high', 'tau'),
                      '3' = c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high'))
  production_par <- list('1' = c('beta0', 'beta1', 'sigma'),
                         '2' = c('beta0', 'beta1_low', 'beta1_high', 'x0', 'delta', 'sigma'))					  
  
  pars_to_get <- list(density = density_par[[dens_model]], production = production_par[[prod_model]])
  
  message('Calculating density slopes . . .')
  # extract density parameters
  pars_dens_top <- as.data.frame(do.call('cbind', extract(fit[['density_top']], pars_to_get[['density']])))
  pars_dens_bottom <- as.data.frame(do.call('cbind', extract(fit[['density_bottom']], pars_to_get[['density']])))
  
  # get fitted density values
  if (dens_model == '1') {
    dens_fitted_top <- sapply(dbh_pred, pdf_pareto, xmin = xmin, alpha = pars_dens_top[,'alpha']) 
    dens_fitted_bottom <- sapply(dbh_pred, pdf_pareto, xmin = xmin, alpha = pars_dens_bottom[,'alpha']) 
  }
  if (dens_model == '2') {
    dens_fitted_top <- sapply(dbh_pred, pdf_2part, xmin = xmin, alpha_low = pars_dens_top[,'alpha_low'], alpha_high = pars_dens_top[,'alpha_high'], tau = pars_dens_top[,'tau'])
    dens_fitted_bottom <- sapply(dbh_pred, pdf_2part, xmin = xmin, alpha_low = pars_dens_bottom[,'alpha_low'], alpha_high = pars_dens_bottom[,'alpha_high'], tau = pars_dens_bottom[,'tau'])
  }
  if (dens_model == '3') {
    dens_fitted_top <- sapply(dbh_pred, pdf_3part, xmin = xmin, alpha_low = pars_dens_top[,'alpha_low'], alpha_mid = pars_dens_top[,'alpha_mid'], alpha_high = pars_dens_top[,'alpha_high'], tau_low = pars_dens_top[,'tau_low'], tau_high = pars_dens_top[,'tau_high'])
    dens_fitted_bottom <- sapply(dbh_pred, pdf_3part, xmin = xmin, alpha_low = pars_dens_bottom[,'alpha_low'], alpha_mid = pars_dens_bottom[,'alpha_mid'], alpha_high = pars_dens_bottom[,'alpha_high'], tau_low = pars_dens_bottom[,'tau_low'], tau_high = pars_dens_bottom[,'tau_high'])
  }
  
  # calculate fitted slopes
  dens_fitted_top <- dens_fitted_top * n[1]
  dens_fittedslopes_top <- map_dfr(as.data.frame(t(dens_fitted_top)), ~ log_slope(dbh_pred, .))
  
  dens_fitted_bottom <- dens_fitted_bottom * n[2]
  dens_fittedslopes_bottom <- map_dfr(as.data.frame(t(dens_fitted_bottom)), ~ log_slope(dbh_pred, .))
  
  message('Getting production values . . .')
  # extract production parameters
  pars_prod_top <- as.data.frame(do.call('cbind', extract(fit[['production_top']], pars_to_get[['production']])))
  pars_prod_bottom <- as.data.frame(do.call('cbind', extract(fit[['production_bottom']], pars_to_get[['production']])))
  
  # get fitted production values
  if (prod_model == '1') {
    prod_fitted_top <- sapply(dbh_pred, powerlaw_log, beta0 = pars_prod_top[,'beta0'], beta1 = pars_prod_top[,'beta1'])
    prod_fitted_bottom <- sapply(dbh_pred, powerlaw_log, beta0 = pars_prod_bottom[,'beta0'], beta1 = pars_prod_bottom[,'beta1'])
  }
  if (prod_model == '2') {
    prod_fitted_top <- sapply(dbh_pred, powerlaw_hinge_log, beta0 = pars_prod_top[,'beta0'], beta1_low = pars_prod_top[,'beta1_low'], beta1_high = pars_prod_top[,'beta1_high'], x0 = pars_prod_top[,'x0'], delta = pars_prod_top[,'delta'])
    prod_fitted_bottom <- sapply(dbh_pred, powerlaw_hinge_log, beta0 = pars_prod_bottom[,'beta0'], beta1_low = pars_prod_bottom[,'beta1_low'], beta1_high = pars_prod_bottom[,'beta1_high'], x0 = pars_prod_bottom[,'x0'], delta = pars_prod_bottom[,'delta'])
  }
  
  message('Calculating total production slopes  . . .')
  
  totalprod_fitted_top <- dens_fitted_top * prod_fitted_top
  totalprod_fittedslopes_top <- map_dfr(as.data.frame(t(totalprod_fitted_top)), ~ log_slope(dbh_pred, .))
  
  totalprod_fitted_bottom <- dens_fitted_bottom * prod_fitted_bottom
  totalprod_fittedslopes_bottom <- map_dfr(as.data.frame(t(totalprod_fitted_bottom)), ~ log_slope(dbh_pred, .))
  
  message('Getting quantiles of ratio slopes (almost done!) . . .')
  # Ratio of density slopes
  ratio_fittedslopes_dens <- dens_fittedslopes_top - dens_fittedslopes_bottom
  # Ratio of total production slopes
  ratio_fittedslopes_totalprod <- totalprod_fittedslopes_top - totalprod_fittedslopes_bottom
  
  # Quantiles for each ratio
  ratio_fittedslopes_dens_quant <- apply(ratio_fittedslopes_dens, 1, quantile, probs = qprobs, na.rm = TRUE)
  ratio_fittedslopes_totalprod_quant <- apply(ratio_fittedslopes_totalprod, 1, quantile, probs = qprobs, na.rm = TRUE)
  
  # Also get quantiles for the fitted values themselves so the lines can be plotted
  dens_fitted_top_quant <- apply(dens_fitted_top, 2, quantile, probs = qprobs, na.rm = TRUE)
  dens_fitted_bottom_quant <- apply(dens_fitted_bottom, 2, quantile, probs = qprobs, na.rm = TRUE)
  
  # Process output into neat data frame and return
  out <- data.frame(dbh = dbh_pred[-1],
                    variable = rep(c('density', 'total production'), each = length(dbh_pred) - 1),
                    rbind(t(ratio_fittedslopes_dens_quant), t(ratio_fittedslopes_totalprod_quant)))
  
  out <- setNames(out, c(scaling_var, 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
}
