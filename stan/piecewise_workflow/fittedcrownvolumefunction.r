fitted_totalvolume <- function(fit, dbh_pred, dens_form = NA, total_prod = NA, x_min = NULL, n_indiv = 1, ll = 1.1, ul = 316, pars_to_get, scaling_var = 'dbh') {
  require(purrr)
  require(pracma)
  

  
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975) 
   

	  pars_dens <- as.data.frame(do.call('cbind', extract(fit, pars_to_get)))

	  if (dens_form == '1') {
		dens_fitted <- sapply(dbh_pred, pdf_pareto, xmin = x_min, alpha = pars_dens[,'alpha']) 
	  }
	  if (dens_form == '2') {
		dens_fitted <- sapply(dbh_pred, pdf_2part, xmin = x_min, alpha_low = pars_dens[,'alpha_low'], alpha_high = pars_dens[,'alpha_high'], tau = pars_dens[,'tau'])
	  }
	  if (dens_form == '3') {
		dens_fitted <- sapply(dbh_pred, pdf_3part, xmin = x_min, alpha_low = pars_dens[,'alpha_low'], alpha_mid = pars_dens[,'alpha_mid'], alpha_high = pars_dens[,'alpha_high'], tau_low = pars_dens[,'tau_low'], tau_high = pars_dens[,'tau_high'])
	  }

      dens_fitted <- dens_fitted * n_indiv
	  dens_fitted_quant <- apply(dens_fitted, 2, quantile, probs = qprobs, na.rm = TRUE)
	  
	  # Use crown volume allometry to get crown volume for each dbh pred.
	  
	  volumes <- exp(-.681 + 2.02 * log(dbh_pred))

	  totalvol_fitted <- sweep(dens_fitted, 2, volumes, '*')
	  
	  
	  # Integrate fitted total production and multiply total production fitted values by the 
	  # ratio of total observed production and integral of fitted production
	  # (Use trapezoidal integration)
	  totalvol_fitted <- t(sapply(1:nrow(totalvol_fitted), function(i) {
		fitted_integral <- trapz(x = dbh_pred, y = totalvol_fitted[i,])
		totalvol_fitted[i,] * total_prod / fitted_integral
	  }))
      totalvol_fitted_quant <- apply(totalvol_fitted, 2, quantile, probs = qprobs, na.rm = TRUE)
 

	out <- data.frame(dbh = dbh_pred,
                    variable = 'total_crownvolume_fitted',
                    as.matrix(t(totalvol_fitted_quant)))

  
  setNames(out, c(scaling_var, 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))  
}
