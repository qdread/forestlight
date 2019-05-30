# Load all the fit info into a big list and put each thing into its own data frame
# (fit info from the cmdstan fits)
# Version created for piecewise. Only done 1995.

# Alternate version created for the ones where density and production are fit separately.

fp <- '~/forestlight/stanoutput/fitinfo'

library(dplyr)
library(purrr)

dens_df <- expand.grid(variable = 'density',
					   dens_model = 1:3,
					   prod_model = as.numeric(NA),
					   fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
					   year = 1995,
					   stringsAsFactors = FALSE)

prod_df <- expand.grid(variable = 'production',
					   dens_model = as.numeric(NA),
					   prod_model = 1:2,
					   fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
					   year = 1995,
					   stringsAsFactors = FALSE)
					   					   
totalprod_df <- expand.grid(variable = 'total_production',
					  dens_model = 1:3,
                      prod_model = 1:2,
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = 1995, 
                      stringsAsFactors = FALSE)

mod_df <- rbind(dens_df, prod_df, totalprod_df)
					  

fit_info_list <- map(1:nrow(mod_df), function(i) {
  load(file.path(fp, paste0('pw_info_'mod_df$variable[i],'_',i,'.r')))
  fit_info
})

# Extract WAIC and LOOIC for density and production from each one.
# ----------------------------------------------------------------
get_ics <- function(x) {
  ic_stats <- c('elpd_', 'p_', '', 'se_elpd_', 'se_p_', 'se_') # Names of stats to extract from each model fit
  wd <- unname(do.call(c, x$waic_dens[paste0(ic_stats, 'waic')]))
  wp <- unname(do.call(c, x$waic_prod[paste0(ic_stats, 'waic')]))
  ld <- unname(do.call(c, x$loo_dens[paste0(ic_stats, 'looic')]))
  lp <- unname(do.call(c, x$loo_prod[paste0(ic_stats, 'looic')]))
  out <- data.frame(variable = c('density', 'production'),
             criterion = c('WAIC','WAIC','LOOIC','LOOIC'),
             rbind(wd, wp, ld, lp))
  setNames(out, nm = c('variable','criterion','elpd','p','ic','se_elpd','se_p','se_ic'))
}

idx <- which(mod_df$variable %in% c('density', 'production'))
fit_ics <- map_dfr(fit_info_list[idx], get_ics)
fit_ics <- cbind(mod_df[idx,][rep(1:nrow(mod_df[idx,]), each=4),], fit_ics)
write.csv(fit_ics, file = '~/forestlight/newpiecewise_ics_by_fg.csv', row.names = FALSE)

# Combine parameter credible intervals into single data frame.
# ------------------------------------------------------------

param_cis <- map_dfr(fit_info_list[idx], 'param_cis')

write.csv(param_cis, file = '~/forestlight/newpiecewise_paramci_by_fg.csv', row.names = FALSE)

# Combine predicted values into single data frame.
# ------------------------------------------------

pred_values <- map_dfr(fit_info_list, 'pred_interval')

write.csv(pred_values, file = '~/forestlight/piecewise_ci_by_fg.csv', row.names = FALSE)

# Combine fitted slopes into single data frame.
# ---------------------------------------------

fitted_slopes <- map_dfr(fit_info_list, 'fitted_slopes')

write.csv(fitted_slopes, file = '~/forestlight/piecewise_fitted_slopes_by_fg.csv', row.names = FALSE)

# Combine R-squared values into single data frame.
# ------------------------------------------------

idx_p <- which(mod_df$variable %in% c('production'))
r2s <- do.call(rbind, map(fit_info_list[idx_p], 'r2s'))
r2s <- cbind(prod_df, r2s)

write.csv(r2s, file = '~/forestlight/newpiecewise_r2_by_fg.csv', row.names = FALSE)
