# Load all the fit info into a big list and put each thing into its own data frame
# (fit info from the cmdstan fits)

fp <- '~/forestlight/stanoutput/fitinfo'

library(dplyr)
library(purrr)

mod_df <- expand.grid(dens_model = c('pareto', 'weibull'),
                      prod_model = c('power', 'exp'),
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = seq(1990, 2010, 5), 
                      stringsAsFactors = FALSE)

fit_info_list <- map(1:nrow(mod_df), function(i) {
  load(file.path(fp, paste0('info_', i, '.r')))
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

fit_ics <- map_dfr(fit_info_list, get_ics)
fit_ics <- cbind(mod_df[rep(1:nrow(mod_df), each=4),], fit_ics)
write.csv(fit_ics, file = '~/forestlight/ics_by_fg.csv', row.names = FALSE)

# Combine parameter credible intervals into single data frame.
# ------------------------------------------------------------

param_cis <- map_dfr(fit_info_list, 'param_cis')

write.csv(param_cis, file = '~/forestlight/paramci_by_fg.csv', row.names = FALSE)

# Combine predicted values into single data frame.
# ------------------------------------------------

pred_values <- map_dfr(fit_info_list, 'pred_interval')

write.csv(pred_values, file = '~/forestlight/ci_by_fg.csv', row.names = FALSE)

# Combine fitted slopes into single data frame.
# ---------------------------------------------

fitted_slopes <- map_dfr(fit_info_list, 'fitted_slopes')

write.csv(fitted_slopes, file = '~/forestlight/fitted_slopes_by_fg.csv', row.names = FALSE)

# Combine R-squared values into single data frame.
# ------------------------------------------------

r2s <- do.call(rbind, map(fit_info_list, 'r2s'))
r2s <- cbind(mod_df, r2s)

write.csv(r2s, file = '~/forestlight/r2_by_fg.csv', row.names = FALSE)
