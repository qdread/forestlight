# Load all the fit info into a big list and put each thing into its own data frame
# (fit info from the cmdstan fits)

fp <- '~/forestlight/stanoutput/fitinfo'

fit_info_list <- map(1:nrow(mod_df), function(i) {
  load(file.path(fp, paste0('info_', i, '.r')))
  fit_info
})

# Extract WAIC and LOOIC for density and production from each one.
# ----------------------------------------------------------------
get_ics <- function(x) {
  wd <- unname(do.call(c, x$waic_dens[1:6]))
  wp <- unname(do.call(c, x$waic_prod[1:6]))
  ld <- unname(do.call(c, x$loo_dens[1:6]))
  lp <- unname(do.call(c, x$loo_prod[1:6]))
  out <- data.frame(variable = c('density', 'production'),
             criterion = c('WAIC','WAIC','LOOIC','LOOIC'),
             rbind(wd, wp, ld, lp))
  setNames(out, nm = c('variable','criterion','elpd','p','ic','se_elpd','se_p','se_ic'))
}

fit_ics <- map(fit_info_list, get_ics)

fit_ics <- cbind(mod_df[rep(1:nrow(mod_df), each=4),], do.call(rbind, fit_ics))
write.csv(fit_ics, file = '~/forestlight/ics_by_fg.csv', row.names = FALSE)

# Combine parameter credible intervals into single data frame.
# ------------------------------------------------------------

param_cis <- do.call(rbind, map(fit_info_list, 'param_cis'))

write.csv(param_cis, file = '~/forestlight/paramci_by_fg.csv', row.names = FALSE)

# Combine predicted values into single data frame.
# ------------------------------------------------

pred_values <- do.call(rbind, map(fit_info_list, 'pred_interval'))

write.csv(pred_values, file = '~/forestlight/ci_by_fg.csv', row.names = FALSE)

