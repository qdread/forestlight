# Load all the fit info into a big list and put each thing into its own data frame
# (fit info from the cmdstan fits)
# Version created for piecewise. Only done 1995.

# Alternate version created for the ones where density and production are fit separately.
# Edited 17 June 2019 to accommodate the ones with total light and crown volume used as scaling.

# DENSITY-PRODUCTION SCALINGS
# ===========================

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
  load(file.path(fp, paste0('pw_info_',mod_df$variable[i],'_',i,'.r')))
  fit_info
})

# Extract WAIC and LOOIC for density and production from each one.
# ----------------------------------------------------------------

get_ics <- function(x) {
	ics <- data.frame(criterion = c('WAIC', 'LOOIC'),
					  IC_value = c(x$waic['waic','Estimate'], x$loo['looic','Estimate']),
					  IC_stderr = c(x$waic['waic','SE'], x$loo['looic','SE']))
	return(ics)
}

idx <- which(mod_df$variable %in% c('density', 'production'))
fit_ics <- map_dfr(fit_info_list[idx], get_ics)
fit_ics <- cbind(mod_df[idx,][rep(1:nrow(mod_df[idx,]), each=2),], fit_ics)
write.csv(fit_ics, file = '~/forestlight/newpiecewise_ics_by_fg.csv', row.names = FALSE)

# Combine parameter credible intervals into single data frame.
# ------------------------------------------------------------

param_cis <- map_dfr(fit_info_list[idx], 'param_cis')

write.csv(param_cis, file = '~/forestlight/newpiecewise_paramci_by_fg.csv', row.names = FALSE)

# Combine predicted values into single data frame.
# ------------------------------------------------

pred_values <- map_dfr(fit_info_list, 'pred_interval')

write.csv(pred_values, file = '~/forestlight/newpiecewise_ci_by_fg.csv', row.names = FALSE)

# Combine fitted slopes into single data frame.
# ---------------------------------------------

fitted_slopes <- map_dfr(fit_info_list, 'fitted_slopes')

write.csv(fitted_slopes, file = '~/forestlight/newpiecewise_fitted_slopes_by_fg.csv', row.names = FALSE)

# Combine R-squared values into single data frame.
# ------------------------------------------------

idx_p <- which(mod_df$variable %in% c('production'))
r2s <- do.call(rbind, map(fit_info_list[idx_p], 'r2s'))
r2s <- cbind(prod_df, r2s)

write.csv(r2s, file = '~/forestlight/newpiecewise_r2_by_fg.csv', row.names = FALSE)


# DENSITY-TOTAL LIGHT SCALINGS
# ============================

mod_df <- rbind(prod_df, totalprod_df)
					  

fit_info_list <- map(1:nrow(mod_df), function(i) {
  load(file.path(fp, paste0('lightpw_info_',mod_df$variable[i],'_',i,'.r')))
  fit_info
})

# Extract WAIC and LOOIC for total light from each one.
# ----------------------------------------------------------------

get_ics <- function(x) {
	ics <- data.frame(criterion = c('WAIC', 'LOOIC'),
					  IC_value = c(x$waic['waic','Estimate'], x$loo['looic','Estimate']),
					  IC_stderr = c(x$waic['waic','SE'], x$loo['looic','SE']))
	return(ics)
}

idx <- which(mod_df$variable %in% c('production'))
fit_ics <- map_dfr(fit_info_list[idx], get_ics)
fit_ics <- cbind(mod_df[idx,][rep(1:nrow(mod_df[idx,]), each=2),], fit_ics)
fit_ics$variable <- 'incoming light individual'
write.csv(fit_ics, file = '~/forestlight/light_piecewise_ics_by_fg.csv', row.names = FALSE)

# Combine parameter credible intervals into single data frame.
# ------------------------------------------------------------

param_cis <- map_dfr(fit_info_list[idx], 'param_cis')
param_cis$variable <- 'incoming light individual'

write.csv(param_cis, file = '~/forestlight/light_piecewise_paramci_by_fg.csv', row.names = FALSE)

# Combine predicted values into single data frame.
# ------------------------------------------------

pred_values <- map_dfr(fit_info_list, 'pred_interval')
pred_values$variable <- factor(pred_values$variable, levels = c('production', 'production_fitted', 'total_production', 'total_production_fitted'), labels = c('incoming_light', 'incoming_light_fitted', 'total_incoming_light', 'total_incoming_light_fitted'))

write.csv(pred_values, file = '~/forestlight/light_piecewise_ci_by_fg.csv', row.names = FALSE)

# Combine fitted slopes into single data frame.
# ---------------------------------------------

fitted_slopes <- map_dfr(fit_info_list, 'fitted_slopes')
fitted_slopes$variable <- factor(fitted_slopes$variable, levels = c('production', 'total_production'), labels = c('incoming_light', 'total_incoming_light'))

write.csv(fitted_slopes, file = '~/forestlight/light_piecewise_fitted_slopes_by_fg.csv', row.names = FALSE)

# Combine R-squared values into single data frame.
# ------------------------------------------------

idx_p <- which(mod_df$variable %in% c('production'))
r2s <- do.call(rbind, map(fit_info_list[idx_p], 'r2s'))
r2s <- cbind(prod_df, r2s)
r2s$variable <- 'incoming light individual'

write.csv(r2s, file = '~/forestlight/light_piecewise_r2_by_fg.csv', row.names = FALSE)


# DENSITY x CROWN VOLUME SCALING
# FITTED VALUES ONLY
# ==============================

mod_df <- expand.grid(variable = 'total_production',
					  dens_model = 1:3,
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = 1995, 
                      stringsAsFactors = FALSE)

fit_info_list <- map(1:nrow(mod_df), function(i) {
  load(paste0('~/forestlight/stanoutput/fitinfo/volumepw_fittedvalues_',i,'.r'))
  fitted_totalvolume_values
})

pred_values <- imap_dfr(fit_info_list, ~ data.frame(mod_df[.y,-1], .x))
write.csv(pred_values, file = '~/forestlight/volume_piecewise_ci_by_fg.csv', row.names = FALSE)


# GROWTH AS DIAMETER PER TIME SCALINGS
# ADDED 24 JULY 2019
# ====================================

mod_df <- prod_df

fit_info_list <- map(1:nrow(mod_df), function(i) {
  load(file.path(fp, paste0('diamgrowthpw_info_',i,'.r')))
  fit_info
})

# Extract WAIC and LOOIC for total light from each one.
# ----------------------------------------------------------------

get_ics <- function(x) {
	ics <- data.frame(criterion = c('WAIC', 'LOOIC'),
					  IC_value = c(x$waic['waic','Estimate'], x$loo['looic','Estimate']),
					  IC_stderr = c(x$waic['waic','SE'], x$loo['looic','SE']))
	return(ics)
}

fit_ics <- map_dfr(fit_info_list, get_ics)
fit_ics <- cbind(mod_df[rep(1:nrow(mod_df), each=2),], fit_ics)
fit_ics$variable <- 'diameter growth individual'
write.csv(fit_ics, file = '~/forestlight/diamgrowth_piecewise_ics_by_fg.csv', row.names = FALSE)

# Combine parameter credible intervals into single data frame.
# ------------------------------------------------------------

param_cis <- map_dfr(fit_info_list, 'param_cis')
param_cis$variable <- 'diameter growth individual'

write.csv(param_cis, file = '~/forestlight/diamgrowth_piecewise_paramci_by_fg.csv', row.names = FALSE)

# Combine predicted values into single data frame.
# ------------------------------------------------

pred_values <- map_dfr(fit_info_list, 'pred_interval')
pred_values$variable <- factor(pred_values$variable, levels = c('production', 'production_fitted'), labels = c('diameter_growth', 'diameter_growth_fitted'))

write.csv(pred_values, file = '~/forestlight/diamgrowth_piecewise_ci_by_fg.csv', row.names = FALSE)

# Combine fitted slopes into single data frame.
# ---------------------------------------------

fitted_slopes <- map_dfr(fit_info_list, 'fitted_slopes')
fitted_slopes$variable <- 'diameter growth individual'

write.csv(fitted_slopes, file = '~/forestlight/diamgrowth_piecewise_fitted_slopes_by_fg.csv', row.names = FALSE)

# Combine R-squared values into single data frame.
# ------------------------------------------------

r2s <- do.call(rbind, map(fit_info_list, 'r2s'))
r2s <- cbind(prod_df, r2s)
r2s$variable <- 'diameter growth individual'

write.csv(r2s, file = '~/forestlight/diamgrowth_piecewise_r2_by_fg.csv', row.names = FALSE)
