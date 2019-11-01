# For each CMDstan model fit, get the following
# Credible intervals of parameters
# Credible and prediction intervals to make graphs
# Information criteria
# Fitted slopes in log space
# Bayesian R-squared


# Alternate version created 30 May 2019 for fits where density and production are done separately.
# Modified 17 June 2019 to do separate extractions for the ones where (a) crown volume and (b) incoming light are used as the y variables instead of production.

# DENSITY - PRODUCTION SCALINGS
# =============================

source('~/forestlight/stancode/model_output_extraction_functions.r')

library(purrr)
library(dplyr)
library(foreach)
library(doParallel)

dens_df <- expand.grid(variable = 'density',
					   dens_model = 1:2,
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
					   					   
mod_df <- expand.grid(variable = 'total_production',
					  dens_model = 1:2,
                      prod_model = 1:2,
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = 1995, 
                      stringsAsFactors = FALSE)

mod_df <- rbind(dens_df, prod_df, mod_df)
					  
min_n <- read.csv('~/forestlight/stanrdump/min_n.csv', stringsAsFactors = FALSE)

total_prod <- read.csv('~/forestlight/stanrdump/production_total.csv', stringsAsFactors = FALSE)

mod_df <- mod_df %>%
  left_join(total_prod) %>%
  rename(total_production = production) %>%
  left_join(min_n)

dbh_pred <- exp(seq(log(1), log(315), length.out = 101))


registerDoParallel(cores = 8)

tmp <- foreach(i = 1:nrow(mod_df)) %dopar% {

	if (mod_df$variable[i] == 'density') {
		fit_info <- extract_density(dens_model = mod_df$dens_model[i],
									fg = mod_df$fg[i],
									year = mod_df$year[i],
									xmin = mod_df$xmin[i],
									n = mod_df$n[i],
									use_subset = FALSE)
	}
	if (mod_df$variable[i] == 'production') {
		fit_info <- extract_production(prod_model = mod_df$prod_model[i],
									fg = mod_df$fg[i],
									year = mod_df$year[i],
									xmin = mod_df$xmin[i],
									n = mod_df$n[i],
									use_subset = FALSE)
	}
	if (mod_df$variable[i] == 'total_production') {
		fit_info <- extract_totalproduction(dens_model = mod_df$dens_model[i],
									prod_model = mod_df$prod_model[i],
									fg = mod_df$fg[i],
									year = mod_df$year[i],
									xmin = mod_df$xmin[i],
									n = mod_df$n[i],
									use_subset = FALSE)
	}

	save(fit_info, file = paste0('~/forestlight/stanoutput/fitinfo/pw_info_',mod_df$variable[i],'_',i,'.r'))
	message('Fit ', i, ' saved')

}

# DENSITY - INCOMING LIGHT SCALINGS
# =================================

source('~/forestlight/stancode/model_output_extraction_functions.r')

library(purrr)
library(dplyr)
library(foreach)
library(doParallel)

prod_df <- expand.grid(variable = 'production',
					   dens_model = as.numeric(NA),
					   prod_model = 1:2,
					   fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
					   year = 1995,
					   stringsAsFactors = FALSE)
					   					   
mod_df <- expand.grid(variable = 'total_production',
					  dens_model = 1:2,
                      prod_model = 1:2,
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = 1995, 
                      stringsAsFactors = FALSE)

mod_df <- rbind(prod_df, mod_df)

# Make sure that the multiplication is done by the number of trees that have light measurements, not the total number.
min_n <- read.csv('~/forestlight/stanrdump/min_n_lighttrees.csv', stringsAsFactors = FALSE)

total_prod <- read.csv('~/forestlight/stanrdump/lightrec_total.csv', stringsAsFactors = FALSE)

mod_df <- mod_df %>%
  left_join(total_prod) %>%
  rename(total_production = light_received) %>%
  left_join(min_n)

dbh_pred <- exp(seq(log(1), log(315), length.out = 101))


registerDoParallel(cores = 8)

tmp <- foreach(i = 1:nrow(mod_df)) %dopar% {

	if (mod_df$variable[i] == 'production') {
		fit_info <- extract_production(prod_model = mod_df$prod_model[i],
									fg = mod_df$fg[i],
									year = mod_df$year[i],
									xmin = mod_df$xmin[i],
									n = mod_df$n[i],
									infix = 'rawlightscaling_',
									use_subset = FALSE)
	}
	if (mod_df$variable[i] == 'total_production') {
		fit_info <- extract_totalproduction(dens_model = mod_df$dens_model[i],
									prod_model = mod_df$prod_model[i],
									fg = mod_df$fg[i],
									year = mod_df$year[i],
									xmin = mod_df$xmin[i],
									n = mod_df$n[i],
									infix = 'rawlightscaling_',
									use_subset = FALSE)
	}

	save(fit_info, file = paste0('~/forestlight/stanoutput/fitinfo/lightpw_info_',mod_df$variable[i],'_',i,'.r'))
	message('Fit ', i, ' saved')

}

# DENSITY x CROWN VOLUME SCALING
# FITTED VALUES ONLY
# ==============================

# this one will change if and when we get individual allometries for crown volume by species. (edit at that time)

source('~/forestlight/stancode/model_output_extraction_functions.r')
source('~/forestlight/stancode/fittedcrownvolumefunction.r')

library(purrr)
library(dplyr)
library(foreach)
library(doParallel)

mod_df <- expand.grid(variable = 'total_production',
					  dens_model = 1:2,
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = 1995, 
                      stringsAsFactors = FALSE)
					  
min_n <- read.csv('~/forestlight/stanrdump/min_n.csv', stringsAsFactors = FALSE)

total_prod <- read.csv('~/forestlight/stanrdump/crownvol_total.csv', stringsAsFactors = FALSE)

mod_df <- mod_df %>%
  left_join(total_prod) %>%
  rename(total_production = crownvolume) %>%
  left_join(min_n)
  

dbh_pred <- exp(seq(log(1), log(315), length.out = 101))

density_par <- list('1' = c('alpha'),
					'2' = c('alpha_low', 'alpha_high', 'tau'))

registerDoParallel(cores = 8)

tmp <- foreach(i = 1:nrow(mod_df)) %dopar% {
	
  require(rstan)
  require(Brobdingnag)

  # Load CSVs as stanfit object
  message('Loading stan fit ', i, ' . . .')
  files <- paste0('fit_density', mod_df$dens_model[i], '_', mod_df$fg[i], '_1995_', 1:3, '.csv')
  fit <- list(density = read_stan_csv(file.path('~/forestlight/stanoutput',files)))

  fitted_totalvolume_values <- fitted_totalvolume(fit = fit[['density']], 
												  dbh_pred = dbh_pred,
												  dens_form = mod_df$dens_model[i],
												  x_min = mod_df$xmin[i],
												  n_indiv = mod_df$n[i],
												  pars_to_get = density_par[[mod_df$dens_model[i]]])
  
	save(fitted_totalvolume_values, file = paste0('~/forestlight/stanoutput/fitinfo/volumepw_fittedvalues_',i,'.r'))
	message('Fit ', i, ' saved')

}
 
 
# GROWTH AS DIAMETER PER TIME SCALING
# INDIVIDUAL PRODUCTION ONLY
# ===================================

source('~/forestlight/stancode/model_output_extraction_functions.r')

library(purrr)
library(dplyr)
library(foreach)
library(doParallel)

prod_df <- expand.grid(variable = 'production',
					   dens_model = as.numeric(NA),
					   prod_model = 1:2,
					   fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
					   year = 1995,
					   stringsAsFactors = FALSE) 
					   
min_n <- read.csv('~/forestlight/stanrdump/min_n.csv', stringsAsFactors = FALSE)	


mod_df <- prod_df %>%
  left_join(min_n)

dbh_pred <- exp(seq(log(1), log(315), length.out = 101))

registerDoParallel(cores = 8)

tmp <- foreach(i = 1:nrow(mod_df)) %dopar% {
		fit_info <- extract_production(prod_model = mod_df$prod_model[i],
									fg = mod_df$fg[i],
									year = mod_df$year[i],
									xmin = mod_df$xmin[i],
									n = mod_df$n[i],
									infix = 'diamgrowthscaling_',
									use_subset = FALSE)
	
	save(fit_info, file = paste0('~/forestlight/stanoutput/fitinfo/diamgrowthpw_info_',i,'.r'))
	message('Fit ', i, ' saved')

}				   

### note: mortality and light extraction of output are in separate scripts.

