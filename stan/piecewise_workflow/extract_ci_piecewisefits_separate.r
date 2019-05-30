# For each CMDstan model fit, get the following
# Credible intervals of parameters
# Credible and prediction intervals to make graphs
# Information criteria
# Fitted slopes in log space
# Bayesian R-squared


# Alternate version created 30 May 2019 for fits where density and production are done separately.

source('~/forestlight/stancode/extraction_functions_piecewise_separate.r')

library(purrr)
library(dplyr)

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
					   					   
mod_df <- expand.grid(variable = 'total_production',
					  dens_model = 1:3,
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

dbh_pred <- exp(seq(log(1.2), log(315), length.out = 101))


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
								total_production = mod_df$total_production[i],
								use_subset = FALSE)
}
if (mod_df$variable[i] == 'total_production') {
	fit_info <- extract_totalproduction(dens_model = mod_df$dens_model[i],
								prod_model = mod_df$prod_model[i],
								fg = mod_df$fg[i],
								year = mod_df$year[i],
								xmin = mod_df$xmin[i],
								n = mod_df$n[i],
								total_production = mod_df$total_production[i],
								use_subset = FALSE)
}

save(fit_info, file = paste0('~/forestlight/stanoutput/fitinfo/pw_info_'mod_df$variable[i],'_',i,'.r'))
