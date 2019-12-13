# For each CMDstan model fit, get the following
# Credible intervals of parameters
# Credible and prediction intervals to make graphs
# Information criteria
# Fitted slopes in log space
# Bayesian R-squared

# New version created 10 Sep: only 1995 for piecewise.
# Edit 10 Aug: specify subset for alltree and for fg 3
# Edit 12 April: Do this in parallel because it is too slow
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

source('~/forestlight/stancode/extraction_functions_piecewise.r')

library(purrr)
library(dplyr)

mod_df <- expand.grid(dens_model = 1:3,
                      prod_model = 1:2,
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = 1995, 
                      stringsAsFactors = FALSE)

min_n <- read.csv('~/forestlight/stanrdump/min_n.csv', stringsAsFactors = FALSE)

total_prod <- read.csv('~/forestlight/stanrdump/production_total.csv', stringsAsFactors = FALSE)

mod_df <- mod_df %>%
  left_join(total_prod) %>%
  rename(total_production = production) %>%
  left_join(min_n)

dbh_pred <- exp(seq(log(1.2), log(315), length.out = 101))

use_subset <- (mod_df$fg[task] %in% c('fg3', 'alltree'))

fit_info <- extract_all_fit(dens_model = mod_df$dens_model[task],
                            prod_model = mod_df$prod_model[task],
                            fg = mod_df$fg[task],
                            year = mod_df$year[task],
                            xmin = mod_df$xmin[task],
                            n = mod_df$n[task],
                            total_production = mod_df$total_production[task],
							use_subset = use_subset)

save(fit_info, file = paste0('~/forestlight/stanoutput/fitinfo/pw_info_',task,'.r'))
