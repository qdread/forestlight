library(rstan)
fp <- 'C:/Users/Q/Dropbox/projects/forestlight/stanoutput'

pareto_fits <- list()
weibull_fits <- list()
fg_names <- c('fg1','fg2','fg3','fg4','fg5','alltree','unclassified')

for (i in 1:length(fg_names)) {
  pareto_fits[[i]] <- read_stan_csv(file.path(fp, paste0('fit_paretoxpower_', fg_names[i], '_1995_', 1:3, '.csv')))
}

for (i in 1:3) {
  weibull_fits[[i]] <- read_stan_csv(file.path(fp, paste0('fit_weibullxexp_', fg_names[i], '_1995_', 1:3, '.csv')))
}

source('stan/extract_ci_stan.r')

diag_plots(pareto_fits[[1]])
pareto_pars <- names(pareto_fits[[2]])[-length(names(pareto_fits[[2]]))]
weib_pars <- names(weibull_fits[[2]])[-length(names(weibull_fits[[2]]))]

mcmc_trace(as.array(pareto_fits[[1]]), pars = pareto_pars)
mcmc_trace(as.array(pareto_fits[[2]]), pars = pareto_pars)
mcmc_trace(as.array(pareto_fits[[3]]), pars = pareto_pars)
mcmc_trace(as.array(pareto_fits[[4]]), pars = pareto_pars)
mcmc_trace(as.array(pareto_fits[[5]]), pars = pareto_pars)
mcmc_trace(as.array(pareto_fits[[6]]), pars = pareto_pars)
mcmc_trace(as.array(pareto_fits[[7]]), pars = pareto_pars)

mcmc_trace(as.array(weibull_fits[[1]]), pars = weib_pars) 
mcmc_trace(as.array(weibull_fits[[2]]), pars = weib_pars) 
mcmc_trace(as.array(weibull_fits[[3]]), pars = weib_pars) 

# Confidence interval plots.
# Load density bin midpoints to get x values at which to predict.
densitybin_5census <- read.csv('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/densitybin_5census.csv', stringsAsFactors = FALSE)
dbh_pred_bins <- densitybin_5census[densitybin_5census$fg == 'all', 'bin_midpoint']

# Minimum x values for Pareto.
min_n <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/stanoutput/min_n.csv', stringsAsFactors = FALSE)

# Get confidence intervals
xmins <- min_n$xmin[min_n$year == 1995]
ns <- min_n$n[min_n$year == 1995]
pareto_cis <- list()
weibull_cis <- list()
for (i in 1:length(pareto_fits)) {
  pareto_cis[[i]] <- dens_prod_ci(pareto_fits[[i]], dbh_pred = dbh_pred_bins, dens_form = 'pareto', prod_form = 'powerlaw', x_min = xmins[i], n_indiv = ns[i])
  weibull_cis[[i]] <- dens_prod_ci(weibull_fits[[i]], dbh_pred = dbh_pred_bins, dens_form = 'weibull', prod_form = 'powerlawexp', x_min = NULL, n_indiv = ns[i])
}

# Make into one data frame.
pareto_ci_df <- data.frame(dens_model = 'pareto', prod_model = 'powerlaw', 
                           fg = rep(fg_names, times = sapply(pareto_cis, nrow)),
                           do.call('rbind', pareto_cis))
weibull_ci_df <- data.frame(dens_model = 'weibull', prod_model = 'powerlawexp', 
                           fg = rep(fg_names, times = sapply(weibull_cis, nrow)),
                           do.call('rbind', weibull_cis))

ci_df <- rbind(pareto_ci_df, weibull_ci_df)

library(cowplot)

area_core <- 42.84

ggplot(dplyr::filter(ci_df, fg == 'alltree', variable == 'density'),
       aes(x = dbh, color = dens_model, group = dens_model)) +
  scale_x_log10() + scale_y_log10() +
  geom_line(aes(y = q50/area_core)) +
  geom_line(aes(y = q025/area_core), linetype = 'dotted') +
  geom_line(aes(y = q975/area_core), linetype = 'dotted')

ggplot(dplyr::filter(ci_df, fg == 'alltree', variable == 'production'),
       aes(x = dbh, color = prod_model, group = prod_model)) +
  scale_x_log10() + scale_y_log10() +
  geom_line(aes(y = q50)) +
  geom_line(aes(y = q025), linetype = 'dotted') +
  geom_line(aes(y = q975), linetype = 'dotted')

ggplot(dplyr::filter(ci_df, fg == 'alltree', variable == 'total_production'),
       aes(x = dbh, color = dens_model, group = dens_model)) +
  scale_x_log10() + scale_y_log10(limits = c(10, 1000), breaks = c(10,100,1000)) +
  geom_line(aes(y = q50/area_core)) +
  geom_line(aes(y = q025/area_core), linetype = 'dotted') +
  geom_line(aes(y = q975/area_core), linetype = 'dotted')
