
# Load stan fits ----------------------------------------------------------

library(rstan)
library(bayesplot)
fp <- 'C:/Users/Q/Dropbox/projects/forestlight/stanoutput'

pareto_fits <- list()
weibull_fits <- list()
fg_names <- c('fg1','fg2','fg3','fg4','fg5','alltree','unclassified')

for (i in 1:length(fg_names)) {
  pareto_fits[[i]] <- read_stan_csv(file.path(fp, paste0('fit_paretoxpower_', fg_names[i], '_1995_', 1:3, '.csv')))
}

for (i in 1:length(fg_names)) {
  prefix <- ifelse(i %in% 5:6, 'ssfit_weibullxexp_', 'fit_weibullxexp_')
  weibull_fits[[i]] <- read_stan_csv(file.path(fp, paste0(prefix, fg_names[i], '_1995_', 1:3, '.csv')))
}


# Run diagnostics ---------------------------------------------------------

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
mcmc_trace(as.array(weibull_fits[[4]]), pars = weib_pars)
mcmc_trace(as.array(weibull_fits[[5]]), pars = weib_pars)
mcmc_trace(as.array(weibull_fits[[6]]), pars = weib_pars)
mcmc_trace(as.array(weibull_fits[[7]]), pars = weib_pars)


# Extract intervals from fits --------------------------------------------

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

# Plot fits only ----------------------------------------------------------

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

# Plot data and fits ------------------------------------------------------
# Get density, production, and total production for all trees in 1995 for plotting.

library(dplyr)

load('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/rawdataobj_22jan.r')
group_names <- c('all','all_classified','fg1','fg2','fg3','fg4','fg5','unclassified')
source('code/allfunctions27july.r')

# Set number of bins       
numbins <- 20

# Bin all trees including unclassified
allyeardbh <- unlist(lapply(alltreedat[2:6], '[', , 'dbh_corr'))
dbhbin_all <- logbin(x = allyeardbh, y = NULL, n = numbins)
dbhbin_all_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_all))

all_dens_1995 <- dbhbin_all_byyear[[2]]

all_prod_1995 <- fakebin_across_years(dat_values = alltreedat[[3]]$production, dat_classes = alltreedat[[3]]$dbh_corr, edges = all_dens_1995, n_census = 1)

totalprodbin_alltree_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_all))

all_totalprod_1995 <- totalprodbin_alltree_byyear[[2]]
# Plot all trees' individual production, density, and total production.

# Set options for error bar widths and dodge amounts for all the plots
error_bar_width <- 0.03
dodge_width <- 0.03

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.
global_diam_xlimits <- c(1, 316) # Maximum and minimum for x axis if it's diameter.

p_all_indivprod <- all_prod_1995 %>%
  filter(!is.na(mean), mean > 0) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_point(color = 'black') +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
  panel_border(colour = 'black')

dat_prod <- dplyr::filter(ci_df, fg == 'alltree', variable == 'production')

p_all_indivprod +
  geom_line(data = dat_prod, aes(y = q50/area_core, x = dbh, color = prod_model, group = prod_model)) +
  geom_line(data = dat_prod, aes(y = q025/area_core, x = dbh, color = prod_model, group = prod_model), linetype = 'dotted') +
  geom_line(data = dat_prod, aes(y = q975/area_core, x = dbh, color = prod_model, group = prod_model), linetype = 'dotted')


p_all_dens <- all_dens_1995 %>%
  filter(!is.na(bin_value), bin_value > 0) %>%
  mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
  mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
  ggplot() +
  geom_point(aes(x = bin_midpoint, y = bin_yvalue), color = 'black') +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  panel_border(colour = 'black')

dat_dens <- dplyr::filter(ci_df, fg == 'alltree', variable == 'density')

# Correction factor for subsample
dens_corr_factor <- sum(all_dens_1995$bin_count)/25000

p_all_dens +
  geom_line(data = dat_dens, aes(y = dens_corr_factor*q50/area_core, x = dbh, color = dens_model, group = dens_model)) +
  geom_line(data = dat_dens, aes(y = dens_corr_factor*q025/area_core, x = dbh, color = dens_model, group = dens_model), linetype = 'dotted') +
  geom_line(data = dat_dens, aes(y = dens_corr_factor*q975/area_core, x = dbh, color = dens_model, group = dens_model), linetype = 'dotted')


p_all_totalprod <- all_totalprod_1995 %>%
  filter(!is.na(bin_value), bin_value > 0) %>%
  mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
  mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
  ggplot() +
  geom_point(aes(x = bin_midpoint, y = bin_yvalue)) +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total production (kg ha'^-1,' y'^-1,')')), limits = c(.1,10000)) +
  panel_border(colour = 'black')

dat_totalprod <- dplyr::filter(ci_df, fg == 'alltree', variable == 'total_production')

p_all_totalprod +
  geom_line(data = dat_totalprod, aes(y = dens_corr_factor*q50/area_core, x = dbh, color = dens_model, group = dens_model)) +
  geom_line(data = dat_totalprod, aes(y = dens_corr_factor*q025/area_core, x = dbh, color = dens_model, group = dens_model), linetype = 'dotted') +
  geom_line(data = dat_totalprod, aes(y = dens_corr_factor*q975/area_core, x = dbh, color = dens_model, group = dens_model), linetype = 'dotted') +
  scale_color_discrete(name = 'Functional forms', labels = c('D Pareto\nP Powerlaw', 'D Weibull\nP Powerlaw + exponential')) +
  theme(legend.position = 'bottom')
