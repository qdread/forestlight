
# Load stan fits ----------------------------------------------------------

library(rstan)
library(bayesplot)
fp <- 'C:/Users/Q/Dropbox/projects/forestlight/stanoutput'

ppow_fits <- list()
wpow_fits <- list()
pexp_fits <- list()
wexp_fits <- list()
fg_names <- c('fg1','fg2','fg3','fg4','fg5','alltree','unclassified')

for (i in 1:length(fg_names)) {
  ppow_fits[[i]] <- read_stan_csv(file.path(fp, paste0('fit_paretoxpower_', fg_names[i], '_1995_', 1:3, '.csv')))
  wpow_fits[[i]] <- read_stan_csv(file.path(fp, paste0('fit_weibullxpower_', fg_names[i], '_1995_', 1:3, '.csv')))
}

for (i in 1:length(fg_names)) {
  prefix <- ifelse(i %in% 5:6, 'ssfit_', 'fit_')
  pexp_fits[[i]] <- read_stan_csv(file.path(fp, paste0(prefix, 'paretoxexp_', fg_names[i], '_1995_', 1:3, '.csv')))
  wexp_fits[[i]] <- read_stan_csv(file.path(fp, paste0(prefix, 'weibullxexp_', fg_names[i], '_1995_', 1:3, '.csv')))
}

# For fits run locally:
ppow_fits <- c(fit_ppow_all, fit_ppow_fg)
wpow_fits <- c(fit_wpow_all, fit_wpow_fg)
pexp_fits <- c(fit_pexp_all, fit_pexp_fg)
wexp_fits <- c(fit_wexp_all, fit_wexp_fg)
fg_names <- c('alltree', 'fg1','fg2','fg3','fg4','fg5','unclassified')

# Run diagnostics ---------------------------------------------------------

source('stan/extract_ci_stan.r')

#diag_plots(ppow_fits[[1]])
par_pars <- c('alpha')
wei_pars <- c('m', 'n')
pow_pars <- c('beta0', 'beta1', 'sigma')
exp_pars <- c('a', 'b', 'c')

ppow_trace <- lapply(ppow_fits, function(x) mcmc_trace(as.array(x), pars = c(par_pars, pow_pars)))
wpow_trace <- lapply(wpow_fits, function(x) mcmc_trace(as.array(x), pars = c(wei_pars, pow_pars)))
pexp_trace <- lapply(pexp_fits, function(x) mcmc_trace(as.array(x), pars = c(par_pars, pow_pars, exp_pars)))
wexp_trace <- lapply(wexp_fits, function(x) mcmc_trace(as.array(x), pars = c(wei_pars, pow_pars, exp_pars)))

ppow_summ <- lapply(ppow_fits, summary)
wpow_summ <- lapply(wpow_fits, summary)
pexp_summ <- lapply(pexp_fits, summary)
wexp_summ <- lapply(wexp_fits, summary)

# Extract intervals from fits --------------------------------------------

# Confidence interval plots.
# Load density bin midpoints to get x values at which to predict.
densitybin_5census <- read.csv('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/densitybin_5census.csv', stringsAsFactors = FALSE)
dbh_pred_bins <- densitybin_5census[densitybin_5census$fg == 'all', 'bin_midpoint']

# Minimum x values for Pareto.
min_n <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/stanoutput/min_n.csv', stringsAsFactors = FALSE)

# Get confidence intervals
ppow_cis <- list()
wpow_cis <- list()
pexp_cis <- list()
wexp_cis <- list()

for (i in 1:length(ppow_fits)) {
  x_min_i <- min_n$xmin[min_n$year == 1995 & min_n$fg == fg_names[i]]
  n_i <- min_n$n[min_n$year == 1995 & min_n$fg == fg_names[i]]
  ppow_cis[[i]] <- dens_prod_ci(ppow_fits[[i]], dbh_pred = dbh_pred_bins, dens_form = 'pareto', prod_form = 'powerlaw', x_min = x_min_i, n_indiv = n_i)
  wpow_cis[[i]] <- dens_prod_ci(wpow_fits[[i]], dbh_pred = dbh_pred_bins, dens_form = 'weibull', prod_form = 'powerlaw', x_min = NULL, n_indiv = n_i)
  pexp_cis[[i]] <- dens_prod_ci(pexp_fits[[i]], dbh_pred = dbh_pred_bins, dens_form = 'pareto', prod_form = 'powerlawexp', x_min = x_min_i, n_indiv = n_i)
  wexp_cis[[i]] <- dens_prod_ci(wexp_fits[[i]], dbh_pred = dbh_pred_bins, dens_form = 'weibull', prod_form = 'powerlawexp', x_min = NULL, n_indiv = n_i)
}

# Make into one data frame.
ppow_ci_df <- data.frame(dens_model = 'pareto', prod_model = 'powerlaw', 
                           fg = rep(fg_names, times = sapply(ppow_cis, nrow)),
                           do.call('rbind', ppow_cis))
wpow_ci_df <- data.frame(dens_model = 'weibull', prod_model = 'powerlaw', 
                         fg = rep(fg_names, times = sapply(wpow_cis, nrow)),
                         do.call('rbind', wpow_cis))
pexp_ci_df <- data.frame(dens_model = 'pareto', prod_model = 'powerlawexp', 
                         fg = rep(fg_names, times = sapply(pexp_cis, nrow)),
                         do.call('rbind', pexp_cis))
wexp_ci_df <- data.frame(dens_model = 'weibull', prod_model = 'powerlawexp', 
                           fg = rep(fg_names, times = sapply(wexp_cis, nrow)),
                           do.call('rbind', wexp_cis))

ci_df <- rbind(ppow_ci_df, wpow_ci_df, pexp_ci_df, wexp_ci_df)

# Plot fits only ----------------------------------------------------------

library(cowplot)

area_core <- 42.84

ggplot(dplyr::filter(ci_df, fg == 'alltree', variable == 'density', prod_model == 'powerlaw'),
       aes(x = dbh, color = dens_model, group = dens_model)) +
  scale_x_log10() + scale_y_log10() +
  geom_line(aes(y = q50/area_core)) +
  geom_line(aes(y = q025/area_core), linetype = 'dotted') +
  geom_line(aes(y = q975/area_core), linetype = 'dotted')

ggplot(dplyr::filter(ci_df, fg == 'alltree', variable == 'production', dens_model == 'pareto'),
       aes(x = dbh, color = prod_model, group = prod_model)) +
  scale_x_log10() + scale_y_log10() +
  geom_line(aes(y = q50)) +
  geom_line(aes(y = q025), linetype = 'dotted') +
  geom_line(aes(y = q975), linetype = 'dotted')

ggplot(dplyr::filter(ci_df, fg == 'alltree', variable == 'total_production'),
       aes(x = dbh, color = interaction(dens_model, prod_model), group = interaction(dens_model, prod_model))) +
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

fakebin_across_years <- function(dat_values, dat_classes, edges, mean_type = 'geometric', n_census = 5) {
  qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  # add some padding just in case
  mins <- edges$bin_min
  mins[1] <- 0
  maxes <- edges$bin_max
  maxes[length(maxes)] <- Inf
  
  binstats <- t(sapply(1:length(mins), function(x) {
    indivs <- dat_values[dat_classes >= mins[x] & dat_classes < maxes[x]]
    if (mean_type == 'geometric') {
      mean_n <- exp(mean(log(indivs)))
      sd_n <- exp(sd(log(indivs)))
      ci_width <- qnorm(0.975) * sd(log(indivs)) / sqrt(length(indivs))
      ci_min <- exp(mean(log(indivs)) - ci_width)
      ci_max <- exp(mean(log(indivs)) + ci_width)
      
    } else {
      mean_n <- mean(indivs)
      sd_n <- sd(indivs)
      ci_width <- qnorm(0.975) * sd(indivs) / sqrt(length(indivs))
      ci_min <- mean_n - ci_width
      ci_max <- mean_n + ci_width
    }
    c(mean = mean_n, 
      sd = sd_n,
      quantile(indivs, probs = qprobs),
      ci_min = ci_min,
      ci_max = ci_max)
  }))
  dimnames(binstats)[[2]] <- c('mean', 'sd', 'q025', 'q25', 'median', 'q75', 'q975', 'ci_min', 'ci_max')
  data.frame(bin_midpoint = edges$bin_midpoint,
             bin_min = edges$bin_min,
             bin_max = edges$bin_max,
             mean_n_individuals = edges$bin_count / n_census,
             binstats)
}



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

fpfig <- 'C:/Users/Q/google_drive/ForestLight/figs/credible_interval_plots'

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.
global_diam_xlimits <- c(1, 316) # Maximum and minimum for x axis if it's diameter.

p_all_indivprod <- all_prod_1995 %>%
  filter(!is.na(mean), mean > 0) %>%
  ggplot() +
  geom_point(aes(x = bin_midpoint, y = mean), color = 'black') +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
  panel_border(colour = 'black')

dat_prod <- dplyr::filter(ci_df, fg == 'alltree', variable == 'production', dens_model == 'pareto')

p_all_indivprod +
  geom_ribbon(data = dat_prod, aes(ymin = q025, ymax = q975, x = dbh, fill = prod_model, group = prod_model), alpha = 0.5) +
  geom_line(data = dat_prod, aes(y = q50, x = dbh, color = prod_model, group = prod_model)) +
  theme(legend.position = 'bottom') + ggtitle('Individual production')
ggsave(file.path(fpfig, 'alltree_1995_production.png'), height=5, width=5, dpi=400)


p_all_dens <- all_dens_1995 %>%
  filter(!is.na(bin_value), bin_value > 0) %>%
  mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
  mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
  ggplot() +
  geom_point(aes(x = bin_midpoint, y = bin_yvalue), color = 'black') +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  panel_border(colour = 'black')

dat_dens <- dplyr::filter(ci_df, fg == 'alltree', variable == 'density', prod_model == 'powerlaw')

p_all_dens +
  geom_ribbon(data = dat_dens, aes(ymin = q025/area_core, ymax = q975/area_core, x = dbh, fill = dens_model, group = dens_model), alpha = 0.5) +
  geom_line(data = dat_dens, aes(y = q50/area_core, x = dbh, color = dens_model, group = dens_model)) +
  theme(legend.position = 'bottom') + ggtitle('Density')
ggsave(file.path(fpfig, 'alltree_1995_density.png'), height=5, width=5, dpi=400)


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
  geom_ribbon(data = dat_totalprod, aes(ymin = q025/area_core, ymax = q975/area_core, x = dbh, fill = interaction(dens_model, prod_model), group = interaction(dens_model, prod_model)), alpha = 0.5) +
  geom_line(data = dat_totalprod, aes(y = q50/area_core, x = dbh, color = interaction(dens_model, prod_model), group = interaction(dens_model, prod_model))) +
  scale_fill_discrete(name = 'Functional forms', labels = c('D Pareto\nP Powerlaw', 'D Weibull\nP Powerlaw', 'D Pareto\nP Powerlaw*Exponential', 'D Weibull\nP Powerlaw*Exponential')) +
  theme(legend.position = 'bottom', legend.text = element_text(size=8)) + ggtitle('Total production') + guides(colour = FALSE)
ggsave(file.path(fpfig, 'alltree_1995_totalproduction.png'), height=5, width=6.5, dpi=400)
