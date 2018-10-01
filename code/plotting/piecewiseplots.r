# Plot of slopes in different segments by different functional groups.

fp <- '~/google_drive/ForestLight/data/data_piecewisefits'
fpfig <- '~/google_drive/ForestLight/figs/piecewiseplots_27sep2018'
ics <- read.csv(file.path(fp, 'piecewise_ics_by_fg.csv'), stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)

# Density model

ggplot(ics %>% filter(prod_model == 1, criterion == 'LOOIC', variable == 'density', !fg %in% 'unclassified'), 
       aes(x = factor(dens_model), y = ic, ymin = ic - se_ic, ymax = ic + se_ic)) +
  facet_wrap(~ fg, scales = 'free_y') +
  geom_pointrange() +
  theme_bw() +
  theme(strip.background = element_blank()) +
  labs(x = 'Number of segments in density function')
ggsave(file.path(fpfig, 'density_model_information_criteria.pdf'), height = 6, width = 9)

# Production model

ggplot(ics %>% filter(dens_model == 1, criterion == 'LOOIC', variable == 'production', !fg %in% 'unclassified'), aes(x = factor(prod_model), y = ic, ymin = ic - se_ic, ymax = ic + se_ic)) +
  facet_wrap(~ fg, scales = 'free_y') +
  geom_pointrange() +
  theme_bw() +
  theme(strip.background = element_blank()) +
  labs(x = 'Number of segments in production function')
ggsave(file.path(fpfig, 'production_model_information_criteria.pdf'), height = 6, width = 9)

# Plot the slopes of each one.

slopes <- read.csv(file.path(fp, 'piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)

# Using 3 segment density and 2 segment production
ggplot(slopes %>% filter(dens_model == 3, prod_model == 2, !fg %in% 'unclassified'), 
       aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  scale_x_log10(name = 'Diameter (cm)', expand = c(0,0)) +
  theme_bw() +
  theme(strip.background = element_blank(), panel.grid = element_blank(), legend.position = 'bottom') +
  labs(y = 'Fitted slope in log-log space') +
  ggtitle('Fitted slopes', '3 segment density model and 2 segment production model')
ggsave(file.path(fpfig, 'fitted_slopes_3partdensity_2partproduction.pdf'), height = 6, width = 9)

ggplot(slopes %>% filter(dens_model == 3, prod_model == 1, !fg %in% 'unclassified'), 
       aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  scale_x_log10(name = 'Diameter (cm)', expand = c(0,0)) +
  theme_bw() +
  theme(strip.background = element_blank(), panel.grid = element_blank(), legend.position = 'bottom') +
  labs(y = 'Fitted slope in log-log space') +
  ggtitle('Fitted slopes', '3 segment density model and 1 segment production model')

# Plot the fitted values on top of the observed histograms.

# Read observed data
fp_obs <- '~/google_drive/ForestLight/data/data_forplotting_aug2018'

for (i in dir(fp_obs, pattern = 'obs_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp_obs, i), stringsAsFactors = FALSE))
}

# Read modeled data (CIs)

for (i in dir(fp, pattern = 'pred_|fitted_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}

source('stan/piecewise_workflow/plottingfunctionspiecewise.r')

# Create plots.

plot_dens(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = 3,
          y_limits = c(0.001, 1000),
          y_breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000))
ggsave(file.path(fpfig, 'fits_3partdensity.pdf'), height = 5, width = 6)

# Specify dodging with a certain width of error bar
plot_prod(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = 1,
          y_limits = c(0.01, 6000),
          y_breaks = c(0.1, 10, 1000),
          error_bar_width = 0.01,
          dodge_width = 0.05)

plot_prod(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = 2,
          y_limits = c(0.01, 6000),
          y_breaks = c(0.1, 10, 1000),
          error_bar_width = 0.01,
          dodge_width = 0.05)
ggsave(file.path(fpfig, 'fits_2partproduction.pdf'), height = 5, width = 6)

plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
               model_fit_density = 3, 
               model_fit_production = 2,
               y_limits = c(0.01, 200),
               y_breaks = c(0.1, 1, 10, 100),
               preddat = fitted_totalprod)
ggsave(file.path(fpfig, 'fits_3by2_totalproduction.pdf'), height = 5, width = 6)

plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
               model_fit_density = 3, 
               model_fit_production = 1,
               y_limits = c(0.01, 200),
               y_breaks = c(0.1, 1, 10, 100))
