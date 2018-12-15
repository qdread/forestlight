# Plot of slopes in different segments by different functional groups.
library(tidyverse)
john_wd <- "/Users/jgradym/Google Drive/ForestLight"
setwd(john_wd)

#fp <- '~/google_drive/ForestLight/data/data_piecewisefits'
fp <- 'data/data_piecewisefits'

#fpfig <- '~/google_drive/ForestLight/figs/piecewiseplots_27sep2018'
fpfig <- 'figs/piecewiseplots_27sep2018'

ics <- read.csv(file.path(fp, 'piecewise_ics_by_fg.csv'), stringsAsFactors = FALSE)
ics$fg <- factor(ics$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))

# Density model

ggplot(ics %>% filter(prod_model == 1, criterion == 'LOOIC', variable == 'density', !fg %in% 'unclassified'), 
       aes(x = factor(dens_model), y = ic, ymin = ic - se_ic, ymax = ic + se_ic)) +
  facet_wrap(~ fg, labeller = label_value, scales = 'free_y') +
  geom_pointrange() +
  theme_bw() +
  theme(strip.background = element_blank(),panel.grid = element_blank()) +
  labs(x = 'Number of Segments in Density Function')
ggsave(file.path(fpfig, 'density_model_information_criteria.pdf'), height = 6, width = 9)

# Production model

ggplot(ics %>% filter(dens_model == 1, criterion == 'LOOIC', variable == 'production', !fg %in% 'unclassified'), aes(x = factor(prod_model), y = ic, ymin = ic - se_ic, ymax = ic + se_ic)) +
  facet_wrap(~ fg, labeller = label_value,scales = 'free_y') +
  geom_pointrange() +
  theme_bw() +
  theme(strip.background = element_blank(), panel.grid = element_blank()) +
  labs(x = 'Number of Segments in Production Function')
ggsave(file.path(fpfig, 'production_model_information_criteria.pdf'), height = 6, width = 9)

# Plot the slopes of each one.

slopes <- read.csv(file.path(fp, 'piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
slopes$fg <- factor(slopes$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))
slopes$variable <- factor(slopes$variable, labels = c("Density", "Individual Growth", "Total Growth"))
colors <- c("sienna4", "yellowgreen", "springgreen4")

# Using 3 segment density and 2 segment production
ggplot(slopes %>% filter(dens_model == 3, prod_model == 1, !fg %in% 'Unclassified'), 
       aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg,labeller = label_value) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = "springgreen4", size = 0.3) +
  geom_hline(yintercept = -2, linetype = 'dashed', col = "sienna4", size = 0.3) +
  geom_hline(yintercept = 2, linetype = 'dashed', col = "yellowgreen", size = 0.3) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_ribbon(alpha = 0.4, size = 0.3) +
  geom_line() +
  scale_x_log10(name = 'Diameter (cm)', expand = c(0,0)) +
  theme_bw() + theme(axis.text = element_text(color = "black"))+
  theme(strip.background = element_blank(), panel.grid = element_blank(), legend.position = 'bottom') +
  labs(y = 'Fitted Slope') +  #coord_fixed(ratio = .1)+
  ggtitle('Fitted slopes \n 3 segment density model and 2 segment production model')+
  theme(plot.title = element_text(hjust=0.5)) #+ theme(legend.title=element_blank())

ggsave(file.path(fpfig, 'fitted_slopes_3partdensity_2partproduction.pdf'), height = 6, width = 9)

ggplot(slopes %>% filter(dens_model == 3, prod_model == 2, !fg %in% 'Unclassified'), 
       aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg,scale = "free_y", labeller = label_value) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = "springgreen4", size = 0.3) +
  geom_hline(yintercept = -2, linetype = 'dashed', col = "sienna4",size = 0.3) +
  geom_hline(yintercept = 2, linetype = 'dashed', col = "yellowgreen",size = 0.3) +
  geom_ribbon(alpha = 0.4, size = 0.3) +
  geom_line() +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_x_log10(name = 'Diameter (cm)', expand = c(0,0)) +
  theme_bw() + theme(axis.text = element_text(color = "black"))+
  theme(strip.background = element_blank(), panel.grid = element_blank(), legend.position = 'bottom') +
  labs(y = 'Slope') +
  ggtitle('Fitted Slopes \n 3 Segment Density Model & 2 Segment Production Model')+
  theme(plot.title = element_text(hjust=0.5)) #+ theme(legend.title=element_blank())

# Plot the fitted values on top of the observed histograms.

# Read observed data
#fp_obs <- '~/google_drive/ForestLight/data/data_forplotting_aug2018'
fp_obs <- 'data/data_forplotting_aug2018'

for (i in dir(fp_obs, pattern = 'obs_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp_obs, i), stringsAsFactors = FALSE))
}

# Read modeled data (CIs)

for (i in dir(fp, pattern = 'pred_|fitted_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}

#source('stan/piecewise_workflow/plottingfunctionspiecewise.r')
source('/Users/jgradym/Documents/GitHub/forestlight/stan/piecewise_workflow/plottingfunctionspiecewise.r')
# Create plots.
#Model fit 1 = pareto, 1 segment
#Model Fit 2  = 2 segments, etc

plot_dens(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = 3,
          x_limits = c(1, 260),
          y_limits = c(0.001, 3000),
          y_labels = c(0.001, 0.1, 10,1000),
          y_breaks = c(0.001, 0.1,  10, 1000))
ggsave(file.path(fpfig, 'fits_3partdensity.pdf'), height = 5, width = 6)

# Specify dodging with a certain width of error bar
# Model fit 1 = power law
# Model fit 2 = power law exp
plot_prod(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = 2,
          x_limits = c(1, 280),
          y_limits = c(0.001, 2000),
          y_breaks = c(0.001,0.1, 10, 1000),
          y_labels = c(0.001,0.1,10,1000),
          error_bar_width = 0.01,
          dodge_width = 0.05)
ggsave(file.path(fpfig, 'fits_2partproduction.pdf'), height = 5, width = 6)

plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
               model_fit_density = 3, 
               model_fit_production = 2,
               x_limits = c(0.9,250),
               y_limits = c(0.03, 200),
               y_breaks = c(0.1, 1, 10, 100),
               y_labels = c(0.1, 1, 10, 100),
               preddat = fitted_totalprod)
ggsave(file.path(fpfig, 'fits_3by2_totalproduction.pdf'), height = 5, width = 6)
