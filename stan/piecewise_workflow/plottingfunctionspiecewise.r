
# Plotting functions for piecewise.
# Define plotting functions -----------------------------------------------
# Plotting functions for piecewise.
library(tidyverse)

# Plot single model fit with multiple functional groups for density
plot_dens <- function(year_to_plot = 1995,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                      model_fit = 'pareto',
                      x_limits,
                      x_breaks = c(1, 3, 10, 30, 100,300),
                      y_limits,
                      y_breaks,
                      y_labels,
                      color_names = c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "ivory"),
                      x_name = 'Diameter (cm)',
                      y_name = expression(paste('Density (trees ha'^-1,'cm'^-1,')')),
                      obsdat = obs_dens,
                      preddat = pred_dens
) {
  
  require(dplyr)
  require(ggplot2)
  
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot, bin_count > 1) %>%
    filter(bin_value > 0)
  
  # Get minimum and maximum observed bin value for each group to be plotted
  # Delete points on the predicted line that are outside of this range
  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(dens_model %in% model_fit, prod_model == 2, fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs) 
  
  ggplot() +
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    geom_line(data = preddat[preddat$fg == "fg5",], aes(x = dbh, y = q50), color = "gray")+ # white circles get gray line
    geom_ribbon(data = preddat[preddat$fg == "fg5",], aes(x = dbh, ymin = q025, ymax = q975), fill = "gray", alpha = 0.4) +
    geom_point(data = obsdat, aes(x = bin_midpoint, y = bin_value, group = fg, fill=fg,size=2), shape=21,color="black") +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks,labels = y_labels) +
    scale_color_manual(values = color_names) +theme_plant+
    scale_fill_manual(values = color_names) 
  
  
}

# Plot single model fit with multiple functional groups for individual production
plot_prod <- function(year_to_plot = 1995,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                      model_fit = 'powerlaw',
                      x_limits,
                      x_breaks = c(1, 3, 10, 30, 100,300),
                      y_limits,
                      y_labels,
                      y_breaks,
                      color_names = c("#BFE046", "#267038", "#27408b", "#87Cefa", "ivory"),
                      x_name = 'Diameter (cm)',
                      y_name = expression(paste('Production (kg y'^-1,')')),
                      average = 'mean',
                      error_quantiles = c('ci_min', 'ci_max'),
                      error_bar_width = 0.03,
                      dodge_width = 0.03,
                      dodge_errorbar = TRUE,
                      obsdat = obs_indivprod,
                      preddat = fitted_indivprod
) {
  
  require(dplyr)
  require(ggplot2)
  
  pos <- if (dodge_errorbar) position_dodge(width = dodge_width) else 'identity'
  
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot, !is.na(mean), mean_n_individuals > 1) %>%
    group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup
  
  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(prod_model %in% model_fit, dens_model == 1, fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs)
  
  ggplot() +
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    #geom_errorbar(data = obsdat, aes_string(x = 'bin_midpoint', ymin = error_quantiles[1], ymax = error_quantiles[2], group = 'fg', color = 'fg', width = 'width'), position = pos) +
    geom_line(data = preddat[preddat$fg == "fg5",], aes(x = dbh, y = q50), color = "gray")+ # white circles get gray line
    geom_ribbon(data = preddat[preddat$fg == "fg5",], aes(x = dbh, ymin = q025, ymax = q975), fill = "gray", alpha = 0.4) +
    geom_point(data = obsdat, aes_string(x = 'bin_midpoint', y = average, group = 'fg', fill = 'fg'),size=4,color="black",shape=21,position = pos) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels=y_labels) +
    scale_color_manual(values = color_names) +theme_plant+
    scale_fill_manual(values = color_names) 
  
  
}

# Plot single model fit with multiple functional groups for total production
# Specify both density and production type

plot_totalprod <- function(year_to_plot = 1995,
                           fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                           model_fit_density = 'pareto',
                           model_fit_production = 'powerlaw',
                           x_limits,
                           x_breaks = c(1, 3, 10, 30, 100,300),
                           y_limits = c(0.03,100),
                           y_breaks = c(0.01,0.1, 1, 10,100, 1000),
                           y_labels,
                           color_names = c("black","#BFE046", "#267038", "#27408b", "#87Cefa", "ivory"),
                           x_name = 'Diameter (cm)',
                           y_name = expression(paste('Total production (kg ha'^-1,' y'^-1,')')),
                           obsdat = obs_totalprod,
                           preddat = fitted_totalprod
) {
  
  require(dplyr)
  require(ggplot2)
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot) %>%
    filter(bin_value > 0)
  
  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(dens_model %in% model_fit_density, prod_model %in% model_fit_production, fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs)
  
  ggplot() +
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    geom_line(data = preddat[preddat$fg == "fg5",], aes(x = dbh, y = q50), color = "gray")+ # white circles get gray line
    geom_ribbon(data = preddat[preddat$fg == "fg5",], aes(x = dbh, ymin = q025, ymax = q975), fill = "gray", alpha = 0.4) +
    geom_point(data = obsdat, aes(x = bin_midpoint, y = bin_value, size=2,group = fg, fill=fg), color = "black",shape=21) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) +
    scale_color_manual(values = color_names) +theme_plant2+
    scale_fill_manual(values = color_names) 
  
  
}  

