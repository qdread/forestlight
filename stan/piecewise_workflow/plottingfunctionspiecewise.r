# Plotting functions for piecewise.

# Define plotting functions -----------------------------------------------

# Plot single model fit with multiple functional groups for density
plot_dens <- function(year_to_plot = 1990,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                      model_fit = 'pareto',
                      x_limits = c(1, 316),
                      x_breaks = c(1, 3, 10, 30, 100),
                      y_limits,
                      y_breaks,
                      color_names = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"
                      ),
                      x_name = 'Diameter (cm)',
                      y_name = 'Density (trees/ha/cm)',
                      obsdat = obs_dens,
                      preddat = pred_dens
) {
  
  require(dplyr)
  require(cowplot)
  
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot) %>%
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
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg), fill = 'gray80') +
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    geom_point(data = obsdat, aes(x = bin_midpoint, y = bin_value, group = fg, color = fg)) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks) +
    scale_color_manual(values = color_names) +
    panel_border(colour = 'black')
  
  
}

# Plot single model fit with multiple functional groups for individual production
plot_prod <- function(year_to_plot = 1990,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                      model_fit = 'powerlaw',
                      x_limits = c(1, 316),
                      x_breaks = c(1, 3, 10, 30, 100),
                      y_limits,
                      y_breaks,
                      color_names = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"
                      ),
                      x_name = 'Diameter (cm)',
                      y_name = 'Production (kg/y)',
                      average = 'mean',
                      error_quantiles = c('ci_min', 'ci_max'),
                      error_bar_width = 0.03,
                      dodge_width = 0.03,
                      dodge_errorbar = TRUE,
                      obsdat = obs_indivprod,
                      preddat = fitted_indivprod
) {
  
  require(dplyr)
  require(cowplot)
  
  pos <- if (dodge_errorbar) position_dodge(width = dodge_width) else 'identity'
  
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot, !is.na(mean), mean > 0) %>%
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
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg), fill = 'gray80') +
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    geom_errorbar(data = obsdat, aes_string(x = 'bin_midpoint', ymin = error_quantiles[1], ymax = error_quantiles[2], group = 'fg', color = 'fg', width = 'width'), position = pos) +
    geom_point(data = obsdat, aes_string(x = 'bin_midpoint', y = average, group = 'fg', color = 'fg'), position = pos) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks) +
    scale_color_manual(values = color_names) +
    panel_border(colour = 'black')
  
  
}

# Plot single model fit with multiple functional groups for total production
# Specify both density and production type
plot_totalprod <- function(year_to_plot = 1990,
                           fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                           model_fit_density = 'pareto',
                           model_fit_production = 'powerlaw',
                           x_limits = c(1, 316),
                           x_breaks = c(1, 3, 10, 30, 100),
                           y_limits,
                           y_breaks,
                           color_names = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"
                           ),
                           x_name = 'Diameter (cm)',
                           y_name = 'Total production (kg/y/ha/cm)',
                           obsdat = obs_totalprod,
                           preddat = fitted_totalprod
) {
  
  require(dplyr)
  require(cowplot)
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
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg), fill = 'gray80') +
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    geom_point(data = obsdat, aes(x = bin_midpoint, y = bin_value, group = fg, color = fg)) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks) +
    scale_color_manual(values = color_names) +
    panel_border(colour = 'black')
  
  
}  
