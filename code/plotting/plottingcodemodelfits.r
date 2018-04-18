# Plot observed data, predicted values of different model fits, and different confidence intervals on model fits
# Edit 03 April: add dodge error bars and a dashed line at various places.
# Edit 11 April: predictions should not extend past bins

# Load data ---------------------------------------------------------------

fp <- 'C:/Users/Q/google_drive/ForestLight/data/data_forplotting_12apr2018' ## CHANGE PATH AS NEEDED

# Read all the csvs in directory.
for (i in dir(fp, pattern = '.csv')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}

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
    filter(dens_model %in% model_fit, prod_model == 'powerlaw', fg %in% fg_names, year == year_to_plot) %>%
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
                      preddat = pred_indivprod
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
    filter(prod_model %in% model_fit, dens_model == 'pareto', fg %in% fg_names, year == year_to_plot) %>%
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
                      preddat = pred_totalprod
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

# Example plot ------------------------------------------------------------

# Density plots for all functional groups with Pareto fit in 1995
plot_dens(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = 'pareto',
          y_limits = c(0.001, 1000),
          y_breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000))

# Density plots for all functional groups with Weibull fit in 1995
plot_dens(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = 'weibull',
          y_limits = c(0.0001, 1000),
          y_breaks = c(0.001, 0.1, 10, 1000))

# Density plot for just one fg with Weibull fit in 1995
plot_dens(year_to_plot = 1995,
          fg_names = c('fg1'),
          model_fit = 'weibull',
          y_limits = c(0.0001, 10),
          y_breaks = c(0.001, 0.1, 10))

# Production plots for all functional groups with power law fit in 1995
plot_prod(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = 'powerlaw',
          y_limits = c(0.01, 6000),
          y_breaks = c(0.1, 10, 1000))

# Production plots for all functional groups with power law fit in 1995
# Specify not to dodge
plot_prod(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = 'powerlaw',
          y_limits = c(0.01, 6000),
          y_breaks = c(0.1, 10, 1000),
          dodge_errorbar = FALSE)

# Specify dodging with a certain width of error bar
plot_prod(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = 'powerlaw',
          y_limits = c(0.01, 6000),
          y_breaks = c(0.1, 10, 1000),
          error_bar_width = 0.01,
          dodge_width = 0.05)

# Production plot for just one fg with power law times expo fit in 1995
# Specify that we are using the geometric mean and ci around it
plot_prod(year_to_plot = 1995,
          fg_names = c('fg1'),
          model_fit = 'powerlawexp',
          y_limits = c(0.01, 7000),
          y_breaks = c(0.1, 10, 1000),
          error_quantiles = c('ci_min', 'ci_max'),
          average = 'mean',
          dodge_errorbar = FALSE,
          color_names = 'black')

# Specify that we are using the median and the .025 and .975 quantiles
plot_prod(year_to_plot = 1995,
          fg_names = c('fg1'),
          model_fit = 'powerlawexp',
          y_limits = c(0.01, 7000),
          y_breaks = c(0.1, 10, 1000),
          error_quantiles = c('q025', 'q975'),
          average = 'median',
          color_names = 'black')

plot_prod(year_to_plot = 1995,
          fg_names = c('fg1'),
          model_fit = 'powerlawexp',
          y_limits = c(0.01, 7000),
          y_breaks = c(0.1, 10, 1000),
          error_quantiles = c('ci_min', 'ci_max'),
          average = 'mean',
          color_names = 'black')


# Total production plot for some functional groups for Weibull and power law in 1995
plot_totalprod(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3'),
          model_fit_density = 'weibull', 
          model_fit_production = 'powerlaw',
          y_limits = c(0.01, 200),
          y_breaks = c(0.1, 1, 10, 100))

plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg4','fg5'),
               model_fit_density = 'weibull', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 500),
               y_breaks = c(0.1, 1, 10, 100))

plot_totalprod(year_to_plot = 1995,
               fg_names = c('all'),
               model_fit_density = 'weibull', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 500),
               y_breaks = c(0.1, 1, 10, 100))

# Add reference line to a plot
# Density plots for all functional groups with Weibull fit in 1995
p <- plot_dens(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = 'weibull',
          y_limits = c(0.0001, 1000),
          y_breaks = c(0.001, 0.1, 10, 1000))

p + geom_abline(intercept = 4, slope = -2, color = 'darkgray', linetype = 'dotted', size = 0.5)
