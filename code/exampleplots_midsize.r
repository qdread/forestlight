# Example plotting code for the mid size trees

# Density plots for all functional groups with Pareto fit in 1995
plot_dens(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = 'pareto',
          y_limits = c(0.001, 1000),
          y_breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
          obsdat = obs_dens_midsize, preddat = pred_dens_midsize)

# Density plots for all functional groups with Weibull fit in 1995
plot_dens(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = 'weibull',
          y_limits = c(0.0001, 1000),
          y_breaks = c(0.001, 0.1, 10, 1000),
          obsdat = obs_dens_midsize, preddat = pred_dens_midsize)


# Production plots for all functional groups with power law fit in 1995
plot_prod(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = 'powerlaw',
          y_limits = c(0.01, 6000),
          y_breaks = c(0.1, 10, 1000),
          obsdat = obs_indivprod_midsize, preddat = pred_indivprod_midsize)


# Total production plot for some functional groups for Weibull/Pareto and power law in 1995

plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3'),
               model_fit_density = 'pareto', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 200),
               y_breaks = c(0.1, 1, 10, 100),
               obsdat = obs_totalprod_midsize, preddat = pred_totalprod_midsize)

plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3'),
               model_fit_density = 'weibull', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 200),
               y_breaks = c(0.1, 1, 10, 100),
               obsdat = obs_totalprod_midsize, preddat = pred_totalprod_midsize)

plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg4','fg5'),
               model_fit_density = 'pareto', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 500),
               y_breaks = c(0.1, 1, 10, 100),
               obsdat = obs_totalprod_midsize, preddat = pred_totalprod_midsize)

plot_totalprod(year_to_plot = 1995,
               fg_names = c('all'),
               model_fit_density = 'pareto', 
               model_fit_production = 'powerlaw',
               y_limits = c(10, 500),
               y_breaks = c(10, 100),
               obsdat = obs_totalprod_midsize, preddat = pred_totalprod_midsize)

plot_totalprod(year_to_plot = 1995,
               fg_names = c('all'),
               model_fit_density = 'weibull', 
               model_fit_production = 'powerlaw',
               y_limits = c(10, 500),
               y_breaks = c(10, 100),
               obsdat = obs_totalprod_midsize, preddat = pred_totalprod_midsize)

# Add reference line to a plot
# Density plots for all functional groups with Weibull fit in 1995
p <- plot_dens(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
               model_fit = 'weibull',
               y_limits = c(0.0001, 1000),
               y_breaks = c(0.001, 0.1, 10, 1000))

p + geom_abline(intercept = 4, slope = -2, color = 'darkgray', linetype = 'dotted', size = 0.5)
