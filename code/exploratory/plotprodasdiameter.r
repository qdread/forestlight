p <- plot_prod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3','fg4','fg5'),
               model_fit = 2,
               x_limits = c(1, 280),
               y_limits = c(0.02, 1),
               y_breaks = c(0.03, 0.1, 0.3, 1),
               y_labels = c(0.03, 0.1, 0.3, 1),
               error_bar_width = 0.01,
               dodge_width = 0.05,
               obsdat = obs_indivdiamgrowth,
               preddat = fitted_indivdiamgrowth,
               plot_abline = FALSE,
               x_name = 'Diameter (cm)',
               y_name = expression(paste('Diameter growth (cm y'^-1,')')))

p + theme(axis.text.x = element_text()) + labs(x = 'Diameter (cm)')
