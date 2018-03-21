# Plot observed data, predicted values of different model fits, and different confidence intervals on model fits


# Load data ---------------------------------------------------------------

fp <- 'C:/Users/Q/google_drive/ForestLight/data/data_20mar2018' ## CHANGE PATH AS NEEDED

# Read all the csvs in directory.
for (i in dir(fp)) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp,i), stringsAsFactors = FALSE))
}

# Define plotting functions -----------------------------------------------

plot_dens <- function(year_to_plot = 1990,
                           fg_names = c('fg1','fg2','fg3','fg4','fg5','all')
                           model_fit = c('pareto', 'weibull'),
                           x_limits = c(1, 316),
                           x_breaks,
                           y_limits,
                           y_breaks,
                           obsdat = obs_dens,
                           preddat = pred_dens
) {
  
  require(dplyr)
  require(cowplot)
  
  preddat <- preddat %>%
    filter(dens_model %in% model_fit, prod_model == 'powerlaw', fg %in% fg_names, year == year_to_plot)
  
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot)
  
  ggplot() +
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975), fill = 'gray50') +
    geom_line(data = preddat, aes(x = dbh, y = q50)) +
    geom_point(data = obsdat, aes(x = bin_midpoint, y = bin_value)) +
    
  
}
  
