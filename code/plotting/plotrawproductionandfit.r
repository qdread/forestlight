# Plot raw production data as a hex density plot with the fits on top.
# qdr, forest light, 24 May 2018

# Load data ---------------------------------------------------------------

fp <- 'C:/Users/Q/google_drive/ForestLight/data/data_forplotting_12apr2018' ## CHANGE PATH AS NEEDED

# Read all the csvs in directory.
for (i in dir(fp, pattern = '.csv')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}

# Load raw data
load('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/rawdataobj_22jan.r')

# Process the raw data to get one single data frame with a lot of rows.
library(dplyr)
library(purrr)

# Get only year, func group, dbh, and production (no more is needed to plot right now)
raw_prod <- do.call(rbind, map2(alltreedat, seq(1985,2010,5), function(x, y) cbind(year = y, x %>% select(fg, dbh_corr, production))))

raw_prod <- raw_prod %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg', fg), 'unclassified'))

# Function to plot production with raw data -------------------------------

plot_prod <- function(year_to_plot = 1990,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
                      full_names = c('fast', 'slow', 'pioneer', 'breeder', 'middle', 'unclassified'),
                      func_names = c('power law', 'power law\ntimes exponential'),
                      x_limits = c(1, 316),
                      x_breaks = c(1, 3, 10, 30, 100),
                      y_limits,
                      y_breaks,
                      x_name = 'Diameter (cm)',
                      y_name = 'Production (kg/y)',
                      color_names = c('green', 'blue'),
                      aspect_ratio = 1,
                      hex_scale = scale_fill_gradient(low = 'gray90', high = 'gray10', guide = FALSE),
                      obsdat = raw_prod,
                      preddat = fitted_indivprod
) {
  
  require(dplyr)
  require(cowplot)
  #require(pracma)
  
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot)

  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(production), max_obs = max(production))
  
  # Add more identical rows of the f. groups with less data so that the fill scales appear the same on all plots.
  #multiples <- Reduce(Lcm, table(obsdat$fg)) / table(obsdat$fg)
  #obsdat <- obsdat[rep(1:nrow(obsdat), ceiling(multiples/(min(multiples)/3))[obsdat$fg]), ]
  
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(dens_model == 'pareto', fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs) %>%
    mutate(prod_model = factor(prod_model, labels = func_names))
  
  labels <- setNames(full_names, fg_names)
  
  ggplot() +
    geom_hex(data = obsdat, aes(x = dbh_corr, y = production)) +
    geom_line(data = preddat, aes(x = dbh, y = q50, group = prod_model, color = prod_model)) +
    geom_line(data = preddat, aes(x = dbh, y = q025, group = prod_model, color = prod_model), size = 0.2, linetype = 'dotted') +
    geom_line(data = preddat, aes(x = dbh, y = q975, group = prod_model, color = prod_model), size = 0.2, linetype = 'dotted') +
    facet_wrap(~ fg, labeller = labeller(fg = labels)) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks) +
    scale_color_manual(values = color_names, name = 'Functional form') +
    hex_scale +
    panel_border(colour = 'black') +
    coord_fixed(ratio = aspect_ratio) +
    theme(legend.position = c(0.8, 0.2), strip.background = element_blank())
  
  
}


# Call function to draw plot ----------------------------------------------

# Original grayscale
hex_scale_1 <- scale_fill_gradient(low = 'gray90', high = 'gray10', guide = FALSE)

# Biased grayscale to emphasize the variation in the hexagons with less data.
# Bias >1 does this.
hex_scale_2 <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=10)(50))

# Biased color scale of red yellow and blue to do the same, using a brewer palette
# the bias is <1 this time because I had to reverse the scale.
hex_scale_3 <- scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(9,'RdYlBu'), bias=0.3)(50)), guide = FALSE)

# Biased and customized color scale
hex_scale_4 <- scale_fill_gradientn(colours = colorRampPalette(c('skyblue', 'goldenrod', 'indianred'), bias = 3)(50), guide = FALSE)

# Edit the hex scale argument to draw this with other color scales.
hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=1)(50), trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))

plot_prod(year_to_plot = 1990,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          full_names = c('fast', 'slow', 'pioneer', 'breeder', 'middle', 'unclassified'),
          x_limits = c(1, 316),
          x_breaks = c(1, 3, 10, 30, 100),
          y_limits = c(5e-03, 1e04),
          y_breaks = c(.001, .1, 10, 1000),
          color_names = c('red', 'purple'),
          hex_scale = hex_scale_log,
          aspect_ratio = 0.4) 
          
          
