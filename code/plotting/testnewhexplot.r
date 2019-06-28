## New Hex Plot Code
gdrive_path <- '~/google_drive/ForestLight'

# Load raw data
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

# Load fitted data
fitted_indivprod <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/fitted_indivprod.csv'), stringsAsFactors = FALSE)

# Process the raw data to get one single data frame with a lot of rows.
library(dplyr)
library(purrr)
library(ggplot2)

theme_plant <- theme(panel.grid = element_blank(), #for Total Production
                     aspect.ratio = .75,
                     axis.text = element_text(size = 19, color = "black"), 
                     axis.ticks.length=unit(0.2,"cm"),
                     axis.title = element_text(size = 19),
                     axis.title.y = element_text(margin = margin(r = 10)),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.x.top = element_text(margin = margin(b = 5)),
                     plot.title = element_text(size = 19, face = "plain", hjust = 10),
                     panel.border = element_rect(color = "black", fill=NA,  size=1),
                     panel.background = element_blank(),
                     legend.position = "none",
                     rect = element_rect(fill = "transparent"),
                     text = element_text(family = 'Helvetica')) 

# Get only year, func group, dbh, and production (no more is needed to plot right now)
raw_prod <- do.call(rbind, map2(alltreedat, seq(1985,2010,5), function(x, y) cbind(year = y, x %>% select(fg, dbh_corr, production))))

raw_prod <- raw_prod %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg', fg), 'unclassified'))

# Function to plot production with raw data -------

plot_rawprod <- function(year_to_plot = 1995,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
                      full_names = c('Fast', 'LL Pioneer', 'Slow', 'SL Breeder', 'Medium', 'Unclassified'),
                      func_names = c('1 segment power law', '2 segment power law'),
                      x_limits = c(1, 300),
                      x_breaks = c(1, 10, 100),
                      y_limits,
                      y_breaks,
                      x_name = 'Diameter (cm)',
                      y_name = expression(paste('Growth (kg yr'^-1,')')), line_types = c('dashed', 'solid'),
                      aspect_ratio = 0.75,
                      hex_scale = hex_scale_log_colors,
                      obsdat = raw_prod,
                      preddat = fitted_indivprod
) {
  
  require(dplyr)
  require(ggplot2)

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
    filter(fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs) %>%
    mutate(prod_model = factor(prod_model, labels = func_names))
  
  labels <- setNames(full_names, fg_names)
  
  ggplot() +
    geom_hex(data = obsdat, aes(x = dbh_corr, y = production)) +
    #geom_line(data = preddat, aes(x = dbh, y = q50, group = prod_model, linetype = prod_model), size=0.25) +
    geom_abline(slope = 2, intercept = -2.1, linetype = "dashed")+
    facet_wrap(~ fg, ncol = 2, labeller = labeller(fg = labels)) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels=signif) +
    scale_linetype_manual(values = line_types) +
    hex_scale + theme_plant+
    coord_fixed(ratio = aspect_ratio) + guides(linetype = 'none')+
    theme(legend.position = c(0.7, 0.15), strip.background = element_blank(),
          strip.text = element_text(size=12))#, legend.text = element_blank())
  
  
}
# Original grayscale
hex_scale_1 <- scale_fill_gradient(low = 'gray90', high = 'gray10', guide = FALSE)

# Biased grayscale to emphasize the variation in the hexagons with less data.
# Bias >1 does this.
#hex_scale_2 <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=10)(50))

# Biased color scale of red yellow and blue to do the same, using a brewer palette
# the bias is <1 this time because I had to reverse the scale.
#hex_scale_3 <- scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(9,'RdYlBu'), bias=0.3)(50)), guide = FALSE)

# Biased and customized color scale
#hex_scale_4 <- scale_fill_gradientn(colours = colorRampPalette(c('forestgreen', 'goldenrod', 'indianred'), bias = 3)(50), guide = FALSE)

# Edit the hex scale argument to draw this with other color scales.
#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=1)(50), trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))
hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000,10000), 
                                             labels = c(1,10,100,1000,10000), limits=c(1,10000))


p <- plot_rawprod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
               full_names = c('Fast', 'Pioneer', 'Slow', 'Breeder', 'Medium', 'Unclassified'),
               x_limits = c(1, 316),
               x_breaks = c(1,10, 100),
               y_limits = c(0.001, 1000),
               y_breaks = c(.001, .1, 10, 1000),
               line_types = c('dashed', 'solid'),
               hex_scale = hex_scale_log_colors)
p 
pdf(file.path(gdrive_path, 'Figures/Growth_hex/growth_hex.pdf'))
p
dev.off()
#pdf(file.path(gdrive_path, "Plots_J/New/Supplementals/raw_growth_heat.pdf"))
#p
#dev.off()  