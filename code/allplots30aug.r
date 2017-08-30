# All plots and data compilation, after running all the scalings

load('C:/Users/Q/Dropbox/projects/forestlight/allscalings_30aug.RData') # Very large. Takes a long time to load. (400mb on hard disk)

fpfig <- 'C:/Users/Q/google_drive/ForestLight/figs/figures_30aug/'

source('code/allplottingfns30aug.r')

library(cowplot)

area_core <- 42.84 # area in hectares of the BCI plot after subtracting secondary forest and the 20-m wide strip around the edge.

# Descriptive names
longnames <- c('All trees', 'Shade-tolerant', 'Intermediate', 'Gap', 'Unclassified', 'Shade-tolerant, quantile', 'Intermediate, quantile', 'Gap, quantile')
# Better order for plotting
plot_order <- c(1:4, 6:8, 5)
# Census years
years <- c(1985, 1990, 1995, 2000, 2005, 2010)

# 1. Density by diameter scaling.

## Generate bins

numbins <- 20 # for visualization

for (i in allyears_names) {
  assign(paste0(i, '_dens_logbin'), lapply(get(i), function(z) logbin(x=z$dbh_corr, y=NULL, n=numbins)))  
}

for (i in c(names1990, names1995)) {
  assign(paste0(i, '_dens_logbin'), logbin(x=get(i)$dbh_corr, y=NULL, n=numbins)) 
}

## Plot

xl2 <- 'Diameter (cm)'
yl2 <- 'Density'

for (i in allyears_names) {
  dens_fit_i <- get(paste0(i, '_dens_fit'))
  logbin_i <- get(paste0(i, '_dens_logbin'))
  dens_plots_i <- list()
  for (j in 1:length(dens_fit_i)) {
    dens_plots_i[[j]] <- plotbinsandfits(dens_fit_i[[j]], logbin_i[[j]], plottitle = paste(i,'at census',j), xl = xl2, yl = yl2, plottype = 'point')
  }
  assign(paste0(i, '_dens_plots'), dens_plots_i)
}

for (i in c(names1990, names1995)) {
  dens_fit_i <- get(paste0(i, '_dens_fit'))
  logbin_i <- get(paste0(i, '_dens_logbin'))
  assign(paste0(i, '_dens_plots'), plotbinsandfits(dens_fit_i, logbin_i, plottitle = i, xl = xl2, yl = yl2, plottype = 'point'))
}

# 1b. Density by diameter scaling, include the Pareto fit and Pareto-cutoff fit
# Also make sure that the axis limits are the same in every panel so that you can compare the different scalings.

dens_lim <- get_axis_limits(paste0(allyears_names, '_dens_logbin'))

xl2 <- 'Diameter (cm)'
yl2 <- expression(paste('Density (individuals ha'^-1,')', sep=''))

for (i in 1:length(allyears_names)) {
  name_i <- allyears_names[i]
  dens_fit_i <- get(paste0(name_i, '_dens_fit'))
  cutoff_fit_i <- get(paste0(name_i, '_paretofits'))
  logbin_i <- get(paste0(name_i, '_dens_logbin'))
  dens_plots_i <- list()
  for (j in 1:length(dens_fit_i)) {
    dens_plots_i[[j]] <- plotbinsandfits_cutoff(dens_fit_i[[j]], cutoff_fit_i[[j]]$fit_cutoff, logbin_i[[j]], plottitle = paste(longnames[i], 'density,', years[j]), xl = xl2, yl = yl2, plottype = 'point', plotarea = area_core, y_min = 0.00001, y_max = dens_lim$y_max/area_core, x_max = dens_lim$x_max, y_values = 10^c(-3,-1,1,3))
  }
  assign(paste0(name_i, '_dens_plots'), dens_plots_i)
}

# Create grid
plot_columns <- list()

for (i in 1:length(allyears_names)) {
  plot_columns[[i]] <- plot_grid(plotlist = get(paste0(allyears_names[plot_order[i]], '_dens_plots')), ncol = 1)
}

dens_grid <- plot_grid(plotlist = plot_columns, ncol = length(plot_columns))
ggsave(file.path(fpfig, 'density_scalings_allyears.png'), dens_grid, height = 30, width = 50, dpi = 300, limitsize = FALSE)

# Create each year as a figure.

plot_years <- list()

for (i in 1:length(years)) {
  plot_years[[i]] <- list()
  for (j in 1:length(allyears_names)) {
    plot_years[[i]][[j]] <- get(paste0(allyears_names[plot_order[j]], '_dens_plots'))[[i]]
  }
  
  grid_year_i <- plot_grid(plotlist = plot_years[[i]], nrow = 2)
  ggsave(file.path(fpfig, paste0('density_scalings_',years[i],'.png')), grid_year_i, height = 8, width = 16, dpi = 300)
}


