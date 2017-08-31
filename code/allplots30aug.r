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

# 1c. Density scaling by diameter, only trees that have light received values for 1990 and 1995.
# With and without cutoffs.

dens_lim90 <- get_axis_limits(paste0(names1990, '_dens_logbin'))

for (i in 1:length(names1990)) {
  name_i <- names1990[i]
  dens_fit_i <- get(paste0(name_i, '_dens_fit'))
  cutoff_fit_i <- get(paste0(name_i, '_paretofits'))
  logbin_i <- get(paste0(name_i, '_dens_logbin'))
  dens_plots_i <- plotbinsandfits_cutoff(dens_fit_i, cutoff_fit_i$fit_cutoff, logbin_i, plottitle = paste(longnames[i], 'density, 1990'), xl = xl2, yl = yl2, plottype = 'point', plotarea = area_core, y_min = dens_lim90$y_min/area_core, y_max = dens_lim90$y_max/area_core, x_max = dens_lim90$x_max, y_values = 10^c(-3,-1,1,3))
  assign(paste0(name_i, '_dens_plots'), dens_plots_i)
}

dens_lim95 <- get_axis_limits(paste0(names1995, '_dens_logbin'))

for (i in 1:length(names1995)) {
  name_i <- names1995[i]
  dens_fit_i <- get(paste0(name_i, '_dens_fit'))
  cutoff_fit_i <- get(paste0(name_i, '_paretofits'))
  logbin_i <- get(paste0(name_i, '_dens_logbin'))
  dens_plots_i <- plotbinsandfits_cutoff(dens_fit_i, cutoff_fit_i$fit_cutoff, logbin_i, plottitle = paste(longnames[i], 'density, 1995'), xl = xl2, yl = yl2, plottype = 'point', plotarea = area_core, y_min = dens_lim95$y_min/area_core, y_max = dens_lim95$y_max/area_core, x_max = dens_lim95$x_max, y_values = 10^c(-2, 0, 2))
  assign(paste0(name_i, '_dens_plots'), dens_plots_i)
}

# 2. Density scaled by amount of light received per individual, in 1990 and 1995 only
# Plot with and without cutoffs.

## Generate bins

numbins <- 20 # for visualization

for (i in c(names1990, names1995)) {
  assign(paste0(i, '_lightdens_logbin'), logbin(x=get(i)$light_received, y=NULL, n=numbins)) 
}

xl4 <- 'Light received (W)'
yl4 <- expression(paste('Density (individuals ha'^-1,')', sep=''))

lightdens_lim90 <- get_axis_limits(paste0(names1990, '_lightdens_logbin'))

for (i in 1:length(names1990)) {
  name_i <- names1990[i]
  dens_fit_i <- get(paste0(name_i, '_lightdens_fit'))
  cutoff_fit_i <- get(paste0(name_i, '_lightparetofits'))
  logbin_i <- get(paste0(name_i, '_lightdens_logbin'))
  dens_plots_i <- plotbinsandfits_cutoff(dens_fit_i, cutoff_fit_i$fit_cutoff, logbin_i, plottitle = paste(longnames[i], 'density, 1990'), xl = xl4, yl = yl4, plottype = 'point', plotarea = area_core, y_min = lightdens_lim90$y_min/area_core, y_max = lightdens_lim90$y_max/area_core, x_max = lightdens_lim90$x_max, x_values = 10^c(0:6), y_values = 10^c(-4, -2, 0, 2))
  assign(paste0(name_i, '_lightdens_plots'), dens_plots_i)
}

lightdens_lim95 <- get_axis_limits(paste0(names1995, '_lightdens_logbin'))

for (i in 1:length(names1995)) {
  name_i <- names1995[i]
  dens_fit_i <- get(paste0(name_i, '_lightdens_fit'))
  cutoff_fit_i <- get(paste0(name_i, '_lightparetofits'))
  logbin_i <- get(paste0(name_i, '_lightdens_logbin'))
  dens_plots_i <- plotbinsandfits_cutoff(dens_fit_i, cutoff_fit_i$fit_cutoff, logbin_i, plottitle = paste(longnames[i], 'density, 1995'), xl = xl4, yl = yl4, plottype = 'point', plotarea = area_core, y_min = lightdens_lim95$y_min/area_core, y_max = lightdens_lim95$y_max/area_core, x_max = lightdens_lim95$x_max, x_values = 10^c(0:6), y_values = 10^c(-4, -2, 0, 2))
  assign(paste0(name_i, '_lightdens_plots'), dens_plots_i)
}

# 3. Individual production by diameter
# try it out using something like hexbin, or points with error bars on them.
# Plot both the cutoffs and no cutoffs.

xl1 <- 'Diameter (cm)'
yl1 <- expression(paste('Individual production (kg y'^-1,')', sep=''))
ptype <- 'hex'

# 1990
# Only trees with light measurements. Use light_received!
 twoslopeplot(dat = alltree_light_90, plottitle = 'All species 1990', xl = xl1, yl =  yl1, plottype = ptype)

 twoslopeplot(dat = shade_light_90, plottitle = 'Shade-tolerant species 1990', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = int_light_90, plottitle = 'Intermediate species 1990', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = gap_light_90, plottitle = 'Gap species 1990', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = shadequant_light_90, plottitle = 'Shade-tolerant species 1990 (by quantile)', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = intquant_light_90, plottitle = 'Intermediate species 1990 (by quantile)', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = gapquant_light_90, plottitle = 'Gap species 1990 (by quantile)', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = unclassified_light_90, plottitle = 'Unclassified species 1990', xl = xl1, yl =  yl1, plottype = ptype)
#
# # 1995
 twoslopeplot(dat = alltree_light_95, plottitle = 'All species 1995', xl = xl1, yl =  yl1, plottype = ptype)
 
 twoslopeplot(dat = shade_light_95, plottitle = 'Shade-tolerant species 1995', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = int_light_95, plottitle = 'Intermediate species 1995', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = gap_light_95, plottitle = 'Gap species 1995', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = shadequant_light_95, plottitle = 'Shade-tolerant species 1995 (by quantile)', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = intquant_light_95, plottitle = 'Intermediate species 1995 (by quantile)', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = gapquant_light_95, plottitle = 'Gap species 1995 (by quantile)', xl = xl1, yl =  yl1, plottype = ptype)
 twoslopeplot(dat = unclassified_light_95, plottitle = 'Unclassified species 1995', xl = xl1, yl =  yl1, plottype = ptype)

 # 4. Individual production by light received
 
 hex_plot <- function(dat, xv, yv, plottitle, xl, yl, textx, texty) {
   fit <- lm(formula = paste0('I(log10(',yv,'))~I(log10(',xv,'))'), data = dat)
   fit_r2 <- round(summary(fit)$r.sq, 2)
   fit_slope <- round(fit$coeff[2], 3)
   label_df <- data.frame(x = textx, y = texty, lab = paste0('slope = ', fit_slope, ', r2 = ', fit_r2))
   names(label_df) <- c(xv, yv, 'lab')
   
   p <- ggplot(dat, aes_string(x = xv, y = yv))
   p <- p + geom_hex() + 
     stat_smooth(method = 'lm', se = FALSE) +
     geom_text(data = label_df, aes(label=lab)) +
     scale_fill_gradientn(colours = colorRampPalette(RColorBrewer::brewer.pal(9, 'YlOrRd'), bias=2)(10)) +
     theme(legend.position = c(0.8, 0.2))
   
   
   p <- p +
     scale_x_log10(name = xl,
                   breaks = scales::trans_breaks("log10", function(x) 10^x),
                   labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
     scale_y_log10(name = yl,
                   breaks = scales::trans_breaks("log10", function(x) 10^x),
                   labels = scales::trans_format("log10", scales::math_format(10^.x))) +
     ggtitle(plottitle) +
     panel_border(colour = 'black')
   return(p)
 }

 scatter_plot <- function(dat, xv, yv, plottitle, xl, yl, textx, texty) {
   fit <- lm(formula = paste0('I(log10(',yv,'))~I(log10(',xv,'))'), data = dat)
   fit_r2 <- round(summary(fit)$r.sq, 2)
   fit_slope <- round(fit$coeff[2], 3)
   label_df <- data.frame(x = textx, y = texty, lab = paste0('slope = ', fit_slope, ', r2 = ', fit_r2))
   names(label_df) <- c(xv, yv, 'lab')
   
   p <- ggplot(dat, aes_string(x = xv, y = yv))
   p <- p + 
     geom_point(alpha = 0.6) +
     stat_smooth(method = 'lm', se = FALSE) +
     geom_text(data = label_df, aes(label=lab))
   
   
   p <- p +
     scale_x_log10(name = xl,
                   breaks = scales::trans_breaks("log10", function(x) 10^x),
                   labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
     scale_y_log10(name = yl,
                   breaks = scales::trans_breaks("log10", function(x) 10^x),
                   labels = scales::trans_format("log10", scales::math_format(10^.x))) +
     ggtitle(plottitle) +
     panel_border(colour = 'black')
   return(p)
 }
 
 xl1b <- 'Light received (W)'
 yl1b <- expression(paste('Individual production (kg y'^-1,')', sep='')) 
 
 prod_light_hexplots90 <- list()
 prod_light_scatterplots90 <- list()
 prod_light_hexplots95 <- list()
 prod_light_scatterplots95 <- list()
 
 for (i in 1:length(names1990)) {
   prod_light_hexplots90[[length(prod_light_hexplots90) + 1]] <- hex_plot(dat = get(names1990[plot_order[i]]),  xv = 'light_received', yv = 'production', plottitle = paste(longnames[plot_order[i]], '1990'), xl = xl1b, yl =  yl1b, textx = 10, texty = 10000)
   prod_light_scatterplots90[[length(prod_light_scatterplots90) + 1]] <- scatter_plot(dat = get(names1990[plot_order[i]]),  xv = 'light_received', yv = 'production', plottitle = paste(longnames[plot_order[i]], '1990'), xl = xl1b, yl =  yl1b, textx = 10, texty = 10000)
 }

 for (i in 1:length(names1995)) {
   prod_light_hexplots95[[length(prod_light_hexplots95) + 1]] <- hex_plot(dat = get(names1990[plot_order[i]]),  xv = 'light_received', yv = 'production', plottitle = paste(longnames[plot_order[i]], '1995'), xl = xl1b, yl =  yl1b, textx = 10, texty = 10000)
   prod_light_scatterplots95[[length(prod_light_scatterplots95) + 1]] <- scatter_plot(dat = get(names1990[plot_order[i]]),  xv = 'light_received', yv = 'production', plottitle = paste(longnames[plot_order[i]], '1995'), xl = xl1b, yl =  yl1b, textx = 10, texty = 10000)
 }
 
 # 5. Individual light received by diameter
 
 xl1c <- 'Diameter (cm)'
 yl1c <- 'Light received (W)'
 
 light_diam_hexplots90 <- list()
 light_diam_scatterplots90 <- list()
 light_diam_hexplots95 <- list()
 light_diam_scatterplots95 <- list()
 
 for (i in 1:length(names1990)) {
   light_diam_hexplots90[[length(light_diam_hexplots90) + 1]] <- hex_plot(dat = get(names1990[plot_order[i]]),  xv = 'dbh_corr', yv = 'light_received', plottitle = paste(longnames[plot_order[i]], '1990'), xl = xl1b, yl =  yl1b, textx = 5, texty = 100000)
   light_diam_scatterplots90[[length(light_diam_scatterplots90) + 1]] <- scatter_plot(dat = get(names1990[plot_order[i]]),  xv = 'dbh_corr', yv = 'light_received', plottitle = paste(longnames[plot_order[i]], '1990'), xl = xl1b, yl =  yl1b, textx = 5, texty = 100000)
 }
 
 for (i in 1:length(names1995)) {
   light_diam_hexplots95[[length(light_diam_hexplots95) + 1]] <- hex_plot(dat = get(names1990[plot_order[i]]),  xv = 'dbh_corr', yv = 'light_received', plottitle = paste(longnames[plot_order[i]], '1995'), xl = xl1b, yl =  yl1b, textx = 5, texty = 100000)
   light_diam_scatterplots95[[length(light_diam_scatterplots95) + 1]] <- scatter_plot(dat = get(names1990[plot_order[i]]),  xv = 'dbh_corr', yv = 'light_received', plottitle = paste(longnames[plot_order[i]], '1995'), xl = xl1b, yl =  yl1b, textx = 5, texty = 100000)
 }
 
 # 6. Binned production by diameter
 
 xl3 <- 'Diameter (cm)'
 yl3 <- expression(paste('Total production (kg y'^-1,' ha'^-1, ')', sep=''))
 
prod_lim <- get_axis_limits(paste0(allyears_names, '_prod_logbin'))
 
 for (i in 1:length(allyears_names)) {
   name_i <- allyears_names[i]
   cutoff_fit_i <- get(paste0(name_i, '_paretofits'))
   logbin_i <- get(paste0(name_i, '_prod_logbin'))
   prodbin_plots_i <- list()
   for (j in 1:length(logbin_i)) {
     prodbin_plots_i[[j]] <- plotlogbin_cutoff(logbin_i[[j]], cutoff=cutoff_fit_i[[j]]$fit_cutoff@coef[2], reg=TRUE, plottitle = paste(longnames[i], 'total production,', years[j]), xl = xl3, yl = yl3, plottype = 'point', plotarea = area_core, y_min = 0.001, y_max = prod_lim$y_max/area_core, x_max = prod_lim$x_max, y_values = 10^c(-2,0,2))
   }
   assign(paste0(name_i, '_prodbin_plots'), prodbin_plots_i)
 }
 
 # Create grid
 plot_columns <- list()
 
 for (i in 1:length(allyears_names)) {
   plot_columns[[i]] <- plot_grid(plotlist = get(paste0(allyears_names[plot_order[i]], '_prodbin_plots')), ncol = 1)
 }
 
prod_grid <- plot_grid(plotlist = plot_columns, ncol = length(plot_columns))
 ggsave(file.path(fpfig, 'production_scalings_allyears.png'), prod_grid, height = 30, width = 50, dpi = 300, limitsize = FALSE)
 
 # Create each year as a figure.
 
 plot_years <- list()
 
 for (i in 1:length(years)) {
   plot_years[[i]] <- list()
   for (j in 1:length(allyears_names)) {
     plot_years[[i]][[j]] <- get(paste0(allyears_names[plot_order[j]], '_prodbin_plots'))[[i]]
   }
   
   grid_year_i <- plot_grid(plotlist = plot_years[[i]], nrow = 2)
   ggsave(file.path(fpfig, paste0('production_scalings_',years[i],'.png')), grid_year_i, height = 8, width = 16, dpi = 300)
 }
 
 # 6b. Binned production scaling by diameter, only trees that have light received values for 1990 and 1995.
 # With and without cutoffs.
 
 prod_lim90 <- get_axis_limits(paste0(names1990, '_prod_logbin'))
 
 for (i in 1:length(names1990)) {
   name_i <- names1990[i]
   cutoff_fit_i <- get(paste0(name_i, '_paretofits'))
   logbin_i <- get(paste0(name_i, '_prod_logbin'))
   prod_plots_i <- plotlogbin_cutoff(logbin_i, cutoff=cutoff_fit_i$fit_cutoff@coef[2], reg=TRUE, plottitle = paste(longnames[i], 'total production, 1990'), xl = xl3, yl = yl3, plottype = 'point', plotarea = area_core, y_min = prod_lim90$y_min/area_core, y_max = prod_lim90$y_max/area_core, x_max = prod_lim90$x_max, y_values = 10^c(-2,0,2))
   assign(paste0(name_i, '_prodbin_plots'), prod_plots_i)
 }
 
 prod_lim95 <- get_axis_limits(paste0(names1995, '_prod_logbin'))
 
 for (i in 1:length(names1995)) {
   name_i <- names1995[i]
   cutoff_fit_i <- get(paste0(name_i, '_paretofits'))
   logbin_i <- get(paste0(name_i, '_prod_logbin'))
   prod_plots_i <- plotlogbin_cutoff(logbin_i, cutoff=cutoff_fit_i$fit_cutoff@coef[2], reg=TRUE, plottitle = paste(longnames[i], 'total production, 1995'), xl = xl3, yl = yl3, plottype = 'point', plotarea = area_core, y_min = prod_lim95$y_min/area_core, y_max = prod_lim95$y_max/area_core, x_max = prod_lim95$x_max, y_values = 10^c(-2,0,2))
   assign(paste0(name_i, '_prodbin_plots'), prod_plots_i)
 }
 
 # 7. Binned light scaling by diameter, only trees that have light received data for '90 and '95
 # With and without cutoffs
 
 xl4 <- 'Diameter (cm)'
 yl4 <-expression(paste("Total energy received (W ha"^-1, ")", sep = ""))
 
 
 light_lim90 <- get_axis_limits(paste0(names1990, '_light_logbin'))
 
 for (i in 1:length(names1990)) {
   name_i <- names1990[i]
   cutoff_fit_i <- get(paste0(name_i, '_lightparetofits'))
   logbin_i <- get(paste0(name_i, '_light_logbin'))
   light_plots_i <- plotlogbin_cutoff(logbin_i, cutoff=cutoff_fit_i$fit_cutoff@coef[2], reg=TRUE, plottitle = paste(longnames[i], 'total light received, 1990'), xl = xl4, yl = yl4, plottype = 'point', plotarea = area_core, y_min = light_lim90$y_min/area_core, y_max = light_lim90$y_max/area_core, x_max = light_lim90$x_max, y_values = 10^c(0,2,4,6))
   assign(paste0(name_i, '_lightbin_plots'), light_plots_i)
 }
 
 light_lim95 <- get_axis_limits(paste0(names1995, '_light_logbin'))
 
 for (i in 1:length(names1995)) {
   name_i <- names1995[i]
   cutoff_fit_i <- get(paste0(name_i, '_lightparetofits'))
   logbin_i <- get(paste0(name_i, '_light_logbin'))
   light_plots_i <- plotlogbin_cutoff(logbin_i, cutoff=cutoff_fit_i$fit_cutoff@coef[2], reg=TRUE, plottitle = paste(longnames[i], 'total light received, 1995'), xl = xl4, yl = yl4, plottype = 'point', plotarea = area_core, y_min = light_lim95$y_min/area_core, y_max = light_lim95$y_max/area_core, x_max = light_lim95$x_max, y_values = 10^c(0,2,4,6))
   assign(paste0(name_i, '_lightbin_plots'), light_plots_i)
 }
 
 # 8. Overlaid histogram plots: density scalings split in a 3x3 group. Both types of splits.
 
 xl2 <- 'Diameter (cm)'
 yl2 <- expression(paste('Density (individuals ha'^-1,')', sep=''))
 
 th1 <- theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10), legend.position = c(0.8, 0.8))
 
 threebin_names <- c('threebins_even_1990', 'threebins_quantile_1990', 'threebins_even_1995', 'threebins_quantile_1995')
 
 for (i in threebin_names) {
   overlay_lim <- get_axis_limits(i)
   overlay_low <- plotlogbin_overlaid(get(i)[[1]], xl2, yl2, 'Low light', NULL, reg = FALSE, cutoff = NA, y_min = overlay_lim$y_min, y_max = overlay_lim$y_max, x_max = overlay_lim$x_max, y_values = -4:3, plottype = 'point') + th1
   overlay_mid <- plotlogbin_overlaid(get(i)[[2]], xl2, yl2, 'Medium light', NULL, reg = FALSE, cutoff = NA, y_min = overlay_lim$y_min, y_max = overlay_lim$y_max, x_max = overlay_lim$x_max, y_values = -4:3, plottype = 'point') + th1
   overlay_high <- plotlogbin_overlaid(get(i)[[3]], xl2, yl2, 'High light', NULL, reg = FALSE, cutoff = NA, y_min = overlay_lim$y_min, y_max = overlay_lim$y_max, x_max = overlay_lim$x_max, y_values = -4:3, plottype = 'point') + th1
   
   pg <- plot_grid(overlay_low, overlay_mid, overlay_high, nrow = 3)
   ggsave(file.path(fpfig, paste0(i,'.png')), pg, height = 9, width = 4, dpi = 400)
 }
 