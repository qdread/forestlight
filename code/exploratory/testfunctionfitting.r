# Test Caroline's function fitting code

load(file = 'C:/Users/Q/google_drive/ForestLight/data/data_06sep/bci_data_object.RData')

source('~/GitHub/FunctionFitting/PowerLawFit_etc_160509.r')
source('code/allfunctions27july.r')

# Test powerlaw fit with the same set xmin as before. Use all tree data.
# Convert to mm, and round to the nearest mm

x_test <- alltreedat[[6]]$dbh_corr
x_test <- round(x_test*10, 0)

x_powerlaw <- powerlaw_fit(x = x_test, xmin = min(x_test), interval = c(1.2, 3), plot = TRUE, lines = TRUE)
x_powerexp <- powerexp_fit(x = x_test, xmin1 = min(x_test), alpha1 = x_powerlaw$alpha, plotting = TRUE)

x_shade <- subset(alltreedat[[6]], tol_wright == 'S')$dbh_corr
x_shade <- round(x_shade*10, 0)

x_shadepowerlaw <- powerlaw_fit(x = x_shade, xmin = min(x_shade), interval = c(1.2, 3), plot = TRUE, lines = TRUE)
x_shadepowerexp <- powerexp_fit(x = x_shade, xmin1 = min(x_shade), alpha1 = x_shadepowerlaw$alpha, plotting = TRUE)

x_gap <- subset(alltreedat[[6]], tol_wright == 'G')$dbh_corr
x_gap <- round(x_gap*10, 0)

x_gappowerlaw <- powerlaw_fit(x = x_gap, xmin = min(x_gap), interval = c(1.2, 3), plot = TRUE, lines = TRUE)
x_gappowerexp <- powerexp_fit(x = x_gap, xmin1 = min(x_gap), alpha1 = x_gappowerlaw$alpha, plotting = TRUE)


pareto_fn <- function(x, alpha, C) {
  C * x ^ -alpha
}

exponential_fn <- function(x, beta, C1, Cx) {
  Cx * C1 * exp(-beta * x)
}

plotbinsandfits_piecewise <- function(pl, plc, bindat, plottitle = 'plot title', xl = 'x label', yl = 'log PDF', plotarea=50, y_min=0.001, y_max=1075, x_max=141, plottype = 'bin', x_values=10^(0:2), y_values=10^c(-3,-1,1,3), x_units = 'cm') {
  
  
  expr1 <- as.character(as.expression(substitute(
    alpha ==a*","~~min==minval*","~~max==maxval, list(a = round(pl$alpha, 2), minval = round(pl$xmin/10, 2), maxval = round(plc$xmin2/10, 2)))))

  
  
  bindat <- transform(bindat, bin_value = bin_value/plotarea) # Trees per hectare.
  bindat <- subset(bindat, bin_value > 0)
  
  n_indivs <- sum(bindat$bin_count)
  
  if (plottype == 'bin') data_geom <- geom_rect(aes(xmin = bin_min, xmax = bin_max, ymin = 0, ymax = bin_value), alpha = 0.5)
  if (plottype == 'point') data_geom <- data_geom <- geom_point(aes(x = bin_midpoint, y = bin_value), size = 1.5)
  
  x_labels <- x_values
  x_min <- 1
  if (x_units == 'mm') {
    x_labels <- x_labels/10
    x_min <- 10
  }
  
  p <- ggplot(bindat) + 
    scale_x_log10(name = xl,
                  limits = c(x_min, x_max),
                  breaks = x_values,
                  labels = as.character(x_labels),
                  expand = c(0,0)) +
    scale_y_log10(name = yl,
                  limits = c(y_min, y_max),
                  breaks = y_values,
                  labels = as.character(y_values),
                  expand = c(0,0))
  if (plottype == 'bin') p <- p + data_geom
  p <- p +
    stat_function(fun = pareto_fn, xlim = log10(c(pl$xmin, x_max)),
                  args = list(alpha = pl$alpha, C = pl$C * n_indivs / plotarea), color = 'forestgreen', size = 2)

    # Only plot the exponential function if the cutoff is within the range of the data
  if (plc$xmin2 < max(bindat$bin_max)) {
  p <- p +
    stat_function(fun = exponential_fn, xlim = log10(c(plc$xmin2, x_max)),
                  args = list(beta = plc$beta, C1 = plc$C1, Cx = plc$Cx * n_indivs / plotarea), color = 'goldenrod', size = 2)
  }
  if (plottype == 'point') p <- p + data_geom
  p <- p +
    geom_text(x = -Inf, y = -Inf, label = expr1, parse = TRUE, hjust = 0, vjust = -1) +
    panel_border(colour = 'black') +
    labs(x = xl, y = yl) +
    ggtitle(plottitle)
  
  return(p)
}

# Plot for shade trees
shade_logbin <- logbin(x = x_shade, y = NULL, n = 20)


plotbinsandfits_piecewise(x_shadepowerlaw, x_shadepowerexp, shade_logbin, plottitle = 'Shade tolerant trees 2010', xl = 'Diameter', yl = 'Density', plottype = 'point', plotarea = 42.84, y_min = 0.00001, y_max = 10000, x_max = 2000, x_values=10^(1:3), y_values = 10^c(-3,-1,1,3), x_units='mm')

ggsave('C:/Users/Q/google_drive/ForestLight/figs/new_cutoff_plots/shade2010.pdf')

gap_logbin <- logbin(x = x_gap, y = NULL, n = 20)


plotbinsandfits_piecewise(x_gappowerlaw, x_gappowerexp, gap_logbin, plottitle = 'Gap trees 2010', xl = 'Diameter', yl = 'Density', plottype = 'point', plotarea = 42.84, y_min = 0.00001, y_max = 10000, x_max = 2000, x_values=10^(1:3), y_values = 10^c(-3,-1,1,3), x_units='mm')

ggsave('C:/Users/Q/google_drive/ForestLight/figs/new_cutoff_plots/gap2010.pdf')



# Fitting both xmin and cutoff.
# xmin should indicate where the power law kicks in.
# Below that, it's light limitation?


x_shade <- subset(alltreedat[[6]], tol_wright == 'S')$dbh_corr
x_shade <- round(x_shade*10, 0)

x_shadepowerlaw <- powerlaw_fit(x = x_shade, interval = c(1.2, 3), plot = TRUE, lines = TRUE)
x_shadepowerexp <- powerexp_fit(x = x_shade, xmin1 = x_shadepowerlaw$xmin, alpha1 = x_shadepowerlaw$alpha, plotting = TRUE)

x_gap <- subset(alltreedat[[6]], tol_wright == 'G')$dbh_corr
x_gap <- round(x_gap*10, 0)

x_gappowerlaw <- powerlaw_fit(x = x_gap, interval = c(1.2, 3), plot = TRUE, lines = TRUE)
x_gappowerexp <- powerexp_fit(x = x_gap, xmin1 = x_gappowerlaw$xmin, alpha1 = x_gappowerlaw$alpha, plotting = TRUE)
