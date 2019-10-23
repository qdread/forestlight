# New plotting/binning functions for the revamped analysis
# 30 Aug 2017

# Get the maximum value and axis scale numbers from a list.

get_axis_limits <- function(list_of_bins) {
  # progressively flatten the list of bins.
  biglist <- lapply(list_of_bins, get)
  if (inherits(biglist[[1]], 'list')) biglist <- lapply(biglist, function(x) do.call('rbind', x))
  allvals <- do.call('rbind', biglist)
  y_min <- 10 ^ floor(log10(min(allvals$bin_value)))
  y_max <- max(allvals$bin_value) * 1.1
  x_max <- max(allvals$bin_max)
  return(list(y_min=y_min, y_max=y_max, x_max=x_max))
}

plotbinsandfits <- function(pl, bindat, plottitle = 'plot title', xl = 'x label', yl = 'log PDF', plottype = 'bin') {
  
  expr1 <- as.character(as.expression(substitute(
    "Pareto:"~~alpha == a, list(a = round(pl$alpha, 2)))))
  
  bindat <- transform(bindat, bin_value = bin_value/sum(bin_count))
  bindat <- subset(bindat, bin_value > 0)
  y_min <- 10^floor(log10(min(bindat$bin_value, na.rm = TRUE)))
  y_max <- max(bindat$bin_value, na.rm = TRUE) * 1.1
  
  if (plottype == 'bin') data_geom <- geom_rect(aes(xmin = bin_min, xmax = bin_max, ymin = 0, ymax = bin_value), alpha = 0.5)
  if (plottype == 'point') data_geom <- geom_point(aes(x = bin_midpoint, y = bin_value), size = 1.5)
  
  p <- ggplot(bindat) + 
    scale_x_log10(name = xl,
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  expand = c(0,0)) +
    scale_y_log10(name = yl,
                  limits = c(y_min, y_max),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  expand=c(0,0)) +
    data_geom +
    geom_line(data = subset(pl$plpdf, y>0), aes(x,y), color = 'forestgreen', size = 2) +
    geom_text(x = -Inf, y = -Inf, label = expr1, parse = TRUE, hjust = 0, vjust = -1.5) +
    panel_border(colour = 'black') +
    labs(x = xl, y = yl) +
    ggtitle(plottitle)
  return(p)
}


# Plot bins and fits with cutoff and set axis limits, giving the option to plot either points or bars.

plotbinsandfits_cutoff <- function(pl, plc, bindat, plottitle = 'plot title', xl = 'x label', yl = 'log PDF', plotarea=50, y_min=0.001, y_max=1075, x_max=141, plottype = 'bin', x_values=10^(0:2), y_values=10^c(-3,-1,1,3)) {
  
  expr1 <- as.character(as.expression(substitute(
    "Pareto:"~~alpha == a, list(a = round(pl$alpha, 2)))))
  expr2 <- as.character(as.expression(substitute(
    "Pareto-cutoff:"~~alpha == a*","~~L == lval, list(a = round(plc@coef[1], 2), 
                                                      lval = round(plc@coef[2], 2)))))
  
  
  bindat <- transform(bindat, bin_value = bin_value/plotarea) # Trees per hectare.
  bindat <- subset(bindat, bin_value > 0)

  n_indivs <- sum(bindat$bin_count)
  
  if (plottype == 'bin') data_geom <- geom_rect(aes(xmin = bin_min, xmax = bin_max, ymin = 0, ymax = bin_value), alpha = 0.5)
  if (plottype == 'point') data_geom <- data_geom <- geom_point(aes(x = bin_midpoint, y = bin_value), size = 1.5)

  p <- ggplot(bindat) + 
    scale_x_log10(name = xl,
                  limits = c(1, x_max),
                  breaks = x_values,
                  labels = as.character(x_values),
                  expand = c(0,0)) +
    scale_y_log10(name = yl,
                  limits = c(y_min, y_max),
                  breaks = y_values,
                  labels = as.character(y_values),
                  expand = c(0,0))
  if (plottype == 'bin') p <- p + data_geom
  p <- p +
    geom_line(data = transform(subset(pl$plpdf, y>0), y=y*n_indivs/plotarea), aes(x,y), color = 'forestgreen', size = 2) +
    stat_function(fun = pareto_cutoff_n, args = list(alpha=plc@coef[1], xmin = 1, L = plc@coef[2], n = n_indivs/plotarea), color = 'goldenrod', size=2)
  if (plottype == 'point') p <- p + data_geom
  p <- p +
    geom_text(x = -Inf, y = -Inf, label = expr1, parse = TRUE, hjust = 0, vjust = -1.5) +
    geom_text(x = -Inf, y = -Inf, label = expr2, parse = TRUE, hjust = 0, vjust = -0.2) +
    panel_border(colour = 'black') +
    labs(x = xl, y = yl) +
    ggtitle(plottitle)
  return(p)
}

# Two-slope plot with different geom

twoslopeplot <- function(dat, plottitle = 'plot title', xl = 'x label', yl = 'y label', binvar='diameter', plottype = 'point') {
  if (binvar == 'mass') dat$yv <- dat$agb_corr else dat$yv <- dat$dbh_corr
  lmsize <- lm(log10(production) ~ log10(yv), data=dat)
  lmsizecomp <- lm(log10(production) ~ log10(yv) + log10(light_received), data=dat)
  dbh_slope <- lmsizecomp$coefficients[2] # extract slope from full model
  # refit model, setting slope from full model and estimating intercept
  adjusted_lm <- lm(log10(production) ~ 1 + offset(dbh_slope * log10(yv)), data=dat) 
  adjusted_intercept <- adjusted_lm$coefficients[1]
  
  p <- ggplot(dat, aes(x = yv, y = production))
  
  if (plottype == 'point') p <- p + geom_point(alpha = 0.6)
  if (plottype == 'hex') p <- p + geom_hex() + 
    scale_fill_gradientn(colours = colorRampPalette(RColorBrewer::brewer.pal(9, 'YlOrRd'), bias=2)(10)) +
    theme(legend.position = c(0.8, 0.2))

  
  p <- p + geom_abline(slope = lmsize$coefficients[2], intercept = lmsize$coefficients[1], color = 'indianred', size = 2) +
    geom_abline(slope = dbh_slope, intercept = adjusted_intercept, color = 'steelblue1', size = 2) +
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


# Plot log bin and cutoff, used for binned production plots.
plotlogbin_cutoff <- function(dat, xl, yl, plottitle, reg = FALSE, cutoff = NA, y_min=0.1, y_max=53, x_max=141, plotarea=50, plottype = 'bin', y_values=10^c(-1,0,1)) {
  
  dat <- transform(dat, bin_value = bin_value/plotarea) # kg per hectare.
  dat <- subset(dat, bin_value > 0)
  
  p <- ggplot(dat, aes(xmin=bin_min, xmax=bin_max, ymin=0, ymax=bin_value))
  if (plottype == 'bin') p <- p + geom_rect(alpha = 0.5)
  p <- p +
    scale_x_log10(name = xl, expand = c(0,0),
                  breaks = c(1,10,100),
                  labels = c(1,10,100), limits=c(1,x_max)) +
    scale_y_log10(name = yl, expand = c(0,0), limits = c(y_min, y_max),
                  breaks = y_values,
                  labels = as.character(y_values)) +
    panel_border(colour = 'black') + 
    ggtitle(plottitle)
  if (reg) {
    p <- p +
      stat_smooth(method = 'lm', se = FALSE, color = 'forestgreen', size = 2,
                  aes(x = bin_midpoint, y = bin_value)) +
      geom_text(x = -Inf, y = -Inf, 
                label = paste('Slope without cutoff:', 
                              round(lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=dat)$coef[2], 2)),
                hjust = 0, vjust = -1.5)
    if (!is.na(cutoff)) {
      p <- p +
        stat_smooth(method = 'lm', se = FALSE, color = 'goldenrod', size = 2,
                    aes(x = bin_midpoint, y = bin_value), data = subset(dat, bin_midpoint <= cutoff)) +
        geom_text(x = -Inf, y = -Inf, 
                  label = paste('Slope with cutoff:', 
                                round(lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=subset(dat, bin_midpoint <= cutoff))$coef[2], 2)),
                  hjust = 0, vjust = -0.2)
    }
  }
  if (plottype == 'point') p <- p + geom_point(aes(x = bin_midpoint, y = bin_value), size = 1.5)
  return(p)
}

# Overlaid histogram plot

plotlogbin_overlaid <- function(dat, xl, yl, plottitle, plotsubtitle=NULL, reg = FALSE, cutoff = NA, y_min=0.1, y_max=53, x_min=1.1, x_max=141, plotarea=50, y_values=-1:3, x_values = 0:2, plottype = 'bin') {
  
  dat <- transform(dat, bin_value = bin_value/plotarea) # kg per hectare.
  
  p <- ggplot(dat, aes(xmin=bin_min, xmax=bin_max, ymin=0, ymax=bin_value))
  if (plottype == 'bin') p <- p + geom_rect(alpha = 1, aes(group=guild, fill=guild), color = 'black')
  if (plottype == 'point') p <- p + geom_point(aes(group = guild, colour = guild, x=bin_midpoint, y=bin_value))
  p <- p +
    scale_x_log10(name = xl, expand = c(0,0),
                  breaks = 10^x_values,
                  labels = as.character(10^x_values), limits=c(x_min,x_max)) +
    scale_y_log10(name = yl, expand = c(0,0), limits = c(y_min, y_max),
                  breaks = 10^y_values,
                  labels = as.character(10^y_values)) +
    panel_border(colour = 'black') + 
    ggtitle(plottitle, plotsubtitle)
  if (reg) {
    p <- p +
      stat_smooth(method = 'lm', se = FALSE, color = 'forestgreen', size = 2,
                  aes(x = bin_midpoint, y = bin_value)) +
      geom_text(x = -Inf, y = -Inf, 
                label = paste('Slope without cutoff:', 
                              round(lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=dat)$coef[2], 2)),
                hjust = 0, vjust = -1.5)
    
    if (!is.na(cutoff)) {
      p <- p +
        stat_smooth(method = 'lm', se = FALSE, color = 'goldenrod', size = 2,
                    aes(x = bin_midpoint, y = bin_value), data = subset(dat, bin_midpoint <= cutoff)) +
        geom_text(x = -Inf, y = -Inf, 
                  label = paste('Slope with cutoff:', 
                                round(lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=subset(dat, bin_midpoint <= cutoff))$coef[2], 2)),
                  hjust = 0, vjust = -0.2)
    }
  }
  return(p)
}
