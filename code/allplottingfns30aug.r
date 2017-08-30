# New plotting/binning functions for the revamped analysis
# 30 Aug 2017

# Get the maximum value and axis scale numbers from a list.

get_axis_limits <- function(list_of_bins) {
  # progressively flatten the list of bins.
  biglist <- lapply(list_of_bins, function(x) do.call('rbind', get(x)))
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

plotbinsandfits_cutoff <- function(pl, plc, bindat, plottitle = 'plot title', xl = 'x label', yl = 'log PDF', plotarea=50, y_min=0.001, y_max=1075, x_max=141, plottype = 'bin', y_values=10^c(-3,-1,1,3)) {
  
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
                  breaks = c(1,10,100),
                  labels = c(1,10,100),
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

