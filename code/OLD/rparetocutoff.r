# Create the cutoff distributions from the pareto_cutoff density function 
# Sample the same number of individuals from gap, intermediate and shade distributions.
# See whether we get as many big shade trees and few big gap trees when we sample the same large number of individuals from each of the distributions.

pareto_cutoff <- function(x, alpha, xmin, L) {
  C <- (1/L) / (expint::gammainc(1-alpha, xmin/L))
  fx <- ( (x/L)^ -alpha ) * exp(-x/L)
  px <- C * fx
  return(px)
}

rparetocutoff <- function(n, alpha, xmin, L) {
  require(distr)
  mydist <- AbscontDistribution(d = function(x) pareto_cutoff(x, alpha=alpha, xmin=xmin, L=L), low1=xmin, withStand = TRUE)
  r(mydist)(n)
}

plotlogbin_cutoff <- function(dat, xl, yl, plottitle, reg = FALSE, cutoff = NA, y_min=0.1, y_max=53, x_max=141, plotarea=50, powers=0:6) {
  
  dat <- transform(dat, bin_value = bin_value/plotarea) # kg per hectare.
  
  p <- ggplot(dat, aes(xmin=bin_min, xmax=bin_max, ymin=0, ymax=bin_value)) + 
    geom_rect(alpha = 0.5) +
    scale_x_log10(name = xl, expand = c(0,0),
                  breaks = c(1,10,100),
                  labels = c(1,10,100), limits=c(1,x_max)) +
    scale_y_log10(name = yl, expand = c(0,0), limits = c(y_min, y_max),
                  breaks = 10^powers,
                  labels = as.character(10^powers)) +
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
  return(p)
}

# For each group, sample many random draws with the correct parameters.

ntrees <- 1e6

shade_sample <- rparetocutoff(n=ntrees, alpha=coefdat$slope2[2], xmin=1, L=coefdat$cutoff[2])
inter_sample <- rparetocutoff(n=ntrees, alpha=coefdat$slope2[3], xmin=1, L=coefdat$cutoff[3])
gap_sample <- rparetocutoff(n=ntrees, alpha=coefdat$slope2[4], xmin=1, L=coefdat$cutoff[4])

quantile(shade_sample, probs=.999)
quantile(inter_sample, probs=.999)
quantile(gap_sample, probs=.999)

# Plot the samples.
shade_bin <- logbin(x=shade_sample, n=20)
inter_bin <- logbin(x=inter_sample, n=20)
gap_bin <- logbin(x=gap_sample, n=20)
plotlogbin_cutoff(dat=shade_bin, xl='Diameter (cm)', yl='Density', plottitle='Random sample from shade distribution', y_min=1, y_max=1e6, x_max=200, plotarea=1)
plotlogbin_cutoff(dat=inter_bin, xl='Diameter (cm)', yl='Density', plottitle='Random sample from intermediate distribution', y_min=1, y_max=1e6, x_max=200, plotarea=1)
plotlogbin_cutoff(dat=gap_bin, xl='Diameter (cm)', yl='Density', plottitle='Random sample from gap distribution', y_min=1, y_max=1e6, x_max=200, plotarea=1)
