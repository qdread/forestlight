# Functions for analyzing and plotting BCI data.
# 27 July 2017.


logbin <- function(x, y = NULL, n) {
  logx <- log10(x)                                           # log transform x value (biomass)
  bin_edges <- seq(min(logx), max(logx), length.out = n + 1) # get edges of bins
  logxbin <- rep(NA, length(logx))                           # create data structure to assign trees to bins
  b <- bin_edges                                             # add a little to the biggest bin temporarily
  b[length(b)] <- b[length(b)] + 1                           # (so that the biggest single tree is put in a bin)
  for (i in 1:length(logx)) {
    logxbin[i] <- sum(logx[i] >= b)                          # assign each tree to a bin
  }
  bin_midpoints <- numeric(n)
  for (i in 1:n) {
    bin_midpoints[i] <- mean(10^(bin_edges[i:(i+1)]))        # backtransform bin edges to linear, and get midpoints
  }
  bin_widths <- diff(10^bin_edges)                           # get linear width of each bin
  bin_factor <- factor(logxbin, levels=1:n)                  # convert bin to factor (required to deal with zeroes if present)
  bin_counts <- table(bin_factor)                            # find number of trees in each bin
  if (!is.null(y)) {
    rawy <- tapply(y, bin_factor, sum)                       # sum y value (production) in each bin
    rawy[is.na(rawy)] <- 0                                   # add zeroes back in if present
    bin_values <- as.numeric(rawy/bin_widths)                # divide production by width for each bin 
  }
  else {
    bin_values <- as.numeric(bin_counts/bin_widths)          # 1-dimensional case.
  }
  
  return(data.frame(bin_midpoint = bin_midpoints,            # return result!
                    bin_value = bin_values,                  # also add bin min and max for bar plot purposes
                    bin_count = as.numeric(bin_counts),
                    bin_min = 10^bin_edges[1:n],
                    bin_max = 10^bin_edges[2:(n+1)]))
  
}

# Function for plotting.

plotlogbin_cutoff <- function(dat, xl, yl, plottitle, plotsubtitle=NULL, reg = FALSE, cutoff = NA, y_min=0.1, y_max=53, x_min=1.1, x_max=141, plotarea=50, y_values=-1:3, x_values = 0:2) {
  
  dat <- transform(dat, bin_value = bin_value/plotarea) # kg per hectare.
  
  p <- ggplot(dat, aes(xmin=bin_min, xmax=bin_max, ymin=0, ymax=bin_value)) + 
    geom_rect(alpha = 0.5) +
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

#gammainc <- function(a, x) gamma(a) * pgamma(x, a, 1, lower = FALSE)

nll_powerlaw <- function(alpha, xmin) {
  C <- (alpha - 1) * ( xmin ^ (alpha - 1) )
  fx <- x ^ ( -alpha )
  px <- C * fx
  -sum(log(px))
} 

# See: http://www.stat.cmu.edu/~cshalizi/2010-10-18-Meetup.pdf
nll_powerlaw_cutoff2 <- function(alpha, xmin, L) {
  C <- (1/L) / (gsl::gamma_inc(1-alpha, xmin/L))
  fx <- ( (x/L)^ -alpha ) * exp(-x/L)
  px <- C * fx
  -sum(log(px))
}

# Also include the pareto with exponential cutoff fits.
pareto_cutoff <- function(x, alpha, xmin, L) {
  C <- (1/L) / (gsl::gamma_inc(1-alpha, xmin/L))
  fx <- ( (x/L)^ -alpha ) * exp(-x/L)
  px <- C * fx
  return(px)
}
pareto_cutoff_n <- function(x, alpha, xmin, L, n) {
  C <- (1/L) / (gsl::gamma_inc(1-alpha, xmin/L))
  fx <- ( (x/L)^ -alpha ) * exp(-x/L)
  px <- C * fx
  return(px*n)
}

AICc <- function(n, k, lhat) {
  2 * k - 2 * -lhat + ( 2 * k * (k + 1) ) / ( n - k )
}

powerlawfit <- function(dat) {
  library(poweRlaw)
  pl_dat <- conpl$new(dat)
  lognorm_dat <- conlnorm$new(dat)
  xmin_pl <- pl_dat$getXmin()
  xmin_lognorm <- lognorm_dat$getXmin()
  pars_pl <- estimate_pars(pl_dat)
  pars_lognorm <- estimate_pars(lognorm_dat)
  pl_dat$setPars(pars_pl)
  lognorm_dat$setPars(pars_lognorm)
  plotdat <- plot(pl_dat)
  plfit_dat <- lines(pl_dat)
  lognormfit_dat <- lines(lognorm_dat)
  pl_pdf <- dist_pdf(m = pl_dat, q = plfit_dat$x, log = FALSE)
  lognorm_pdf <- dist_pdf(m = lognorm_dat, q = lognormfit_dat$x, log = FALSE)
  
  # bootstrap confidence interval of Pareto fit
  # discard 500 burnin iterations
  n_boot <- 1499
  n_burn <- 500
  pl_boot <- bootstrap(m = pl_dat, xmins = pl_dat$getXmin(), no_of_sims = n_boot, threads = 3)
  boot_ci <- quantile(pl_boot$bootstraps$pars[-(1:n_burn)], probs = c(0.025, 0.975))
  
  return(list(plotdat = plotdat, 
              plfit = plfit_dat, 
              lognormfit = lognormfit_dat, 
              plpdf = data.frame(x = plfit_dat$x, y = pl_pdf),
              lognormpdf = data.frame(x = lognormfit_dat$x, y = lognorm_pdf),
              xmin = xmin_pl, 
              alpha = pars_pl$pars,
              xmin_lognorm = xmin_lognorm,
              pars_lognorm = pars_lognorm$pars,
              boot_ci = as.numeric(boot_ci)
  ))
}

plotbinsandfits_cutoff <- function(pl, plc, bindat, plottitle = 'plot title', xl = 'x label', yl = 'log PDF', plotarea=50, y_min=0.001, y_max=1075, x_min = 1, x_max=141, y_values = c(-3,-1,1,3), x_values = 1:3) {
  
  expr1 <- as.character(as.expression(substitute(
    "Pareto:"~~alpha == a, list(a = round(pl$alpha, 2)))))
  expr2 <- as.character(as.expression(substitute(
    "Pareto-cutoff:"~~alpha == a*","~~L == lval, list(a = round(plc@coef[1], 2), 
                                                      lval = round(plc@coef[2], 2)))))
  
  
  bindat <- transform(bindat, bin_value = bin_value/plotarea) # Trees per hectare.
  bindat <- subset(bindat, bin_value > 0)
  #y_min <- 10^floor(log10(min(bindat$bin_value, na.rm = TRUE)))
  #y_max <- max(bindat$bin_value, na.rm = TRUE) * 1.1
  
  n_indivs <- sum(bindat$bin_count)
  #y_min <- 0.001
  #y_max <- 1075
  #x_max <- 141
  
  p <- ggplot(bindat) + 
    scale_x_log10(name = xl,
                  limits = c(1, x_max),
                  breaks = 10^x_values,
                  labels = 10^x_values,
                  expand = c(0,0)) +
    scale_y_log10(name = yl,
                  limits = c(y_min, y_max),
                  breaks = 10^y_values,
                  labels = as.character(10^y_values),
                  expand = c(0,0)) +
    # scale_y_log10(name = yl,
    #               limits = c(y_min, y_max),
    #               breaks = scales::trans_breaks("log10", function(x) 10^x),
    #               labels = scales::trans_format("log10", scales::math_format(10^.x)),
    #               expand=c(0,0)) +
    geom_rect(aes(xmin = bin_min, xmax = bin_max, ymin = 0, ymax = bin_value), alpha = 0.5) +
    geom_line(data = transform(subset(pl$plpdf, y>0), y=y*n_indivs/plotarea), aes(x,y), color = 'forestgreen', size = 2) +
    stat_function(fun = pareto_cutoff_n, args = list(alpha=plc@coef[1], xmin = 1, L = plc@coef[2], n = n_indivs/plotarea), color = 'goldenrod', size=2) +
    geom_text(x = -Inf, y = -Inf, label = expr1, parse = TRUE, hjust = 0, vjust = -1.5) +
    geom_text(x = -Inf, y = -Inf, label = expr2, parse = TRUE, hjust = 0, vjust = -0.2) +
    panel_border(colour = 'black') +
    labs(x = xl, y = yl) +
    ggtitle(plottitle)
  return(p)
}

boot_mle <- function(xboot, nboot, L_init = 1000) {
  boot_stats <- list()
  pb <- txtProgressBar(0, nboot, style = 3)
  for (i in 1:nboot) {
    setTxtProgressBar(pb, i)
    x <<- sample(xboot, size = length(xboot), replace = TRUE) # Sets x in global environment.
    fit1_i <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
    fit2_i <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=L_init), fixed = list(xmin=min(x)), method='BFGS')
    boot_stats[[i]] <- c(slope1=fit1_i@coef[1], slope2=fit2_i@coef[1], cutoff=fit2_i@coef[2])
  }
  close(pb)
  do.call('rbind', boot_stats)
}

twoslopeplot <- function(dat, plottitle = 'plot title', xl = 'x label', yl = 'y label', binvar='diameter') {
  if (binvar == 'mass') dat$yv <- dat$agb_corr else dat$yv <- dat$dbh_corr
  lmsize <- lm(log10(production) ~ log10(yv), data=dat)
  lmsizecomp <- lm(log10(production) ~ log10(yv) + log10(light_received), data=dat)
  dbh_slope <- lmsizecomp$coefficients[2] # extract slope from full model
  # refit model, setting slope from full model and estimating intercept
  adjusted_lm <- lm(log10(production) ~ 1 + offset(dbh_slope * log10(yv)), data=dat) 
  adjusted_intercept <- adjusted_lm$coefficients[1]
  
  ggplot(dat, aes(x = yv, y = production)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = lmsize$coefficients[2], intercept = lmsize$coefficients[1], color = 'indianred', size = 2) +
    geom_abline(slope = dbh_slope, intercept = adjusted_intercept, color = 'steelblue1', size = 2) +
    scale_x_log10(name = xl,
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    scale_y_log10(name = yl,
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggtitle(plottitle) +
    panel_border(colour = 'black')
}

AICc <- function(n, k, lhat) {
  2 * k - 2 * -lhat + ( 2 * k * (k + 1) ) / ( n - k )
}

# Yet another log binning function that just takes existing bin edges and puts the trees into it.
# Redo log bin with same width.

logbin_setedges <- function(x, y = NULL, edges) {
  logx <- log10(x)                                           # log transform x value (biomass)
  bin_edges <- log10(c(edges$bin_min, edges$bin_max[length(edges$bin_max)]))
  n <- length(edges$bin_min)
  logxbin <- rep(NA, length(logx))                           # create data structure to assign trees to bins
  b <- bin_edges                                             # add a little to the biggest bin temporarily
  b[length(b)] <- b[length(b)] + 1                           # (so that the biggest single tree is put in a bin)
  for (i in 1:length(logx)) {
    logxbin[i] <- sum(logx[i] >= b)                          # assign each tree to a bin
  }
  bin_midpoints <- edges$bin_midpoint
  bin_widths <- diff(10^bin_edges)                           # get linear width of each bin
  bin_factor <- factor(logxbin, levels=1:n)                  # convert bin to factor (required to deal with zeroes if present)
  bin_counts <- table(bin_factor)                            # find number of trees in each bin
  if (!is.null(y)) {
    rawy <- tapply(y, bin_factor, sum)                       # sum y value (production) in each bin
    rawy[is.na(rawy)] <- 0                                   # add zeroes back in if present
    bin_values <- as.numeric(rawy/bin_widths)                # divide production by width for each bin 
  }
  else {
    bin_values <- as.numeric(bin_counts/bin_widths)          # 1-dimensional case.
  }
  
  return(data.frame(bin_midpoint = bin_midpoints,            # return result!
                    bin_value = bin_values,                  # also add bin min and max for bar plot purposes
                    bin_count = as.numeric(bin_counts),
                    bin_min = 10^bin_edges[1:n],
                    bin_max = 10^bin_edges[2:(n+1)]))
  
}

plotlogbin_overlaid <- function(dat, xl, yl, plottitle, plotsubtitle=NULL, reg = FALSE, cutoff = NA, y_min=0.1, y_max=53, x_min=1.1, x_max=141, plotarea=50, y_values=-1:3, x_values = 0:2) {
  
  dat <- transform(dat, bin_value = bin_value/plotarea) # kg per hectare.
  
  p <- ggplot(dat, aes(xmin=bin_min, xmax=bin_max, ymin=0, ymax=bin_value)) + 
    geom_rect(alpha = 1, aes(group=guild, fill=guild), color = 'black') +
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
