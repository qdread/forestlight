# Alternative density plots with different axis labels.

source('~/GitHub/forestlight/code/load50ha.r')
source('~/GitHub/forestlight/code/powerlawnoboot.r')

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

plotbinsandfits_cutoff <- function(pl, plc, bindat, plottitle = 'plot title', xl = 'x label', yl = 'log PDF', plotarea=50, y_min=0.001, y_max=1075, x_max=141) {
  
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
                  breaks = c(1,10,100),
                  labels = c(1,10,100),
                  expand = c(0,0)) +
    scale_y_log10(name = yl,
                  limits = c(y_min, y_max),
                  breaks = 10^c(-3,-1,1,3),
                  labels = as.character(10^c(-3,-1,1,3)),
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


# Fit power laws
numbins <- 20 # Can be edited if desired. Only for visualization purposes.

alltreedat <- subset(bcicensusdat, !is.na(dbh) & production67 > 0 & comp_idx > 0)
shadedat <- subset(bcicensusdat, !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'S')
intdat <- subset(bcicensusdat, !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'I')
gapdat <- subset(bcicensusdat, !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'G')


# Pareto fits
bci_dens_fit_all <- powerlawfit(alltreedat$dbh, doboot = F)
bci_dens_fit_shade <- powerlawfit(shadedat$dbh, doboot = F)
bci_dens_fit_inter <- powerlawfit(intdat$dbh, doboot = F)
bci_dens_fit_gap <- powerlawfit(gapdat$dbh, doboot = F)

# Log bins
bci_dens_logbin_all <- with(alltreedat, logbin(x=dbh, y=NULL, n = numbins))
bci_dens_logbin_shade <- with(shadedat, logbin(x=dbh, y=NULL, n = numbins))
bci_dens_logbin_inter <- with(intdat, logbin(x=dbh, y=NULL, n = numbins))
bci_dens_logbin_gap <- with(gapdat, logbin(x=dbh, y=NULL, n = numbins))

# Also include the pareto with exponential cutoff fits.
pareto_cutoff <- function(x, alpha, xmin, L) {
  C <- (1/L) / (expint::gammainc(1-alpha, xmin/L))
  fx <- ( (x/L)^ -alpha ) * exp(-x/L)
  px <- C * fx
  return(px)
}
pareto_cutoff_n <- function(x, alpha, xmin, L, n) {
  C <- (1/L) / (expint::gammainc(1-alpha, xmin/L))
  fx <- ( (x/L)^ -alpha ) * exp(-x/L)
  px <- C * fx
  return(px*n)
}

xl2 <- 'Diameter (cm)'
yl2 <- expression(paste('Density (individuals ha'^-1,')', sep=''))

# For the 3 groups, concatenate the bin values in individuals per hectare to figure out y axis labels
raw_numbers <- c(bci_dens_logbin_shade$bin_value, bci_dens_logbin_inter$bin_value, bci_dens_logbin_gap$bin_value)/50
y_min <- 10^floor(log10(min(raw_numbers)))
y_max <- max(raw_numbers) * 1.1
y_values <- c(0.001, 0.1, 1, 10, 100, 1000)

x_max <- max(bci_dens_logbin_shade$bin_max) 

pdall <- plotbinsandfits_cutoff(bci_dens_fit_all, fit2, bci_dens_logbin_all,
                       plottitle = 'All species', xl = xl2, yl = yl2, y_max = 2000, x_max = 250)
pdshade <- plotbinsandfits_cutoff(bci_dens_fit_shade, fit2shade, bci_dens_logbin_shade,
                       plottitle = 'Shade-tolerant species', xl = xl2, yl = yl2)
pdint <- plotbinsandfits_cutoff(bci_dens_fit_inter, fit2int, bci_dens_logbin_inter,
                       plottitle = 'Intermediate species', xl = xl2, yl = yl2)
pdgap <- plotbinsandfits_cutoff(bci_dens_fit_gap, fit2gap, bci_dens_logbin_gap,
                       plottitle = 'Gap species', xl = xl2, yl = yl2)


##########################################
# Total production plots with cutoff line, and standardized axis labeling

plotlogbin_cutoff <- function(dat, xl, yl, plottitle, reg = FALSE, cutoff = NA, y_min=0.1, y_max=53, x_max=141, plotarea=50) {
  
  dat <- transform(dat, bin_value = bin_value/plotarea) # kg per hectare.
  
  p <- ggplot(dat, aes(xmin=bin_min, xmax=bin_max, ymin=0, ymax=bin_value)) + 
    geom_rect(alpha = 0.5) +
    scale_x_log10(name = xl, expand = c(0,0),
                  breaks = c(1,10,100),
                  labels = c(1,10,100), limits=c(1,x_max)) +
    scale_y_log10(name = yl, expand = c(0,0), limits = c(y_min, y_max),
                  breaks = 10^c(-1,0,1),
                  labels = as.character(10^c(-1,0,1))) +
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

numbins <- 20 # Can be edited if desired. ***NOT JUST FOR LOOKS***

# Log bins
bci_prod_logbin_all <- with(alltreedat, 
                            logbin(x=dbh, y=production67, n=numbins))
bci_prod_logbin_shade <- with(shadedat, 
                              logbin(x=dbh, y=production67, n=numbins))
bci_prod_logbin_inter <- with(intdat, 
                              logbin(x=dbh, y=production67, n=numbins))
bci_prod_logbin_gap <- with(gapdat, 
                            logbin(x=dbh, y=production67, n=numbins))

xl3 <- 'Diameter (cm)'
yl3 <- expression(paste('Total production (kg y'^-1,' ha'^-1, ')', sep=''))

# determine the axis limits
# For the 3 groups, concatenate the bin values in individuals per hectare to figure out y axis labels
raw_numbers <- c(bci_prod_logbin_shade$bin_value, bci_prod_logbin_inter$bin_value, bci_prod_logbin_gap$bin_value)/50
y_min <- 10^floor(log10(min(raw_numbers)))
y_max <- max(raw_numbers) * 1.1
y_values <- c(0.1, 1, 10)

y_min <- 10^floor(log10(min(dat$bin_value, na.rm = TRUE)))
y_max <- max(dat$bin_value, na.rm = TRUE) * 1.1

ppall <- plotlogbin_cutoff(bci_prod_logbin_all, xl3, yl3, 
           'All species', reg = TRUE, cutoff = fit2@coef[2], y_max = 120, x_max = 250)
ppshade <- plotlogbin_cutoff(bci_prod_logbin_shade, xl3, yl3,  
           'Shade-tolerant species', reg = TRUE, cutoff = fit2shade@coef[2])
ppint <- plotlogbin_cutoff(bci_prod_logbin_inter, xl3, yl3, 
           'Intermediate species', reg = TRUE, cutoff = fit2int@coef[2])
ppgap <- plotlogbin_cutoff(bci_prod_logbin_gap, xl3, yl3, 
           'Gap species', reg = TRUE, cutoff = fit2gap@coef[2])

# Slope confidence intervals
slope_cutoff <- function(dat, cutoff) {
  lm1 <- lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=dat)
  lm2 <- lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=subset(dat, bin_midpoint <= cutoff))
  ci1 <- confint(lm1)
  ci2 <- confint(lm2)
  data.frame(model=c('no cutoff','below cutoff'), est=c(lm1$coef[2], lm2$coef[2]), cimin=c(ci1[2,1], ci2[2,1]), cimax=c(ci1[2,2],ci2[2,2]))
}

slopeall <- slope_cutoff(bci_prod_logbin_all, fit2@coef[2])
slopeshade <- slope_cutoff(bci_prod_logbin_shade, fit2shade@coef[2])
slopeint <- slope_cutoff(bci_prod_logbin_inter, fit2int@coef[2])
slopegap <- slope_cutoff(bci_prod_logbin_gap, fit2gap@coef[2])
slopedat <- cbind(guild = rep(c('all','shade','intermediate','gap'),each=2), rbind(slopeall,slopeshade,slopeint,slopegap))

pslopes <- ggplot(slopedat, aes(x=guild, y=est, ymin=cimin, ymax=cimax)) +
  geom_hline(yintercept = 0, linetype='dotted', color = 'slateblue') +
  geom_pointrange() +
  facet_wrap(~ model) +
  panel_border(colour='black') +
  theme(strip.background = element_blank()) +
  labs(y = 'slope estimate')

library(gridExtra)

pdf('C:/Users/Q/google_drive/ForestLight/docs/density_production_plots_cutoffs.pdf', height=9, width=9)
  grid.arrange(pdall, pdshade, pdint, pdgap, nrow=2)
  grid.arrange(ppall, ppshade, ppint, ppgap, nrow=2)
  pslopes
dev.off()

ggsave('C:/Users/Q/google_drive/ForestLight/figs/density_scaling_shade.png', pdshade, height=4.5, width=4.5, dpi=400)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/density_scaling_inter.png', pdint, height=4.5, width=4.5, dpi=400)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/density_scaling_gap.png', pdgap, height=4.5, width=4.5, dpi=400)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/prod_scaling_shade.png', ppshade, height=4.5, width=4.5, dpi=400)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/prod_scaling_inter.png', ppint, height=4.5, width=4.5, dpi=400)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/prod_scaling_gap.png', ppgap, height=4.5, width=4.5, dpi=400)

#################################
# added 09 May: edit axis labels for the scatter plots as well.

twoslopeplot <- function(dat, plottitle = 'plot title', xl = 'x label', yl = 'y label', binvar='diameter', x_max=141, y_min=1e-4, y_max=1000) {
  if (binvar == 'mass') dat$xv <- dat$agb else dat$xv <- dat$dbh
  lmsize <- lm(log10(production67) ~ log10(xv), data=dat)
  lmsizecomp <- lm(log10(production67) ~ log10(xv) + log10(comp_idx), data=dat)
  dbh_slope <- lmsizecomp$coefficients[2] # extract slope from full model
  # refit model, setting slope from full model and estimating intercept
  adjusted_lm <- lm(log10(production67) ~ 1 + offset(dbh_slope * log10(xv)), data=dat) 
  adjusted_intercept <- adjusted_lm$coefficients[1]
  
  ggplot(dat, aes(x = xv, y = production67)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = lmsize$coefficients[2], intercept = lmsize$coefficients[1], color = 'indianred', size = 2) +
    geom_abline(slope = dbh_slope, intercept = adjusted_intercept, color = 'steelblue1', size = 2) +
    scale_x_log10(name = xl,
                  limits = c(1, x_max),
                  breaks = c(1,10,100),
                  labels = c(1,10,100),
                  expand = c(0,0)) +
    scale_y_log10(name = yl,
                  limits = c(y_min, y_max),
                  breaks = 10^c(-3,-1,1,3),
                  labels = as.character(10^c(-3,-1,1,3)),
                  expand = c(0,0)) +
    ggtitle(plottitle) +
    panel_border(colour = 'black')
}

xl1 <- 'Diameter (cm)'
yl1 <- expression(paste('Individual production (kg y'^-1,' ha'^-1, ')', sep=''))

x_max <- max(shadedat$dbh)
y_min <- 10^floor(log10(min(shadedat$production67)))
y_max <- 1000
y_values <- c(0.001, 0.1, 10, 1000)

twoslopeplot(dat = alltreedat, 
             plottitle = 'All species', 
             xl = xl1, 
             yl =  yl1)

pipshade <- twoslopeplot(dat = shadedat, 
             plottitle = 'Shade-tolerant species', 
             xl = xl1, 
             yl =  yl1)

pipint <- twoslopeplot(dat = intdat, 
             plottitle = 'Intermediate species', 
             xl = xl1, 
             yl =  yl1)

pipgap <- twoslopeplot(dat = gapdat, 
             plottitle = 'Gap species', 
             xl = xl1, 
             yl =  yl1)

ggsave('C:/Users/Q/google_drive/ForestLight/figs/twoslope_shade.png', pipshade, height=4.5, width=4.5, dpi=400)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/twoslope_inter.png', pipint, height=4.5, width=4.5, dpi=400)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/twoslope_gap.png', pipgap, height=4.5, width=4.5, dpi=400)
