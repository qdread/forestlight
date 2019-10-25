# More plots using Nadja's modeled relative irradiance values
# QDR 26 Jun 2017

#######################################################################################
# 1. Single plot of total intercepted PAR (y) vs log diameter bins (x).
# Calculate total intercepted PAR by using allometry to get area of each tree's crown and multiply by the % par of each tree.

library(cowplot)
# Use bcicensusdat df generated by the code in the Rmd.

# Function to get tree height and crown dimensions from dbh
# Use same parameters for all trees, taken from Bohlman and O'Brien

# Function to get a rough approximation of insolation by latitude.

insolation <- function(lat) {
  lat <- lat * pi/180 # to radians
  y <- sin(lat)
  0.25 * 1367 * (1 - 0.482 * (3*y^2 - 1)/2)
}

# Insolation at BCI, 9.2 degrees N
(insol_bci <- insolation(9.2))

tp <- function(dbh) {
  h <- exp(.438 + .595 * log(dbh))    # Height
  cd <- exp(-.157 + .702 * log(dbh))  # Crown depth
  cr <- exp(-.438 + .658 * log(dbh))  # Crown radius
  cV <- exp(-.681 + 2.02 * log(dbh))  # Crown volume
  data.frame(h=h, cd=cd, cr=cr, cV=cV)
}

crowndim <- tp(bcicensusdat$dbh) 
bcicensusdat$crownarea <- pi * crowndim$cr^2
bcicensusdat <- transform(bcicensusdat, light_received = light * crownarea * insol_bci)

# Classification of light into 3 groups.
light_groups <- cut(bcicensusdat$light, breaks = 3)
light_groupcodes <- factor(light_groups, labels = c('Low light','Intermediate light','High light'))

bcicensusdat$light_group <- light_groupcodes

alltreedat <- subset(bcicensusdat, !is.na(dbh) & production34 > 0 & !is.na(light))
shadedat <- subset(bcicensusdat, !is.na(dbh) & production34 > 0 & !is.na(light) & tol_wright == 'S')
intdat <- subset(bcicensusdat, !is.na(dbh) & production34 > 0 & !is.na(light) & tol_wright == 'I')
gapdat <- subset(bcicensusdat, !is.na(dbh) & production34 > 0 & !is.na(light) & tol_wright == 'G')

################
# Run binning algorithm.

numbins <- 20 # Can be edited if desired. ***NOT JUST FOR LOOKS***

# Log bins
bci_par_logbin_all <- with(alltreedat, 
                            logbin(x=dbh, y=light_received, n=numbins))
bci_par_logbin_shade <- with(shadedat, 
                              logbin(x=dbh, y=light_received, n=numbins))
bci_par_logbin_inter <- with(intdat, 
                              logbin(x=dbh, y=light_received, n=numbins))
bci_par_logbin_gap <- with(gapdat, 
                            logbin(x=dbh, y=light_received, n=numbins))

# Energy-equivalence slopes
bci_par_lm_all <- lm(log10(bin_value) ~ log10(bin_midpoint), data = bci_par_logbin_all)
bci_par_lm_shade <- lm(log10(bin_value) ~ log10(bin_midpoint), data = bci_par_logbin_shade)
bci_par_lm_inter <- lm(log10(bin_value) ~ log10(bin_midpoint), data = bci_par_logbin_inter)
bci_par_lm_gap <- lm(log10(bin_value) ~ log10(bin_midpoint), data = bci_par_logbin_gap)

####
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

# Determine axis limits
raw_numbers <- c(bci_par_logbin_shade$bin_value, bci_par_logbin_inter$bin_value, bci_par_logbin_gap$bin_value)/50
y_min <- 10^floor(log10(min(raw_numbers)))
y_max <- max(raw_numbers) * 1.1
y_values <- 1:5


max(bci_par_logbin_all$bin_value/50)
#y_min <- 10^floor(log10(min(dat$bin_value, na.rm = TRUE)))
#y_max <- max(dat$bin_value, na.rm = TRUE) * 1.1

xl4 <- 'Diameter (cm)'
yl4 <-'Relative energy received per hectare'

ppall <- plotlogbin_cutoff(bci_par_logbin_all, xl4, yl4, 
                           'All species', reg = FALSE, cutoff = NA, y_min = 10, y_max = 1e5, x_max = 250, y_values = y_values)
ppshade <- plotlogbin_cutoff(bci_par_logbin_shade, xl4, yl4,  
                             'Shade-tolerant species', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values)
ppint <- plotlogbin_cutoff(bci_par_logbin_inter, xl4, yl4, 
                           'Intermediate species', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values)
ppgap <- plotlogbin_cutoff(bci_par_logbin_gap, xl4, yl4, 
                           'Gap species', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values)


#######################################################################################
# 2. 3x3 plot, class light availability into 3 groups (horiz) and shade tolerance groups (vert)
# In each panel, do a density scaling.
# Or do as a 3 panel plot in which there are stacked bars. Each panel is the light group and the stacks are color-coded by the shade tolerance class.

numbins <- 10 # only for visualizing

shade_lowlight_bin <-  with(subset(shadedat, light_group == 'Low light'), logbin(x=dbh, y=NULL, n = numbins))
shade_intlight_bin <-  with(subset(shadedat, light_group == 'Intermediate light'), logbin(x=dbh, y=NULL, n = numbins))
shade_highlight_bin <-  with(subset(shadedat, light_group == 'High light'), logbin(x=dbh, y=NULL, n = numbins))

int_lowlight_bin <-  with(subset(intdat, light_group == 'Low light'), logbin(x=dbh, y=NULL, n = numbins))
int_intlight_bin <-  with(subset(intdat, light_group == 'Intermediate light'), logbin(x=dbh, y=NULL, n = numbins))
int_highlight_bin <-  with(subset(intdat, light_group == 'High light'), logbin(x=dbh, y=NULL, n = numbins))

gap_lowlight_bin <-  with(subset(gapdat, light_group == 'Low light'), logbin(x=dbh, y=NULL, n = numbins))
gap_intlight_bin <-  with(subset(gapdat, light_group == 'Intermediate light'), logbin(x=dbh, y=NULL, n = numbins))
gap_highlight_bin <-  with(subset(gapdat, light_group == 'High light'), logbin(x=dbh, y=NULL, n = numbins))

# Determine axis limits
raw_numbers <- c(shade_lowlight_bin$bin_value, shade_intlight_bin$bin_value, shade_highlight_bin$bin_value,
                 int_lowlight_bin$bin_value, int_intlight_bin$bin_value, int_highlight_bin$bin_value,
                 gap_lowlight_bin$bin_value, gap_intlight_bin$bin_value, gap_highlight_bin$bin_value)/50
y_min <- 10^floor(log10(min(raw_numbers)))
y_max <- max(raw_numbers) * 1.1
y_values <- -4:3

xl2 <- 'Diameter (cm)'
yl2 <- expression(paste('Density (individuals ha'^-1,')', sep=''))

plowshade <- plotlogbin_cutoff(shade_lowlight_bin, xl2, yl2,  
                             'Shade-tolerant species', 'Low light', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) 
pintshade <- plotlogbin_cutoff(shade_intlight_bin, xl2, yl2,  
                               'Shade-tolerant species', 'Intermediate light', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values)
phighshade <- plotlogbin_cutoff(shade_highlight_bin, xl2, yl2,  
                               'Shade-tolerant species', 'High light', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) 

plowint <- plotlogbin_cutoff(int_lowlight_bin, xl2, yl2,  
                               'Intermediate species', 'Low light', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) 
pintint <- plotlogbin_cutoff(int_intlight_bin, xl2, yl2,  
                               'Intermediate species', 'Intermediate light', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) 
phighint <- plotlogbin_cutoff(int_highlight_bin, xl2, yl2,  
                                'Intermediate species', 'High light', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) 

plowgap <- plotlogbin_cutoff(gap_lowlight_bin, xl2, yl2,  
                               'gap species', 'Low light', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) 
pintgap <- plotlogbin_cutoff(gap_intlight_bin, xl2, yl2,  
                               'gap species', 'Intermediate light', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) 
phighgap <- plotlogbin_cutoff(gap_highlight_bin, xl2, yl2,  
                                'gap species', 'High light', reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) 

p_labels <- expand.grid(c('shade trees', 'int. trees', 'gap trees'), c('high light', 'int. light', 'low light'))
p_labels <- paste(p_labels[,1], p_labels[,2], sep='\n')

nineplots <- plot_grid(phighshade, phighint, phighgap, pintshade, pintint, pintgap, plowshade, plowint, plowgap, align = 'hv', nrow = 3)

#######################################################################################
# 3. Density scaling, by shade tolerance group, but with PAR as the scaling variable rather than diameter.
# abundance (y) by relative par (x)

# We need to bin again, this time using light received as the binning variable. 
# The density scaling is thus based on light received, not size.

numbins <- 20 # Only for looks
alltree_par_bin <- with(alltreedat, logbin(x=light_received, y=NULL, n = numbins))
shade_par_bin <-  with(shadedat, logbin(x=light_received, y=NULL, n = numbins))
int_par_bin <- with(intdat, logbin(x=light_received, y=NULL, n = numbins))
gap_par_bin <-  with(gapdat, logbin(x=light_received, y=NULL, n = numbins))

# Determine axis limits
raw_numbers <- c(shade_par_bin$bin_value, int_par_bin$bin_value, gap_par_bin$bin_value)/50
y_min <- 10^floor(log10(min(raw_numbers)))
y_max <- max(raw_numbers) * 1.1
y_values <- -6:2
x_values <- 1:5

xl5 <- 'Light received (W)'
yl5 <- expression(paste('Density (individuals ha'^-1,')', sep=''))

pparall <- plotlogbin_cutoff(alltree_par_bin, xl5, yl5,  
                               plottitle = 'All species', plotsubtitle = NULL, reg = TRUE, cutoff = NA, y_min = 5e-7, y_max = 30, x_max = 1e6, y_values = y_values, x_values = x_values)
pparshade <- plotlogbin_cutoff(shade_par_bin, xl5, yl5,  
                               plottitle = 'Shade-tolerant species', plotsubtitle = NULL, reg = TRUE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 4e5, y_values = y_values, x_values = x_values) 
pparint <- plotlogbin_cutoff(int_par_bin, xl5, yl5,  
                               plottitle = 'Intermediate species', plotsubtitle = NULL, reg = TRUE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 4e5, y_values = y_values, x_values = x_values) 
ppargap <- plotlogbin_cutoff(gap_par_bin, xl5, yl5,  
                             plottitle = 'Gap species', plotsubtitle = NULL, reg = TRUE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 4e5, y_values = y_values, x_values = x_values) 


############################################################
# 3b. Same as above, but really a power law should be fit.

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


# Pareto fits
bci_dens_fit_all <- powerlawfit(alltreedat$light_received)
bci_dens_fit_shade <- powerlawfit(shadedat$light_received)
bci_dens_fit_inter <- powerlawfit(intdat$light_received)
bci_dens_fit_gap <- powerlawfit(gapdat$light_received)

xl2 <- 'Light received (W)'
yl2 <- 'Density'

plotbinsandfits(bci_dens_fit_all, alltree_par_bin,
                plottitle = 'All species', xl = xl2, yl = yl2)

plotbinsandfits(bci_dens_fit_shade, shade_par_bin,
                plottitle = 'Shade-tolerant species', xl = xl2, yl = yl2)

plotbinsandfits(bci_dens_fit_inter, int_par_bin,
                plottitle = 'Intermediate species', xl = xl2, yl = yl2)

plotbinsandfits(bci_dens_fit_gap, gap_par_bin,
                plottitle = 'Gap species', xl = xl2, yl = yl2)

# Hypothesis testing

# Density
slopedat <- data.frame(guild = c('all','shade','intermediate','gap'),
                       slope = c(bci_dens_fit_all$alpha, bci_dens_fit_shade$alpha, bci_dens_fit_inter$alpha, bci_dens_fit_gap$alpha),
                       cimin = c(bci_dens_fit_all$boot_ci[1], bci_dens_fit_shade$boot_ci[1], bci_dens_fit_inter$boot_ci[1], bci_dens_fit_gap$boot_ci[1]),
                       cimax = c(bci_dens_fit_all$boot_ci[2], bci_dens_fit_shade$boot_ci[2], bci_dens_fit_inter$boot_ci[2], bci_dens_fit_gap$boot_ci[2]))

ggplot(slopedat, aes(x = guild, y = slope, ymin = cimin, ymax = cimax)) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  ggtitle('Density by diameter')


############################################
# Even more correct way to do it, with cutoffs.

nll_powerlaw <- function(alpha, xmin) {
  C <- (alpha - 1) * ( xmin ^ (alpha - 1) )
  fx <- x ^ ( -alpha )
  px <- C * fx
  -sum(log(px))
} 

# See: http://www.stat.cmu.edu/~cshalizi/2010-10-18-Meetup.pdf
nll_powerlaw_cutoff2 <- function(alpha, xmin, L) {
  C <- (1/L) / (expint::gammainc(1-alpha, xmin/L))
  fx <- ( (x/L)^ -alpha ) * exp(-x/L)
  px <- C * fx
  -sum(log(px))
}

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

# Pareto fits
bci_dens_fit_all <- powerlawfit(alltreedat$light_received)
bci_dens_fit_shade <- powerlawfit(shadedat$light_received)
bci_dens_fit_inter <- powerlawfit(intdat$light_received)
bci_dens_fit_gap <- powerlawfit(gapdat$light_received)

library(stats4)

x <- alltreedat$light_received
fit1 <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
fit2 <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1000), fixed = list(xmin=min(x)), method='BFGS')
aicall1 <- AICc(n = length(x), k = 1, lhat = fit1@min)
aicall2 <- AICc(n = length(x), k = 2, lhat = fit2@min)


x <- subset(alltreedat, tol_wright=='S')$light_received
fit1shade <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
fit2shade <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1000), fixed = list(xmin=min(x)), method='BFGS')
aicshade1 <- AICc(n = length(x), k = 1, lhat = fit1shade@min)
aicshade2 <- AICc(n = length(x), k = 2, lhat = fit2shade@min)

x <- subset(alltreedat, tol_wright=='I')$light_received
fit1int <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
fit2int <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1000), fixed = list(xmin=min(x)), method='BFGS')
aicint1 <- AICc(n = length(x), k = 1, lhat = fit1int@min)
aicint2 <- AICc(n = length(x), k = 2, lhat = fit2int@min)

x <- subset(alltreedat, tol_wright=='G')$light_received
fit1gap <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
fit2gap <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1000), fixed = list(xmin=min(x)), method='BFGS')
aicgap1 <- AICc(n = length(x), k = 1, lhat = fit1gap@min)
aicgap2 <- AICc(n = length(x), k = 2, lhat = fit2gap@min)

# Output coefficients
coefdat <- data.frame(guild = c('all','shade-tolerant','intermediate','gap'),
                      slope1 = c(fit1@coef, fit1shade@coef, fit1int@coef, fit1gap@coef),
                      slope2 = c(fit2@coef[1], fit2shade@coef[1], fit2int@coef[1], fit2gap@coef[1]),
                      cutoff = c(fit2@coef[2], fit2shade@coef[2], fit2int@coef[2], fit2gap@coef[2]))
coefdat$logcutoff <- log10(coefdat$cutoff)
coefdat$deltaaicc <- c(aicall1-aicall2, aicshade1-aicshade2, aicint1-aicint2, aicgap1-aicgap2)

# Bootstrap CIs

boot_mle <- function(xboot, nboot) {
  boot_stats <- list()
  pb <- txtProgressBar(0, nboot, style = 3)
  for (i in 1:nboot) {
    setTxtProgressBar(pb, i)
    x <<- sample(xboot, size = length(xboot), replace = TRUE) # Sets x in global environment.
    fit1_i <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
    fit2_i <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1000), fixed = list(xmin=min(x)), method='BFGS')
    boot_stats[[i]] <- c(slope1=fit1_i@coef[1], slope2=fit2_i@coef[1], cutoff=fit2_i@coef[2])
  }
  close(pb)
  do.call('rbind', boot_stats)
}

boot_shade <- boot_mle(xboot = subset(alltreedat, tol_wright=='S')$light_received, nboot = 99)
boot_int <- boot_mle(xboot = subset(alltreedat, tol_wright=='I')$light_received, nboot = 99)
boot_gap <- boot_mle(xboot = subset(alltreedat, tol_wright=='G')$light_received, nboot = 99)
boot_all <- boot_mle(xboot = alltreedat$light_received, nboot = 99)

shade_ci <- apply(boot_shade, 2, quantile, probs = c(0.025, 0.975))
int_ci <- apply(boot_int, 2, quantile, probs = c(0.025, 0.975))
gap_ci <- apply(boot_gap, 2, quantile, probs = c(0.025, 0.975))
all_ci <- apply(boot_all, 2, quantile, probs = c(0.025, 0.975))

# Make sure confidence interval has converged
nboot <- 99
bootshadeci <- matrix(0, ncol=2, nrow=length(10:nboot))
for (i in 10:nboot) {
  bootshadeci[i-9, ] <- quantile(boot_shade[1:i, 3], probs=c(0.025,0.975))
}

plot(1:90, bootshadeci[,2], type='l', ylim=c(0,10000))
lines(1:90, bootshadeci[,1]) # it converges.

save(boot_shade, boot_int, boot_gap, boot_all, file = 'C:/Users/Q/Dropbox/projects/forestlight/bootout_bylight.r')

load('C:/Users/Q/Dropbox/projects/forestlight/bootout_bylight.r')

# Create the plots.
xl5 <- 'Light received (W)'
yl5 <- expression(paste('Density (individuals ha'^-1,')', sep=''))

# For the 3 groups, concatenate the bin values in individuals per hectare to figure out y axis labels
raw_numbers <- c(shade_par_bin$bin_value, int_par_bin$bin_value, gap_par_bin$bin_value)/50
y_min <- 10^floor(log10(min(raw_numbers)))
y_max <- max(raw_numbers) * 1.1
y_values <- -6:2
x_values <- 1:5


pdall <- plotbinsandfits_cutoff(bci_dens_fit_all, fit2, alltree_par_bin,
                                plottitle = 'All species', xl = xl5, yl = yl5, y_min = 5e-7, y_max = 30, x_min = 1, x_max = 1e6, y_values = y_values, x_values = x_values)
pdshade <- plotbinsandfits_cutoff(bci_dens_fit_shade, fit2shade, shade_par_bin,
                                  plottitle = 'Shade-tolerant species', xl = xl5, yl = yl5, y_min=y_min, y_max=y_max, x_min=1, x_max=1e6, y_values=y_values, x_values = x_values)
pdint <- plotbinsandfits_cutoff(bci_dens_fit_inter, fit2int, int_par_bin,
                                plottitle = 'Intermediate species', xl = xl5, yl = yl5, y_min=y_min, y_max=y_max, x_min=1, x_max=1e6, y_values=y_values, x_values = x_values)
pdgap <- plotbinsandfits_cutoff(bci_dens_fit_gap, fit2gap, gap_par_bin,
                                plottitle = 'Gap species', xl = xl5, yl = yl5, y_min=y_min, y_max=y_max, x_min=1, x_max=1e6, y_values=y_values, x_values = x_values)

# Hypothesis testing: show slopes and possibly AICs.
