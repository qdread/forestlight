# Mass & diameter energy equivalence with some simulated data
# Edited 17 april to sim with different bin widths

# Details of simulation
rmin <- 1
rmax <- 100
binwidth <- 1

# Constants
n0 <- 100000
b0 <- 0.1
k <- 1

simdat <- list()

for (binwidth in c(0.1, 1, 10)) {

r <- seq(from=rmin, to=rmax, by=binwidth)
n <- round(n0 * r^-2)
b <- b0 * r^2
m <- k * r^(8/3)

radii <- rep(r, n)
masses <- rep(m, n)
productions <- rep(b, n)

# Energy equivalence by radius
plot(r, n*b, ylim=c(0.1, max(n*b)), log='xy')

# Bin linearly by mass into 100 bins
masscuts <- seq(0, max(m), length.out=length(r) + 1)
masslower <- masscuts[-(length(r) + 1)]
massupper <- masscuts[-1]
massupper[length(r)] <- massupper[length(r)] + 1

mb <- rep(0, length(r))

for (i in 1:length(r)) {
  mb[i] <- which(m[i] >= masslower & m[i] < massupper)
}

mbfact <- as.numeric(factor(mb, labels=1:length(unique(mb))))

massbins <- rep(mbfact, n)
prodbinnedbymass <- tapply(productions, massbins, sum)
densbinnedbymass <- table(massbins)
densbinnedbyradius <- table(radii)

simdat[[length(simdat) + 1]] <- data.frame(m = m[unique(mb)], n = as.numeric(densbinnedbymass), binwidth=binwidth)

}

simdat <- do.call('rbind', simdat)

ggplot(simdat, aes(x=m, y=n, group=binwidth, color=factor(binwidth))) +
  geom_line() +
  geom_abline(slope=-3/4, intercept=5, color='black',size=1.5) +
  theme_bw() +
  scale_x_log10() + scale_y_log10()





plot(m[unique(mb)], prodbinnedbymass, log='xy')

# Added 11 April: density, in addition to production.

pdf('C:/Users/Q/Google Drive/ForestLight/figs/dens_by_rad_and_mass.pdf', height=5, width=5)
plot(r, n, log='xy', main='Density by radius\n-2 slope line plotted') # Density by radius (exact)
abline(a=5, b=-2, col='red')
plot(m, n, log='xy', main='Density by mass\n-3/4 slope line plotted')
abline(a=5, b=-3/4, col='red')
dev.off()

# Correct binning

# Bin by radius and by mass
pdf('C:/Users/Q/Google Drive/ForestLight/figs/dens_by_rad_and_mass.pdf', height=5, width=5)
plot(r, as.numeric(densbinnedbyradius), log='xy', main='Density by radius\n-2 slope line plotted')
abline(a=5, b=-2, col='red')
plot(m[unique(mb)], as.numeric(densbinnedbymass), log='xy', main='Density by mass\n-3/4 slope line plotted', type='l')
plot(m[unique(mb)], as.numeric(densbinnedbymass1), log='xy', main='Density by mass\n-3/4 slope line plotted', type='l')

abline(a=5, b=-3/4, col='red')
dev.off()

# Pareto fit to density data.
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
  #n_boot <- 1499
  #n_burn <- 500
  #pl_boot <- bootstrap(m = pl_dat, xmins = pl_dat$getXmin(), no_of_sims = n_boot)
  #boot_ci <- quantile(pl_boot$bootstraps$pars[-(1:n_burn)], probs = c(0.025, 0.975))
  
  return(list(plotdat = plotdat, 
              plfit = plfit_dat, 
              lognormfit = lognormfit_dat, 
              plpdf = data.frame(x = plfit_dat$x, y = pl_pdf),
              lognormpdf = data.frame(x = lognormfit_dat$x, y = lognorm_pdf),
              xmin = xmin_pl, 
              alpha = pars_pl$pars,
              xmin_lognorm = xmin_lognorm,
              pars_lognorm = pars_lognorm$pars
              #boot_ci = as.numeric(boot_ci)
              ))
}

pfitrad <- powerlawfit(dat = radii)
pfitmass <- powerlawfit(dat = masses)


# Power law fit using a different package
library(igraph)

power.law.fit(x = radii)
power.law.fit(x = masses)

# Use log binning algorithm on the data
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

prod_logbin_mass <- logbin(x = masses, y = productions, n = 99)
prod_logbin_radius <- logbin(x = radii, y = productions, n = 99)


# Better simulation.

# Sample r directly from a power law.
library(actuar)

set.seed(8720604)

x <- rpareto1(1e5, 1, 1) # one million trees, yo. 
r <- sort(x)

b0 <- 0.1
k <- 1

b <- b0 * r^2
m <- k * r^(8/3)

fit_r <- powerlawfit(r) # Returns alpha=2
fit_m <- powerlawfit(m) # Returns alpha=1.375
fit_b <- powerlawfit(b)

prodbin_r <- logbin(r, b, 50)
prodbin_m <- logbin(m, b, 50)

# Bin production by radius.
qcut <- function(x, n_bins) {
  cuts <- seq(0, max(x), length.out=n_bins + 1)
  binlower <- cuts[-length(cuts)]
  binupper <- cuts[-1]
  binupper[length(binupper)] <- binupper[length(binupper)] + 1
  
  bins <- rep(0, length(x))
  
  for (i in 1:length(x)) {
    bins[i] <- which(x[i] >= binlower & x[i] < binupper)
  }
  
  as.numeric(factor(bins, labels=1:length(unique(bins))))
  
}

rbin <- qcut(r, 100)
mbin <- qcut(m, 100)

b_r <- tapply(b, rbin, sum)
b_m <- tapply(b, mbin, sum)

