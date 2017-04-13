# Mass & diameter energy equivalence with some simulated data

# Details of simulation
rmin <- 1
rmax <- 100
binwidth <- 0.1

# Constants
n0 <- 100000
b0 <- 0.1
k <- 1

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

