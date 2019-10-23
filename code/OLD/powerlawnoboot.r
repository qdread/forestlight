# Power law without bootstrap

powerlawfit <- function(dat, doboot=T) {
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
  if (doboot) {
    n_boot <- 1499
    n_burn <- 500
    pl_boot <- bootstrap(m = pl_dat, xmins = pl_dat$getXmin(), no_of_sims = n_boot, threads = 3)
    boot_ci <- quantile(pl_boot$bootstraps$pars[-(1:n_burn)], probs = c(0.025, 0.975))
  }
  else {
    boot_ci <- c(NA,NA)
  }
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