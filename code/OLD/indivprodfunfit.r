# Function fitting code for individual production
# Version 1: powerlaw multiplied by exponential
# Version 2: piecewise: first fit powerlaw, then fit a cutoff point with an exponential below the cutoff.
# Modify this from Caroline's code
# QDR 10 Oct 2017

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_04oct'
load(file.path(fpdata, 'rawdataobj.r'))


source('~/GitHub/FunctionFitting/PowerLawFit_etc_160509.r')

library(dplyr)

dat <- shadedat[[6]] %>%
  transmute(logprod = log10(production), logdbh = log10(dbh_corr)) %>%
  arrange(logdbh, logprod)

# Fit line to log(production) ~ log(diameter)
shade2010prodlm <- lm(I(log10(production)) ~ I(log10(dbh_corr)), data = shadedat[[6]])

library(cowplot)

ggplot(shadedat[[6]], aes(x=dbh_corr,y=production)) +
  geom_hex() + scale_x_log10() + scale_y_log10() + scale_fill_gradient(low='gray95',high='black')

nll_powerlaw_cutoff2 <- function(alpha, xmin, L) {
  C <- (1/L) / (gsl::gamma_inc(1-alpha, xmin/L))
  fx <- ( (x/L)^ -alpha ) * exp(-x/L)
  px <- C * fx
  -sum(log(px))
}

# Fit exponential.
shade2010prodexp <- nls(logprod ~ I((2.316 * logdbh + (-1.653)) * (a * exp(b * logdbh) + c)), data = dat, start = list(a = -1, b = -1, c = 1))

shade2010prodexp_fitall <- nls(logprod ~ I((a1 * logdbh + b1) * (a * exp(b * logdbh) + c)), data = dat, start = list(a = 0.2, b = -10, c = 1, a1=2, b1=-1.5))

# Try to fit all parameters with better algorithm
library(minpack.lm)

predict_fn <- function(pars, xval) (pars$a1 * xval + pars$b1) * (pars$a * exp(pars$b * xval) + pars$c)
resid_fn <- function(p, obs, x) obs - predict_fn(pars=p, xval=x)

start_values <- list(a = 0.2, b = -10, c = 1, a1=2, b1=-1.5)

shade2010prodexp_fitall <- nls.lm(par = start_values, fn = resid_fn, obs = dat$logprod, x = dat$logdbh)
  


# Test starting values
expf <- function(x, a, b, c) a * exp(b * x) + c

curve(expf(x, -7.7, -4.5, 0.42))

# Powerlaw times exponential:
powerlaw_exp_fn <- function(x, pars1, pars2) (pars1[1] + pars1[2] * log10(x)) * (pars2['a'] * exp(pars2['b'] * log10(x)) + pars2['c'])
powerlaw_exp_fn <- function(x, pars) (pars$b1 + pars$a1 * log10(x)) * (pars$a * exp(pars$b * log10(x)) + pars$c)

# nls way
shadedat[[6]] %>%
  ggplot(aes(x = dbh_corr, y = production)) +
  geom_hex() +
  stat_function(fun = powerlaw_exp_fn, 
                args = list(pars1 = shade2010prodlm$coefficients, pars2 = shade2010prodexp$m$getPars()),
                geom = 'line',
                color = 'red') +
  scale_fill_gradient(low='gray95', high='black') +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Production per unit light') +
  panel_border(colour = 'black')

# minpack nls.lm way
shadedat[[6]] %>%
  ggplot(aes(x = dbh_corr, y = production)) +
  geom_hex() +
  stat_function(fun = powerlaw_exp_fn, 
                args = list(pars = shade2010prodexp_fitall$par),
                geom = 'line',
                color = 'red') +
  scale_fill_gradient(low='gray95', high='black') +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Production per unit light') +
  panel_border(colour = 'black')

####
# Katharine Mullen's function to get log likelihood from nls.lm obj.
logLik.nls.lm <- function(object, REML = FALSE, ...) 
{ 
  res <- object$fvec 
  N <- length(res) 
  val <-  -N * (log(2 * pi) + 1 - log(N) + log(sum(res^2)))/2 
  ## the formula here corresponds to estimating sigma^2. 
  attr(val, "df") <- 1L + length(coef(object)) 
  attr(val, "nobs") <- attr(val, "nall") <- N 
  class(val) <- "logLik" 
  val 
} 
2 * (5) - 2 * logLik.nls.lm(shade2010prodexp_fitall)



# Function fitting for all censuses and guilds ----------------------------

# Individual production
fit_production <- function(dat) {
  require(dplyr)
  require(minpack.lm)
  
  # Log-transform
  dat <- dat %>%
    transmute(logprod = log10(production), logdbh = log10(dbh_corr)) %>%
    arrange(logdbh, logprod)
  
  # Straight power law fit
  powerlaw_fit <- lm(logprod ~ logdbh, data = dat)
  
  # Power law times exponential fit
  predict_fn <- function(pars, xval) (pars$a1 * xval + pars$b1) * (pars$a * exp(pars$b * xval) + pars$c)
  resid_fn <- function(p, obs, x) obs - predict_fn(pars=p, xval=x)
  
  start_values <- list(a = 0.2, b = -10, c = 1, a1=powerlaw_fit$coef[2], b1=powerlaw_fit$coef[1]) # Might change.
  
  powerlawexp_fit <- nls.lm(par = start_values, fn = resid_fn, obs = dat$logprod, x = dat$logdbh)
  
  return(list(powerlaw_fit = powerlaw_fit, powerlawexp_fit = powerlawexp_fit))
}

alltree_fits <- lapply(alltreedat, fit_production)
shade_fits <- lapply(shadedat, fit_production)
gap_fits <- lapply(gapdat, fit_production)

# Calculate AIC
getaics <- function(fits) {
  return(c(AIC_powerlaw = AIC(fits$powerlaw_fit),
           AIC_powerlaw_exp = as.numeric(10 - 2*logLik.nls.lm(fits$powerlawexp_fit))))
}

alltree_aic <- lapply(alltree_fits, getaics)
shade_aic <- lapply(shade_fits, getaics)
gap_aic <- lapply(gap_fits, getaics)

# Plot all censuses and guilds
plot_production <- function(dat) {
  require(cowplot)
  powerlaw_exp_fn <- function(x, pars) (pars$b1 + pars$a1 * log10(x)) * (pars$a * exp(pars$b * log10(x)) + pars$c)
  powerlaw_fn <- function(x, pars) pars$b1 + pars$a1 * log10(x)
  
  fit <- fit_production(dat)
  
  dat %>%
    ggplot(aes(x = dbh_corr, y = production)) +
    geom_hex() +
    stat_function(fun = powerlaw_fn,
                  args = list(pars = fit$powerlawexp_fit$par),
                  geom = 'line',
                  color = 'blue',
                  linetype = 'dotted',
                  size = 1.5) +
    stat_function(fun = powerlaw_exp_fn, 
                  args = list(pars = fit$powerlawexp_fit$par),
                  geom = 'line',
                  color = 'red',
                  size = 1.5) +
    scale_fill_gradient(low='gray95', high='black') +
    scale_x_log10(name = 'Diameter (cm)') + 
    scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
    panel_border(colour = 'black')
}


alltree_indivprod_plots <- lapply(alltreedat, plot_production)
shade_indivprod_plots <- lapply(shadedat, plot_production)
gap_indivprod_plots <- lapply(gapdat, plot_production)

