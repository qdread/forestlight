# Fit density and production with Bayesian method for shade trees and gap trees for each year.
# This generates a predictive interval for each relationship.
# For production, fit powerlaw and powerlaw * exponential. (easy)
# For density, fit powerlaw, powerlaw * exponential, and piecewise powerlaw+exponential. (hard)

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_04oct'
load(file.path(fpdata, 'rawdataobj.r'))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

estimate_model_pred <- function(x, y, x_pred, stanmodel) {
  N <- length(x)
  N_pred <- length(x_pred)
  data <- list(N = N, N_pred = N_pred, x = x, y = y, x_pred = x_pred)
  fit <- sampling(stanmodel, data = data, chains = 3, iter = 2000, warmup = 1000)
  return(fit)
}

summ2plotdata <- function(summ, rawdata, xname) {
  yrows <- grep('y_pred', dimnames(summ$summary)[[1]])
  res <- data.frame(x = seq(min(rawdata[,xname]),max(rawdata[,xname]),length.out = length(yrows)),
                    pred_median = summ$summary[yrows,'50%'],
                    pred_025 = summ$summary[yrows,'2.5%'],
                    pred_975 = summ$summary[yrows,'97.5%'])
  names(res)[1] <- xname
  res
}

plotplawwithdata <- function(fitdata, rawdata, xvar, yvar, xname, yname) {
  require(cowplot)
  ggplot(fitdata) +
    geom_hex(aes_string(x=xvar, y=xvar), data=rawdata) +
    scale_fill_gradient(low = 'gray95', high='gray10') +
    scale_color_manual(values = c('indianred','dodgerblue')) +
    geom_line(aes_string(x=xvar, y='pred_median', group = 'model', color = 'model'), size=1) + 
    geom_line(aes_string(x=xvar, y='pred_025', group = 'model', color = 'model'), linetype = 'dotted', size=0.75) +
    geom_line(aes_string(x=xvar, y='pred_975', group = 'model', color = 'model'), linetype = 'dotted', size=0.75) +
    scale_x_log10(name = xname) + scale_y_log10(name = yname)
}

plotplawnodata <- function(fitdata, xvar, xname, yname) {
  require(cowplot)
  ggplot(fitdata) +
    scale_color_manual(values = c('indianred','dodgerblue')) +
    scale_fill_manual(values = c('indianred','dodgerblue')) +
    geom_ribbon(aes_string(x=xvar, ymin='pred_025', ymax='pred_975', group = 'model', fill = 'model'), alpha=0.25) + 
    geom_line(aes_string(x=xvar, y='pred_median', group = 'model', color = 'model')) + 
    scale_x_log10(name = xname) + scale_y_log10(name = yname)
}

powerlaw_trans_model <- stan_model(file = 'stan/model_powerlaw_logtrans.stan')
powerlawexp_trans_model <- stan_model(file = 'stan/model_powerlawexp_logtrans.stan')

# Fit to gap 2010
gap_fit_powerlaw_trans <- with(gapdat[[6]], estimate_model_pred(x = dbh_corr, y = production, x_pred = seq(min(dbh_corr),max(dbh_corr),length.out = 10), stanmodel = powerlaw_trans_model))
gap_fit_powerlawexp_trans <- with(gapdat[[6]], estimate_model_pred(x = dbh_corr, y = production, x_pred = seq(min(dbh_corr),max(dbh_corr),length.out = 50), stanmodel = powerlawexp_trans_model))

gap_plaw_trans_summary <- summary(gap_fit_powerlaw_trans)
gap_exp_trans_summary <- summary(gap_fit_powerlawexp_trans)

# Generate plot data for gap 2010
gap_plaw_plotdat <- summ2plotdata(summ = gap_plaw_trans_summary, rawdata = gapdat[[6]], xname = 'dbh_corr')
gap_exp_plotdat <- summ2plotdata(summ = gap_exp_trans_summary, rawdata = gapdat[[6]], xname = 'dbh_corr')

gap_all_trans_sum_dat <- rbind(cbind(model = 'powerlaw', gap_plaw_plotdat),
                               cbind(model = 'powerlaw_exp', gap_exp_plotdat))

# Plot gap 2010
plotplawwithdata(fitdata = gap_all_trans_sum_dat, rawdata = gapdat[[6]], xvar = 'dbh_corr', yvar = 'production', xname = 'Diameter (cm)', yname = 'Production (kg y-1)')
plotplawnodata(fitdata = gap_all_trans_sum_dat, xvar = 'dbh_corr', xname = 'Diameter (cm)', yname = 'Production (kg y-1)')


# Fit to shade 2010
shade_fit_powerlaw_trans <- with(shadedat[[6]], estimate_model_pred(x = dbh_corr, y = production, x_pred = seq(min(dbh_corr),max(dbh_corr),length.out = 10), stanmodel = powerlaw_trans_model))
shade_fit_powerlawexp_trans <- with(shadedat[[6]], estimate_model_pred(x = dbh_corr, y = production, x_pred = seq(min(dbh_corr),max(dbh_corr),length.out = 50), stanmodel = powerlawexp_trans_model))

shade_plaw_trans_summary <- summary(shade_fit_powerlaw_trans)
shade_exp_trans_summary <- summary(shade_fit_powerlawexp_trans)

# Generate plot data for shade 2010
shade_plaw_plotdat <- summ2plotdata(summ = shade_plaw_trans_summary, rawdata = shadedat[[6]], xname = 'dbh_corr')
shade_exp_plotdat <- summ2plotdata(summ = shade_exp_trans_summary, rawdata = shadedat[[6]], xname = 'dbh_corr')

shade_all_trans_sum_dat <- rbind(cbind(model = 'powerlaw', shade_plaw_plotdat),
                               cbind(model = 'powerlaw_exp', shade_exp_plotdat))

# Plot shade 2010
plotplawwithdata(fitdata = shade_all_trans_sum_dat, rawdata = shadedat[[6]], xvar = 'dbh_corr', yvar = 'production', xname = 'Diameter (cm)', yname = 'Production (kg y-1)')
plotplawnodata(fitdata = shade_all_trans_sum_dat, xvar = 'dbh_corr', xname = 'Diameter (cm)', yname = 'Production (kg y-1)')
