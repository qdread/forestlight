# Total production plot for some functional groups for Weibull and power law in 1995
plot_totalprod(year_to_plot = 1995,
               fg_names = c('all'),
               model_fit_density = 'weibull', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 400),
               y_breaks = c(0.1, 1, 10, 100))

plot_dens(year_to_plot = 1995,
          fg_names = c('all'),
          model_fit = 'weibull',
          y_limits = c(0.0001, 1300),
          y_breaks = c(0.001, 0.1, 10, 1000))


# Does integral of total production equal sum of bins

dbh_pred <- exp(seq(log(1.2),log(315),length.out=50))
dbh_binwidth <- diff(dbh_pred)
sum(mids(tppred$q50) * dbh_binwidth)

dbh_obs <- c(tp1$bin_min[1], tp1$bin_max)
dbh_obs_binwidth <- diff(dbh_obs)
sum(tp1$bin_value * dbh_obs_binwidth)

sum(alltreedat[[3]]$production)/area_core # Total production agrees with the observed production values.
# The density times production is systematically lower.

plot_totalprod(year_to_plot = 1995,
               fg_names = c('all'),
               model_fit_density = 'weibull', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 400),
               y_breaks = c(0.1, 1, 10, 100))

newtotalprod <- logbin_setedges(alltreedat[[3]]$dbh_corr, alltreedat[[3]]$production, 
                                edges = data.frame(bin_min = dbh_pred[-length(dbh_pred)],
                                                   bin_max = dbh_pred[-1],
                                                   bin_midpoint = mids(dbh_pred)))

plot_totalprod(year_to_plot = 1995,
               fg_names = c('all'),
               model_fit_density = 'weibull', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 400),
               y_breaks = c(0.1, 1, 10, 100),
               obsdat = cbind(fg='all',year=1995, newtotalprod) %>% mutate(bin_value = bin_value/area_core))

myrpareto <- function(n, xmax, xmin, alpha) {
  u <- runif(n, xmin, xmax)
  exp(log(xmin) - (1-u)/alpha)
}

myrpareto(100, 316, 1.1, 1)
