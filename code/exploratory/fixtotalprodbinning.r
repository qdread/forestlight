library(dplyr)
d1 <- obs_dens %>%
  filter(fg == 'all', year == 1995)
p1 <- obs_indivprod %>%
  filter(fg == 'all', year == 1995)
tp <- d1$bin_value * p1$median
tp1 <- obs_totalprod %>%
  filter(fg == 'all', year == 1995)

tp1$value2 <- tp

ggplot(tp1, aes(x=bin_midpoint,y=bin_value)) + geom_point() + 
  geom_point(aes(y=value2),color='red') + 
  scale_x_log10() + scale_y_log10()

# Create bin value by average production * number of individuals / bin width
# Number of individuals is integral of density in the bin.

d_binwidth <- d1$bin_max - d1$bin_min
d_nbin <- d1$bin_value * d_binwidth
sumprod <- p1$median * d_nbin

dpred <- pred_dens %>%
  filter(fg == 'all', dens_model == 'weibull', prod_model == 'powerlaw', year == 1995)
ppred <- pred_indivprod %>%
  filter(fg == 'all', dens_model == 'weibull', prod_model == 'powerlaw', year == 1995)
tppred <- pred_totalprod %>%
  filter(fg == 'all', dens_model == 'weibull', prod_model == 'powerlaw', year == 1995)

allpred <- data.frame(dbh = dpred[,'dbh'], dens = dpred[,'q50'], prod = ppred[,'q50'], totalprod = tppred[,'q50'])
allpred$product <- allpred$dens * allpred$prod


# All trees, weibull, powerlaw, 1995
pars <- structure(c(0.41908453964862, 0.582521657743589, 3.85282689498638,
                    2.33520958923449, 0.751139810725117, -296715.800756511), .Names = c("m",
                                                                                        "n", "beta0", "beta1", "sigma", "lp__"))
pars

dbhs <- exp(seq(log(1.2), log(315), length.out = 50))

ll <- 1.1
ul <- 316
trunc_proportion <- sum(1 - pweibull(q = c(ll,ul), shape = pars['m'], scale = pars['n']))
dens_pred <- sapply(dbhs, dweibull, shape = pars['m'], scale = pars['n']) 
dens_pred <- dens_pred/trunc_proportion
n_ind <- 114058
area_core <- 42.84
dens_pred_inds <- dens_pred*n_ind/area_core

cutoffs <- pweibull(q = c(ll,ul), shape = pars['m'], scale = pars['n'])

weibcdf <- sapply(dbhs, pweibull, shape = pars['m'], scale = pars['n'])
weibpdf <- sapply(dbhs, dweibull, shape = pars['m'], scale = pars['n'])
weibcdf <- (weibcdf - cutoffs[1])/(cutoffs[2] - cutoffs[1]) # Cumulative distribution function of truncated

n_ind_bin <- diff(weibcdf)/diff(dbhs) * n_ind/area_core

mids <- function(a) a[-length(a)] + diff(a)/2
totprodpred_new <- n_ind_bin * mids(ppred$q50)

x <- data.frame(dbh = mids(dbhs), newprod = totprodpred_new, oldprod = mids(tppred$q50))
ggplot(x, aes(x=dbh)) + geom_point(aes(y=newprod)) + geom_point(aes(y=oldprod),color='red') + scale_x_log10() + scale_y_log10()
