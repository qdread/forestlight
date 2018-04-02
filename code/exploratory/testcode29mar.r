truncated_weibull <- function(x, m, n, ll, ul) {
  trunc_pts <- pweibull(q = c(ll,ul), shape = m, scale = n)
  out <- dweibull(x = x, shape = m, scale = n)
  out <- out/diff(trunc_pts)
}

integrate(truncated_weibull, lower=1.2, upper=1.2687, m=0.4213, n=0.6076, ll=1.1, ul=316)
integrate(truncated_weibull, lower=1.2687, upper=1.341, m=0.4213, n=0.6076, ll=1.1, ul=316)

plot_dens(year_to_plot = 1995,
          fg_names = c('all'),
          model_fit = 'weibull',
          y_limits = c(0.0001, 1200),
          y_breaks = c(0.001, 0.1, 10))
plot_prod(year_to_plot = 1995,
          fg_names = c('all'),
          model_fit = 'powerlaw',
          y_limits = c(0.01, 10000),
          y_breaks = c(0.1, 10, 1000))
p1 <- plot_totalprod(year_to_plot = 1995,
               fg_names = c('all'),
               model_fit_density = 'weibull', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 500),
               y_breaks = c(0.1, 1, 10, 100))

dens1995 <- obs_dens %>% filter(year==1995, fg=='all')
prod1995 <- obs_indivprod %>% filter(year==1995, fg=='all')
totalprod1995 <- obs_totalprod %>% filter(year==1995, fg=='all')

production_product <- dens1995$bin_value * prod1995$median

dat <- data.frame(dbh = totalprod1995$bin_midpoint, method = rep(c('binning', 'product'), each=20), production = c(totalprod1995$bin_value, production_product))
ggplot(dat, aes(x=dbh,y=production,color=method)) +
  geom_point() +
  scale_x_log10() + scale_y_log10(name='production per hectare') +
  theme_bw()
sum(totalprod1995$bin_value)  
sum(production_product)
sum(alltreedat[[3]]$production)/42.84

predtotalprod1995 <- pred_totalprod %>% filter(year==1995, fg=='all', dens_model=='weibull')
mids <- function(a) a[-length(a)] + diff(a)/2
mids(predtotalprod1995$q50) / diff(predtotalprod1995$dbh)

# Mess with CI locally ----------------------------------------------------


fp <- 'C:/Users/Q/Dropbox/projects/forestlight/'

source('stan/extract_ci_stan_corrected.r')

library(purrr)
library(tidyr)
library(dplyr)
library(rstan)


min_n <- read.csv(file.path(fp, 'stanoutput/min_n.csv'), stringsAsFactors = FALSE)
dbh_pred <- exp(seq(log(1.1), log(316), length.out = 101))

load(file.path(fp, 'fit_wpow_alltree_1995.r'))
ci1 <- dens_prod_ci(fit, dbh_pred, dens_form = 'weibull', prod_form = 'powerlaw', x_min = 1.1, n_indiv = 114058)
ci1prod <- ci1 %>%
  filter(variable=='total_production') %>%
  mutate_at(vars(starts_with('q')), funs(./42.84))

p1 + geom_line(data=ci1prod, aes(x=dbh,y=q50))
