# Relationship between diameter and diameter growth per time (instead of mass growth per time)
# QDR / Forestlight / 24 July 2019

library(tidyverse)
library(rstan)

#Q Dog
gdrive_path <- '~/google_drive/ForestLight'
github_path <- '~/Documents/GitHub/forestlight'

#Grady_2
gdrive_path <- '/Users/johngrady/Google Drive/ForestLight'
github_path <- '/Users/johngrady/Documents/GitHub/forestlight'

#Grady
gdrive_path <- '/Users/jgradym/Google Drive/ForestLight'
github_path <- '/Users/jgradym/Documents/GitHub/forestlight'


# Load data ---------------------------------------------------------------

# Load raw data
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

# We need to annualize dbh increment to a 1 year difference since the raw dbh increment is 5 year difference.

annual_increment <- function(dbh_old, dbh_new, census_interval = 5, new_interval = 1){
  rate <-  (dbh_new / dbh_old)^(1/census_interval) - 1
  dbh_oldplus1year <- dbh_old * (1 + rate)^new_interval
  return(dbh_oldplus1year - dbh_old)
}

dat <- alltreedat[[3]] %>%
  select(sp, fg, dbh_corr, dbhlastcensus) %>%
  mutate(dbh_increment = annual_increment(dbhlastcensus, dbh_corr),
         fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))


# Quick diagnostic plot ---------------------------------------------------

ggplot(dat, aes(x = dbh_corr, y = dbh_increment)) +
  geom_hex() +
  scale_x_log10() + 
  scale_y_log10() +
  facet_wrap(~ fg)

# Fit 1 segment and 2 segment power law -----------------------------------

# Compile stan models
productionmodel1 <- stan_model(file.path(github_path, 'stan/piecewise_workflow/production1.stan'), model_name = 'production_model_1')
productionmodel2 <- stan_model(file.path(github_path, 'stan/piecewise_workflow/production2.stan'), model_name = 'production_model_2')

options(mc.cores = 3)

# Create data for each functional group
data_by_fg <- dat %>%
  group_by(fg) %>%
  do(dat = list(N = nrow(.),
                x = .$dbh_corr,
                y = .$dbh_increment,
                x_min = min(.$dbh_corr),
                x_max = max(.$dbh_corr)))

data_all <- dat %>%
  mutate(fg = 'all') %>%
  group_by(fg) %>%
  do(dat = list(N = nrow(.),
                x = .$dbh_corr,
                y = .$dbh_increment,
                x_min = min(.$dbh_corr),
                x_max = max(.$dbh_corr)))

data_list <- c(data_by_fg$dat, data_all$dat)

# Fit 1 segment production model

# test 
fit1_fg1 <- sampling(productionmodel1, data = data_list[[1]], chains = 3, iter = 6000, warmup = 5000, seed = 1111, algorithm = 'HMC', control = list(max_treedepth = 20, adapt_delta = 0.9), pars = c('beta0', 'beta1', 'sigma'))
fit2_fg1 <- sampling(productionmodel2, data = data_list[[1]], chains = 3, iter = 6000, warmup = 5000, seed = 1111, algorithm = 'HMC', control = list(max_treedepth = 20, adapt_delta = 0.9), pars =  c('x0', 'beta0', 'beta1_low', 'beta1_high', 'delta', 'sigma'))

# This works fine but takes too long so we will move it to CMDSTAN on cluster.


# Create cmdstan data dumps -----------------------------------------------

walk2(data_list, c(data_by_fg$fg, 'alltree'), ~ with(.x, stan_rdump(list = names(.x), file = file.path('~/Dropbox/projects/forestlight/stanrdump_final', paste0('dump_diamgrowthscaling_',.y,'_1995.r')))))

     