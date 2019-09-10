# Try to correct for Jensen's inequality

biasfn <- function(b, v) exp(b*(b-1)*v/2) # Marshall et al. 2019, SI 1

biasfn(0.25, 2)
biasfn(1, 2)
biasfn(1.2, 5)
biasfn(2, 1)
biasfn(0, 1)

biasfn2 <- function(b, sigma) exp(b*(b-1)*(sigma^2)/2)

library(tidyverse)

gdrive_path <- '~/google_drive/ForestLight'
github_path <- '~/Documents/GitHub/forestlight'

# Load raw data
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

# Load binned values for total production
obs_totalprod <- read.csv(file.path(gdrive_path, 'data/data_forplotting_aug2018/obs_totalprod.csv'), stringsAsFactors = FALSE)

#Get the variance in the bins.
source(file.path(github_path, 'code/allfunctions27july.r'))
fp_obs <- file.path(gdrive_path, 'data/data_forplotting_aug2018')
binedgedata <- read.csv(file.path(fp_obs,'obs_dens.csv'),stringsAsFactors = FALSE) %>% filter(fg == 'all', year == 1995) 
area_core <- 42.84

binbounds <- c(binedgedata$bin_min, binedgedata$bin_max[20])

dat <- alltreedat[[3]] %>%
  select(fg, dbh_corr, production) %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = binbounds, include.lowest = TRUE))

dat %>%
  group_by(fg, dbh_bin) %>%
  mutate(var_dbh = var(dbh_corr),
         var_log_dbh = var(log10(dbh_corr)))

# Retrieve fitted values for individual mass, individual production (parameters)
params <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/newpiecewise_paramci_by_fg.csv'), stringsAsFactors = FALSE)

# Retrieve fitted slopes for individual mass, individual production
fittedslopes <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/newpiecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)

# Bin by the fitted slope values
dbh_pred <- exp(seq(log(1.2), log(315), length.out = 101))

dat <- alltreedat[[3]] %>%
  select(fg, dbh_corr, production) %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = dbh_pred, include.lowest = TRUE))

dat %>%
  group_by(fg, dbh_bin) %>%
  mutate(var_dbh = var(dbh_corr),
         var_log_dbh = var(log10(dbh_corr)))


# Calculate bias correction factor: exp(slope * (slope-1) * variance^2) for each bin.


# Multiply individual mass * individual production * correction factor.
