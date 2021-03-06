# Production versus light: Fit all functional groups.
# Do with CMDstan 

# New version created 25 June 2018 for "clean" final workflow.
# Edited 12 Feb 2019: Create an additional data dump for fitting scaling models to unnormalized production, using light per area as the scaling variable.

# Data dumps --------------------------------------------------------------

# Load data (changing file path if necessary)
fpdata <- '~/Dropbox/projects/forestlight/stanrdump_final'
load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')

library(dplyr)
library(rstan)


dat90 <- alltree_light_90 %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

dat95 <- alltree_light_95 %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

make_logistic_data <- function(x, n_sub) {
  if (n_sub < nrow(x)) {
    sample_rows <- sample(nrow(x), n_sub)
  } else {
    sample_rows <- 1:nrow(x)
    n_sub <- nrow(x)
  }
  with(x[sample_rows, ],
       list(N = n_sub, x = light_received/crownarea, y = production/crownarea))
}

# Do with no subsampling.
dat90_fg <- dat90 %>%
  group_by(fg) %>%
  do(data = make_logistic_data(., n_sub = nrow(.)))
dat90_all <- make_logistic_data(dat90, n_sub = nrow(dat90))

dat95_fg <- dat95 %>%
  group_by(fg) %>%
  do(data = make_logistic_data(., n_sub = nrow(.)))
dat95_all <- make_logistic_data(dat90, n_sub = nrow(dat95))

# Create stan rdumps
fgnames <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified')

library(purrr)
pmap(dat90_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('dump_light', fg, '1990.r', sep = '_')))))
pmap(dat95_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('dump_light', fg, '1995.r', sep = '_')))))

with(dat90_all, stan_rdump(names(dat90_all), file=file.path(fpdata, 'dump_light_alltree_1990.r')))
with(dat95_all, stan_rdump(names(dat90_all), file=file.path(fpdata, 'dump_light_alltree_1995.r')))



# Data dump with unnormalized production ----------------------------------

make_light_scaling_data <- function(x, n_sub, LL = 1.08, UL = 412.2) {
  if (n_sub < nrow(x)) {
    sample_rows <- sample(nrow(x), n_sub)
  } else {
    sample_rows <- 1:nrow(x)
    n_sub <- nrow(x)
  }
  with(x[sample_rows, ],
       list(N = n_sub, x = light_received/crownarea, y = production, x_min = min(light_received/crownarea), x_max = max(light_received/crownarea), LL = LL, UL = UL))
}

# Do with subsampling.
n_sub <- 25000
set.seed(574)

dat90_fg <- dat90 %>%
  group_by(fg) %>%
  do(data = make_light_scaling_data(., n_sub = n_sub))
dat90_all <- make_light_scaling_data(dat90, n_sub = n_sub)

dat95_fg <- dat95 %>%
  group_by(fg) %>%
  do(data = make_light_scaling_data(., n_sub = n_sub))
dat95_all <- make_light_scaling_data(dat90, n_sub = n_sub)

# Create stan rdumps
fgnames <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified')

library(purrr)
pwalk(dat90_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('ssdump_lightscaling', fg, '1990.r', sep = '_')))))
pwalk(dat95_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('ssdump_lightscaling', fg, '1995.r', sep = '_')))))

with(dat90_all, stan_rdump(names(dat90_all), file=file.path(fpdata, 'ssdump_lightscaling_alltree_1990.r')))
with(dat95_all, stan_rdump(names(dat90_all), file=file.path(fpdata, 'ssdump_lightscaling_alltree_1995.r')))


# Minima, n individuals, and total production -----------------------------

# Get minimum x values and numbers of individuals that are used for plotting.

valall <- data.frame(fg = 'alltree', year = 1995, xmin = with(alltree_light_95, min(light_received/crownarea)), n = nrow(alltree_light_95))
valfg <- alltree_light_95 %>%
  group_by(fg) %>%
  summarize(xmin = min(light_received/crownarea),
            n = n())

min_n <- data.frame(fg = c('alltree', 'fg1','fg2','fg3','fg4','fg5','unclassified'),
                    year = 1995,
                    xmin = c(valall$xmin, valfg$xmin),
                    n = c(valall$n, valfg$n))

write.csv(min_n, '~/Dropbox/projects/forestlight/stanrdump_final/min_n_lighttrees.csv', row.names = FALSE)

# Total production for all trees
library(purrr)
alltreeprod <- sum(alltree_light_95$production)
fgprod <- alltree_light_95 %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg',fg))) %>%
  group_by(fg) %>%
  summarize(production = sum(production))

prodtotals <- data.frame(year = 1995, 
                         rbind(c(fg = 'alltree', production = alltreeprod), fgprod))

write.csv(prodtotals, file = '~/Dropbox/projects/forestlight/stanrdump_final/production_total_lighttrees.csv', row.names = FALSE)


# Data dump: diameter scaling, incoming light energy as y variable --------

make_rawlight_scaling_data <- function(x, n_sub, LL = 1.1, UL = 316) {
  if (n_sub < nrow(x)) {
    sample_rows <- sample(nrow(x), n_sub)
  } else {
    sample_rows <- 1:nrow(x)
    n_sub <- nrow(x)
  }
  with(x[sample_rows, ],
       list(N = n_sub, x = dbh_corr, y = light_received, x_min = min(dbh_corr), x_max = max(dbh_corr), LL = LL, UL = UL))
}

# Do with subsampling.
n_sub <- 25000
set.seed(574)

dat90_fg <- dat90 %>%
  group_by(fg) %>%
  do(data = make_rawlight_scaling_data(., n_sub = n_sub))
dat90_all <- make_rawlight_scaling_data(dat90, n_sub = n_sub)

dat95_fg <- dat95 %>%
  group_by(fg) %>%
  do(data = make_rawlight_scaling_data(., n_sub = n_sub))
dat95_all <- make_rawlight_scaling_data(dat90, n_sub = n_sub)

# Create stan rdumps
fgnames <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified')

library(purrr)
pwalk(dat90_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('ssdump_rawlightscaling', fg, '1990.r', sep = '_')))))
pwalk(dat95_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('ssdump_rawlightscaling', fg, '1995.r', sep = '_')))))

with(dat90_all, stan_rdump(names(dat90_all), file=file.path(fpdata, 'ssdump_rawlightscaling_alltree_1990.r')))
with(dat95_all, stan_rdump(names(dat90_all), file=file.path(fpdata, 'ssdump_rawlightscaling_alltree_1995.r')))

# Total incoming light energy for all trees
library(purrr)
alltreelightrec <- sum(alltree_light_95$light_received)
fglightrec <- alltree_light_95 %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', fg)) %>%
  group_by(fg) %>%
  summarize(light_received = sum(light_received))

lightrectotals <- data.frame(year = 1995, 
                         rbind(c(fg = 'alltree', light_received = alltreelightrec), fglightrec))

write.csv(lightrectotals, file = '~/Dropbox/projects/forestlight/stanrdump_final/lightrec_total.csv', row.names = FALSE)
