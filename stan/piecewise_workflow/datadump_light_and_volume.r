# Dumps for individual light and individual crown volume
# For reproducing fits in Figure 4
# QDR / Forestlight / 14 June 2019

# Load data (changing file path if necessary)
fpdata <- '~/Dropbox/projects/forestlight/stanrdump_final'
load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')

fgnames <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified')

library(dplyr)
library(rstan)
library(purrr)

dat90 <- alltree_light_90 %>%
  select(dbh_corr, production, light_received, crownarea, crownvolume, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

dat95 <- alltree_light_95 %>%
  select(dbh_corr, production, light_received, crownarea, crownvolume, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

# Data dump: diameter scaling, incoming light energy as y variable --------

# No subsampling.
make_rawlight_scaling_data <- function(x, LL = 1.1, UL = 316) {
  with(x,
       list(N = nrow(x), x = dbh_corr, y = light_received, x_min = min(dbh_corr), x_max = max(dbh_corr), LL = LL, UL = UL))
}

dat90_fg <- dat90 %>%
  group_by(fg) %>%
  do(data = make_rawlight_scaling_data(.))
dat90_all <- make_rawlight_scaling_data(dat90)

dat95_fg <- dat95 %>%
  group_by(fg) %>%
  do(data = make_rawlight_scaling_data(.))
dat95_all <- make_rawlight_scaling_data(dat95)

# Create stan rdumps

pwalk(dat90_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('dump_rawlightscaling', fg, '1990.r', sep = '_')))))
pwalk(dat95_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('dump_rawlightscaling', fg, '1995.r', sep = '_')))))

with(dat90_all, stan_rdump(names(dat90_all), file=file.path(fpdata, 'dump_rawlightscaling_alltree_1990.r')))
with(dat95_all, stan_rdump(names(dat95_all), file=file.path(fpdata, 'dump_rawlightscaling_alltree_1995.r')))

# Total incoming light energy for all trees
alltreelightrec <- sum(alltree_light_95$light_received)
fglightrec <- alltree_light_95 %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', fg)) %>%
  group_by(fg) %>%
  summarize(light_received = sum(light_received))

lightrectotals <- data.frame(year = 1995, 
                             rbind(c(fg = 'alltree', light_received = alltreelightrec), fglightrec))

write.csv(lightrectotals, file = '~/Dropbox/projects/forestlight/stanrdump_final/lightrec_total.csv', row.names = FALSE)


# Data dump: diameter scaling, individual crown volume as y variab --------

# No subsampling.
make_volume_scaling_data <- function(x, LL = 1.1, UL = 316) {
  with(x,
       list(N = nrow(x), x = dbh_corr, y = crownvolume, x_min = min(dbh_corr), x_max = max(dbh_corr), LL = LL, UL = UL))
}

voldat90_fg <- dat90 %>%
  group_by(fg) %>%
  do(data = make_volume_scaling_data(.))
voldat90_all <- make_volume_scaling_data(dat90)

voldat95_fg <- dat95 %>%
  group_by(fg) %>%
  do(data = make_volume_scaling_data(.))
voldat95_all <- make_volume_scaling_data(dat95)

# Create stan rdumps

pwalk(voldat90_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('dump_volumescaling', fg, '1990.r', sep = '_')))))
pwalk(voldat95_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('dump_volumescaling', fg, '1995.r', sep = '_')))))

with(voldat90_all, stan_rdump(names(dat90_all), file=file.path(fpdata, 'dump_volumescaling_alltree_1990.r')))
with(voldat95_all, stan_rdump(names(dat95_all), file=file.path(fpdata, 'dump_volumescaling_alltree_1995.r')))

# Total crown volume for all trees
alltreecrownvol <- sum(alltree_light_95$crownvolume)
fgcrownvol <- alltree_light_95 %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
  group_by(fg) %>%
  summarize(crownvolume = sum(crownvolume))

crownvolumetotals <- data.frame(year = 1995, 
                             rbind(c(fg = 'alltree', crownvolume = alltreecrownvol), fgcrownvol))

write.csv(crownvolumetotals, file = '~/Dropbox/projects/forestlight/stanrdump_final/crownvol_total.csv', row.names = FALSE)
