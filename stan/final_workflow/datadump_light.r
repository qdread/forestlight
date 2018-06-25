# Production versus light: Fit all functional groups.
# Do with CMDstan 

# New version created 25 June for "clean" final workflow.

# Data dumps --------------------------------------------------------------

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/Dropbox/projects/forestlight/stanrdump_final'
load('C:/Users/Q/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')

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


