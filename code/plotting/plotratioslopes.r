# Plot ratio slopes

slopes_light <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/lightpiecewise/lightpiecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE)
slopes_diam <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE)

library(tidyverse)

ratio_slopes_diam <- slopes_diam %>%
  filter(dens_model == 3, prod_model == 1, fg %in% c('fg1','fg2','fg3','fg4'), variable %in% c('density','total_production')) %>%
  select(fg, dbh, variable, q50) %>%
  group_by(dbh, variable) %>%
  spread(fg, q50) %>%
  mutate(pioneertobreeder = fg2 - fg4, fasttoslow = fg1 - fg3)

ratio_slopes_diam %>%
  select(dbh, variable, pioneertobreeder, fasttoslow) %>%
  gather(fgs, ratio, -dbh, -variable) %>%
  ggplot(aes(x = dbh, y = ratio, color = fgs)) +
  geom_line() +
  facet_wrap(~ variable) +
  scale_x_log10() +
  theme_bw()

ratio_slopes_light <- slopes_light %>%
  filter(dens_model == '3', prod_model == 1, fg %in% c('fg1','fg2','fg3','fg4'), variable %in% c('density','total_production')) %>%
  select(fg, lightperarea, variable, q50) %>%
  group_by(lightperarea) %>%
  spread(fg, q50) %>%
  mutate(pioneertobreeder = fg2 - fg4, fasttoslow = fg1 - fg3)

ratio_slopes_light %>%
  select(lightperarea, variable, pioneertobreeder, fasttoslow) %>%
  gather(fgs, ratio, -lightperarea, -variable) %>%
  ggplot(aes(x = lightperarea, y = ratio, color = fgs)) +
    geom_line() +
    facet_wrap(~ variable) +
    scale_x_log10() +
    theme_bw()


