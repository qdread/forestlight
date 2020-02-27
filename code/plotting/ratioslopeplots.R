library(tidyverse)

ratioslopes <- read_csv('~/google_drive/ForestLight/data/data_piecewisefits/ratio_slope_ci.csv')

ggplot(ratioslopes, aes(x=dbh, y=q50, ymin=q025, ymax=q975, group=ratio, color=ratio, fill = ratio)) +
  facet_wrap(~ variable) +
  geom_ribbon(alpha = 0.5) +
  geom_line(size = 1) +
  scale_x_log10() +
  theme_bw()

min(abs(ratioslopes$dbh - 15))

ratioslopes %>% filter(abs(dbh-15) == min(abs(dbh-15)))


# plot ratios again -------------------------------------------------------


