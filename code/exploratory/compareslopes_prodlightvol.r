# Look at the total and individual light data compared to production data

library(tidyverse)
library(forestscaling)

load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')

# Plot indiv prod, indiv light and indiv crown volume

p1 <- ggplot(alltreedat[[3]], aes(x = dbh_corr, y = production)) +
  geom_hex() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 2, intercept = 0, color = 'red')
p2 <- ggplot(alltreedat[[3]], aes(x = dbh_corr, y = light_received)) +
  geom_hex() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 2, intercept = 0, color = 'red')
p3 <- ggplot(alltreedat[[3]], aes(x = dbh_corr, y = crownvolume)) +
  geom_hex() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 2, intercept = 0, color = 'red')

gridExtra::grid.arrange(p1, p2, p3, nrow = 1)

# Create log bins for total production, incoming light, and crown volume.

with(alltreedat[[3]], logbin(x = dbh_corr, y = ))