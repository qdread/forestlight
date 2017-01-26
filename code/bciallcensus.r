# Load bci full census data.

fp <- 'C:/Users/Q/Dropbox/projects/forestlight/bcidata/'
load(file.path(fp, 'bci.full6.rdata'))
load(file.path(fp, 'bci.full7.rdata'))

# Get agb growth rate from 6 to 7.
# measured in tonnes
# 2005-2010 (5 years)

bci.full7$production67 <- pmax((bci.full7$agb - bci.full6$agb)/5, 0, na.rm = T)

n_bins <- 100

library(dplyr)
library(ggplot2)

bcicensusdat <- bci.full7 %>%
  filter(DFstatus == 'alive') %>%
  mutate(agb_bin = cut(agb, breaks = n_bins))

levels(bcidat$agb_bin) <- with(bcidat, binlab2n(agb_bin, islog=F))

# Calculate minimum distance to a larger tree for all the trees
distance_to_bigger <- rep(0, nrow(bcicensusdat))
pb <- txtProgressBar(0, nrow(bcicensusdat), style=3)

for (i in 1:nrow(bcicensusdat)) {
  neighbors <- bcicensusdat[-i, c('gx','gy','dbh')]
  distances <- ((neighbors$gx - bcicensusdat$gx[i])^2 + (neighbors$gy - bcicensusdat$gy[i])^2)^(1/2)
  is_bigger <- neighbors$dbh > bcicensusdat$dbh[i]
  bigger_dists <- distances[is_bigger]
  if (length(bigger_dists) > 0) distance_to_bigger[i] <- min(bigger_dists)
  setTxtProgressBar(pb, i)
}

close(pb)

ggplot(bcicensusdat, aes(x = agb_bin, y = production67)) +
  stat_summary(geom = 'bar', fun.y = 'sum') 