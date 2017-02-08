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
  mutate(agb_bin = cut(agb, breaks = n_bins),
         dbh = dbh/10)  # mm to cm

levels(bcicensusdat$agb_bin) <- with(bcicensusdat, binlab2n(agb_bin, islog=F))

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

# The Rouvinen competition index, validated by Contreras, is determined to be pretty good for trees. We will use a 20m radius for the trees.
# Might also need to account for edge effects.

# First calculate how much of a 20m radius circle centered on each tree is contained within the plot.

# function for area of a circular segment with given radius and height of circle

aseg <- function(r, h) r^2 * acos((r-h)/r) - (r-h) * sqrt(2*r*h - h^2)

# Use Monte Carlo method because I suck at geometry.

r <- 20

bcicensusdat <- bcicensusdat %>%
  mutate(overmax_x = gx + 20 > 1000,
         undermin_x = gx - 20 < 0,
         overmax_y = gy + 20 > 500,
         undermin_y = gy - 20 < 0,
         out_x = overmax_x | undermin_x,
         out_y = overmax_y | undermin_y)

bcicensusdat$circle_area <- 1

pb <- txtProgressBar(0, nrow(bcicensusdat), style=3)

for (i in 1:nrow(bcicensusdat)) {
  if (bcicensusdat$out_x[i] | bcicensusdat$out_y[i]) {
    # generate lattice points centered around the point.
    x <- bcicensusdat$gx[i]
    y <- bcicensusdat$gy[i]
    xcoords <- seq(x-r, x+r, by=1)
    ycoords <- seq(y-r, y+r, by=1)
    mat <- expand.grid(x=xcoords,y=ycoords)
    incirc <- sqrt((mat$x - x)^2 + (mat$y - y)^2) <= r
    inrect <- (mat$x >= 0 & mat$x <= 1000) & (mat$y >= 0 & mat$y <= 500)
    bcicensusdat$circle_area[i] <- sum(inrect & incirc)/sum(incirc)
  }
  setTxtProgressBar(pb, i)
}

close(pb)

# Calculate competition index for all the trees
comp_idx <- rep(0, nrow(bcicensusdat))
pb <- txtProgressBar(0, nrow(bcicensusdat), style=3)

for (i in 1:nrow(bcicensusdat)) {
  # Identify neighbors within the 20m radius.
  neighbors <- bcicensusdat[-i, c('gx','gy','dbh')]
  distances <- ((neighbors$gx - bcicensusdat$gx[i])^2 + (neighbors$gy - bcicensusdat$gy[i])^2)^(1/2)
  neighbors <- neighbors[distances <= r, ]
  neighbors$d <- distances[distances <= r]
  neighbors <- neighbors[complete.cases(neighbors), ]
  
  comp_idx[i] <- sum( (neighbors$dbh / bcicensusdat$dbh[i]) * atan(neighbors$dbh / neighbors$d) ) # Rouvinen index.
  
  setTxtProgressBar(pb, i)
}

close(pb)

# Correct competition index for circle overlap.
bcicensusdat$raw_comp_idx <- comp_idx
bcicensusdat <- subset(bcicensusdat, circle_area >= 0.5) # This gets rid of a very few trees at the corners that we don't know a lot about their neighbors.
bcicensusdat$comp_idx <- bcicensusdat$raw_comp_idx / bcicensusdat$circle_area

ggplot(bcicensusdat, aes(x = agb_bin, y = production67)) +
  stat_summary(geom = 'bar', fun.y = 'sum') 