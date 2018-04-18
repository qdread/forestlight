# Load bci full census data.

fp <- 'C:/Users/Q/Dropbox/projects/forestlight/bcidata/'
fp <- '~/data/bcidata'
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

# levels(bcicensusdat$agb_bin) <- with(bcicensusdat, binlab2n(agb_bin, islog=F))

# Calculate minimum distance to a larger tree for all the trees
# distance_to_bigger <- rep(0, nrow(bcicensusdat))
# pb <- txtProgressBar(0, nrow(bcicensusdat), style=3)
# 
# for (i in 1:nrow(bcicensusdat)) {
#   neighbors <- bcicensusdat[-i, c('gx','gy','dbh')]
#   distances <- ((neighbors$gx - bcicensusdat$gx[i])^2 + (neighbors$gy - bcicensusdat$gy[i])^2)^(1/2)
#   is_bigger <- neighbors$dbh > bcicensusdat$dbh[i]
#   bigger_dists <- distances[is_bigger]
#   if (length(bigger_dists) > 0) distance_to_bigger[i] <- min(bigger_dists)
#   setTxtProgressBar(pb, i)
# }
# 
# close(pb)

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
  
  comp_idx[i] <- sum( (neighbors$dbh / bcicensusdat$dbh[i]) * atan(neighbors$dbh / neighbors$d), na.rm=TRUE) # Rouvinen index.
  
  setTxtProgressBar(pb, i)
}

close(pb)

# Correct competition index for circle overlap.
bcicensusdat$raw_comp_idx <- comp_idx
bcicensusdat <- subset(bcicensusdat, circle_area >= 0.5) # This gets rid of a very few trees at the corners that we don't know a lot about their neighbors.
bcicensusdat$comp_idx <- bcicensusdat$raw_comp_idx / bcicensusdat$circle_area

write.csv(bcicensusdat, file = file.path(fp, 'bci_compidx.csv'), row.names = FALSE)

# ggplot(bcicensusdat, aes(x = agb_bin, y = production67)) +
#   stat_summary(geom = 'bar', fun.y = 'sum') 

# Load bci competition index data calculated on the cluster.

fp <- 'C:/Users/Q/Dropbox/projects/forestlight/bcidata/'
bcicensusdat <- read.csv(file.path(fp, 'bci_compidx.csv'), stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)

# Plot competition index by bin

bcicensusdat$bin10 <- cut(bcicensusdat$agb, breaks = 10)

ggplot(bcicensusdat, aes(x = bin10, y = comp_idx)) +
  stat_summary(geom = 'bar', fun.y = 'mean') +
  stat_summary(geom = 'errorbar', width = 0) +
  theme_bw() + theme(panel.grid = element_blank())

pciagb <- ggplot(bcicensusdat, aes(x = agb, y = comp_idx)) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = 'BCI 50-hectare plot', subtitle = 'Rouvinen competition index versus aboveground biomass (log-log)')
ggsave('C:/Users/Q/Google Drive/ForestLight/figs/bci50hectare_comp_vs_biomass.png', height=6, width=6)

# Run the biomass binning algorithm on this dataset.


biomass_diameter_bin <- function(x, bin_width = 1, n_per_bin = 10, grp = 'biomass') {
  
  # Set the grouping variable.
  if (grp == 'biomass') x <- transform(x, grpvar = agb) else x <- transform(x, grpvar = dbh)
  
  # First, create small bins. Default bin width is 1 kg but it will work for others.
  grpvarrange <- range(x$grpvar)
  round_range <- plyr::round_any(grpvarrange, bin_width, floor)
  bin_min <- seq(round_range[1], round_range[2], by = bin_width) # Corrected to deal with diameter cutoff at 1.
  bin_max <- bin_min + bin_width
  bins_small <- rep(0, nrow(x))
  for (i in 1:nrow(x)) {
    # Note: changed this to have the bottom open, not top
    bins_small[i] <- bin_width * which(x$grpvar[i] >= bin_min & x$grpvar[i] < bin_max)
  }
  
  x <- x %>% mutate(bins_small = bins_small) %>% arrange(grpvar)
  
  
  
  # Next, group the bins until at least n_per_bin individuals are reached in each bin. 
  # *n_per_bin set to 10 by default*
  
  bin_n_indiv <- rep(1, nrow(x))
  bin_id <- 1
  bin_indiv <- 0
  cumulative_total_indiv <- 1
  max_of_bin <- bin_width
  
  for (i in 1:nrow(x)) {
    bin_indiv <- bin_indiv + 1 # Increment number of individuals in bin by 1
    # If tree i has crossed into the next bin, 
    # check whether n_per_bin individuals have been reached. If so, make a new bin.
    if (x$bins_small[i] > max_of_bin & bin_indiv >= n_per_bin) {
      bin_n_indiv[(cumulative_total_indiv -  1):(cumulative_total_indiv - 1 + bin_indiv)] <- bin_id
      cumulative_total_indiv <- cumulative_total_indiv + bin_indiv
      bin_id <- bin_id + 1
      bin_indiv <- 0
      max_of_bin <- plyr::round_any(x$bins_small[i], bin_width, ceiling)
    }
    # If we've reached the final row without exactly reaching a multiple of n_per_bin, 
    # put the remaining largest trees into their own bin.
    if (i == nrow(x) & bin_indiv > 0) {
      bin_n_indiv[(cumulative_total_indiv - 1):nrow(x)] <- bin_id
    }
  }
  
  # For each of these bins containing a minimum of n individuals, do the following steps:
  
  x %>% 
    mutate(bin_n_indiv = bin_n_indiv) %>% 
    group_by(bin_n_indiv) %>% 
    summarize(production_sum = sum(production67),                  # Gets raw biomass production sum for the bin.
              median_grp = median(grpvar),                     # Finds median biomass or diameter of bin.
              mid_grp = (plyr::round_any(min(grpvar), bin_width, floor) + plyr::round_any(max(grpvar), bin_width, ceiling))/2, # Find biomass or diameter midpoint of the bin.
              n_smallbins = 1 + max(bins_small) - min(bins_small),       # Determines how many virtual bins are contained within the bin.
              production_perbin = production_sum / n_smallbins,  # Divides raw sum by number of bins.
              CI_mean = mean(comp_idx),                            # Gets mean comp index.
              CI_se = sd(comp_idx)/sqrt(length(comp_idx)),              # Gets standard error of comp index.
              n_individuals = n())                             # Returns number of individual trees in the bin.
  
  
}

# Transform biomass into kg.
bcicensusdat <- mutate(bcicensusdat, agb = agb*1000, production67 = production67*1000)

bciall_binned <- biomass_diameter_bin(bcicensusdat, bin_width = 1, n_per_bin = 10, grp = 'biomass')

sc_ylog <- scale_y_log10(name = 'Biomass production (kg/y)', breaks = c(0.1, 1, 10), limits = c(0.1, 55))
sc_xbiomass <- scale_x_log10(name = 'Biomass (kg)', breaks = c(0.1,1,10,100,1000))
sc_xdiam <- scale_x_continuous(name = 'Diameter (cm)')
th_scatter <- theme_bw() + theme(panel.grid = element_blank())


ggplot(massbci, aes(x=mid_grp, y=production_perbin)) +
  geom_point() +
  sc_ylog + sc_xbiomass + th_scatter +
  ggtitle('BCI, binned by biomass, log y axis')