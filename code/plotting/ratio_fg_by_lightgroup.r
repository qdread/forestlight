# Plots of density, growth, total growth, and ratios of different functional groups, split up by light quantile
# QDR / Forestlight / 10 June 2019

# Modified 11 June 2019: directly regress pca axis 1 on light and size. (see bottom of script)

# Load all the raw data.
library(tidyverse)

#### change these file paths if needed ####
#### John please do not use the setwd() command in shared scripts!!! ####
#### Just change the paths to the location of google drive and the forestlight github repository on your machine. ####
#### thanks!!! ####
gdrive_path <- '~/google_drive'
github_path <- '~/Documents/GitHub/forestlight'

load(file.path(gdrive_path, 'ForestLight/data/rawdataobj_alternativecluster.r'))
source(file.path(github_path, 'code/allfunctions27july.r'))

# Calculate the DBH bin bounds (same as the ones used for all other dbh analysis)
n <- 20
dbh_range <- log(range(unlist(map(alltreedat, 'dbh_corr'))))
log_bin_edges <- seq(dbh_range[1], dbh_range[2], length.out = n + 1)
bin_edges <- exp(log_bin_edges) # get edges of bins
bin_midpoints <- exp(log_bin_edges[-length(log_bin_edges)] + diff(log_bin_edges)/2)

# Split up into quantiles by light per unit crown area.

# This is the function to determine max incoming light per area at BCI.
insolation <- function(lat) {
  lat <- lat * pi/180 # to radians
  y <- sin(lat)
  0.25 * 1367 * (1 - 0.482 * (3*y^2 - 1)/2)
}

# Insolation at BCI, 9.2 degrees N
(insol_bci <- insolation(9.2))

# These are the cutoffs of the 3 light environments.
# The light values are all proportions. Split evenly on log scale
cutoffs <- c(0, exp(log(insol_bci) * c(1/3, 2/3, 1)))


dat95 <- alltree_light_95 %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)),
         light = light * insol_bci) %>%
  mutate(light_group = cut(light, breaks = cutoffs)) %>%
  mutate(light_group = factor(light_group, labels = c('low', 'medium', 'high'))) %>%
  select(sp, light_group, fg, X1, X2, dbh_corr, production, light, crownarea, crownvolume)

# Create ratios within each of the bins by counting individuals per bin.
dat95_binned <- dat95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = bin_edges, include.lowest = TRUE)) %>%
  group_by(light_group, dbh_bin, fg) %>%
  summarize(abundance = n(),
            production = sum(production))

dat95_abundance_ratios <- dat95_binned %>%
  select(-production) %>%
  group_by(light_group) %>%
  spread(fg, abundance, fill = 0)  %>%
  mutate(fast_slow_ratio = fg1/fg3,
         pioneer_breeder_ratio = fg2/fg4)

dat95_production_ratios <- dat95_binned %>%
  select(-abundance) %>%
  group_by(light_group) %>%
  spread(fg, production, fill = 0)  %>%
  mutate(fast_slow_ratio = fg1/fg3,
         pioneer_breeder_ratio = fg2/fg4)

ratiosbylight1995 <- rbind(data.frame(variable = 'density', dat95_abundance_ratios),
                           data.frame(variable = 'total growth', dat95_production_ratios)) %>%
  mutate(bin_midpoint = bin_midpoints[as.numeric(dbh_bin)])

### Make plots.

ggplot(ratiosbylight1995 %>% filter(is.finite(fast_slow_ratio), fast_slow_ratio > 0), aes(x = bin_midpoint, y = fast_slow_ratio)) +
  facet_grid(variable ~ light_group) +
  geom_point() +
  scale_x_log10(name = 'Diameter (cm)') + scale_y_log10('Fast to slow ratio') +
  theme_bw()
ggsave(file.path(gdrive_path, 'ForestLight/figs/ratio_bylightenvironment_fasttoslow.pdf'), height = 6, width = 9)

ggplot(ratiosbylight1995 %>% filter(is.finite(pioneer_breeder_ratio), pioneer_breeder_ratio > 0), aes(x = bin_midpoint, y = pioneer_breeder_ratio)) +
  facet_grid(variable ~ light_group) +
  geom_point() +
  scale_x_log10(name = 'Diameter (cm)') + scale_y_log10('Pioneer to breeder ratio') +
  theme_bw()         
ggsave(file.path(gdrive_path, 'ForestLight/figs/ratio_bylightenvironment_pioneertobreeder.pdf'), height = 6, width = 9)


# Regressions with axis X1 ------------------------------------------------

# Slow-fast axis by log dbh
# Note: this is probably not the greatest statistically speaking because they treat all individuals' PCA scores as independent

fast_by_size <- lm(X1 ~ log10(dbh_corr), data = dat95)
summary(fast_by_size)
confint(fast_by_size) # Coefficient on size is between 1.04 and 1.09

fast_by_size_light <- lm(X1 ~ log10(dbh_corr) + log10(light), data = dat95)
summary(fast_by_size_light)
confint(fast_by_size_light) # Coefficient on size increases because coefficient on light is actually negative

# What we see is that the larger the tree, the faster the score tends to be. 
# This is because of the large number of small saplings in the slow grouping.
# Adding light to the regression makes the trend even stronger. The more light, the slower the score tends to be.

library(lavaan)
simplemodel <- '
  X1 ~ light + size
  light ~ size
'

semdat <- dat95 %>%
  transmute(size = log10(dbh_corr), light = log10(light), X1 = X1)

sem1 <- sem(model = simplemodel, data = semdat)
