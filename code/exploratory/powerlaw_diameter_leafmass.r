library(actuar)

n <- 100000
alpha <- 1 # Slope of -2

x <- rpareto(n, alpha, 1)

library(ggplot2)
library(dplyr)

ggplot(data.frame(x=x), aes(x)) + geom_histogram() + scale_x_log10() + scale_y_log10()

library(poweRlaw)
m <- conpl$new()
m$setXmin(1)
m$setPars(2)

diams <- dist_rand(m, n)
leafmasses <- diams^(3/2)
productions <- diams^2

# Create log bin.

source('~/Documents/GitHub/forestlight/code/allfunctions27july.r')

diam_bin <- logbin(diams, n = 20)
leafmass_bin <- logbin(leafmasses, n = 20)
totalprod_bin_bydiameter <- logbin(diams, productions, n = 20)
totalprod_bin_byleafmass <- logbin(leafmasses, productions, n = 20)


ggplot(diam_bin %>% filter(bin_count>0), aes(x = bin_midpoint, y = bin_value)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope=-2, intercept=5)

ggplot(leafmass_bin %>% filter(bin_count>0), aes(x = bin_midpoint, y = bin_value)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = -3/2, intercept = 5)

ggplot(totalprod_bin_bydiameter %>% filter(bin_count>0), aes(x = bin_midpoint, y = bin_value)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() 

ggplot(totalprod_bin_byleafmass %>% filter(bin_count>0), aes(x = bin_midpoint, y = bin_value)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = -1/4, intercept = 5)

ggplot(data.frame(x=leafmasses,y=productions), aes(x,y)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 3, intercept = 1)


# Test this out with data -------------------------------------------------

load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')

dat <- alltreedat[[3]] %>% 
  select(dbh_corr, agb_corr, production) %>%
  mutate(leafmass = dbh_corr ^ (3/2))

dat_diam_bin <- logbin(dat$dbh_corr, n = 20)
dat_leafmass_bin <- logbin(dat$leafmass, n = 20)
dat_totalprod_bin_bydiameter <- logbin(dat$dbh_corr, dat$production, n = 20)
dat_totalprod_bin_byleafmass <- logbin(dat$leafmass, dat$production, n = 20)


ggplot(dat_diam_bin %>% filter(bin_count>0), aes(x = bin_midpoint, y = bin_value)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope=-2, intercept=5)

ggplot(dat_leafmass_bin %>% filter(bin_count>0), aes(x = bin_midpoint, y = bin_value)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = -(3/2), intercept = 5.5)

ggplot(dat_totalprod_bin_bydiameter %>% filter(bin_count>0), aes(x = bin_midpoint, y = bin_value)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 0.25, intercept = 3)

ggplot(dat_totalprod_bin_byleafmass %>% filter(bin_count>0), aes(x = bin_midpoint, y = bin_value)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() 
