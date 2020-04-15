library(tidyverse)
library(brms)

lai_depth <- read_csv('data/LAI_Depth_Kitajima2005.csv') %>%
  mutate(log_LAI = log10(LAI),
         log_depth = log10(Depth))

# Fit Bayesian random-intercept, random-slope model
lai_fit <- brm(log_LAI ~ log_depth + (log_depth | Species), data = lai_depth, family = 'gaussian', chains = 3, iter = 10000, warmup = 9000)


summary(lai_fit)

# Fit Bayesian random-intercept, random-slope model, no Cecropia
lai_fit2 <- brm(log_LAI ~ log_depth + (log_depth | Species), data = lai_depth %>% filter(!Species %in% 'Cecropia'), family = 'gaussian', chains = 3, iter = 10000, warmup = 9000)


summary(lai_fit2)

curve(1 - exp(-0.5 * x), from = 0.1, to = 15, log = 'x', ylim = c(0,1))
curve(1 - exp(-0.5 * (1.58 * x ^ 0.35)), from = 0.1, to = 15, log = 'x', ylim = c(0,1), col = 'red', add  =TRUE)
curve(1 - exp(-0.5 * (2.29 * x ^ 0.37)), from = 0.1, to = 15, log = 'x', col = 'green', add = TRUE)


