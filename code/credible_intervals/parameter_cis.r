# Plot the confidence intervals of parameters

param_ci <- read.csv('C:/Users/Q/google_drive/ForestLight/data/paramci_by_fg.csv', stringsAsFactors = FALSE)

param_ci[param_ci$parameter=='alpha', c(6,9:13)] <- param_ci[param_ci$parameter=='alpha', c(6,9:13)] + 1 # correct alpha

library(dplyr)
library(ggplot2)

param_names <- c(
  alpha = 'density Pareto slope',
  beta0 = 'production power law intercept',
  beta1 = 'production power law slope',
  a = 'production exponential intercept',
  b = 'production exponential slope',
  c = 'production exponential, extra constant',
  m = 'density Weibull shape',
  n = 'density Weibull scale'
)

param_ci <- mutate(param_ci, parameter = factor(parameter, levels = c('alpha', 'm', 'n', 'beta0', 'beta1', 'a', 'b', 'c')))

ggplot(param_ci %>% filter(year == 1995, dens_model == 'pareto', prod_model == 'powerlaw'), 
       aes(x = fg, y = q50, ymin = q025, ymax = q975)) +
  geom_pointrange() +
  facet_wrap(~ parameter, scales = 'free', labeller = labeller(parameter=param_names)) +
  theme_classic() +
  theme(strip.background = element_blank(), panel.border = element_rect(color = 'black', fill = 'transparent')) +
  ggtitle('Pareto density, power law production', '1995')

ggplot(param_ci %>% filter(year == 1995, dens_model == 'weibull', prod_model == 'powerlaw'), 
       aes(x = fg, y = q50, ymin = q025, ymax = q975)) +
  geom_pointrange() +
  facet_wrap(~ parameter, scales = 'free', labeller = labeller(parameter=param_names)) +
  theme_classic() +
  theme(strip.background = element_blank(), panel.border = element_rect(color = 'black', fill = 'transparent')) +
  ggtitle('Weibull density, power law production', '1995')

ggplot(param_ci %>% filter(year == 1995, dens_model == 'pareto', prod_model == 'powerlawexp'), 
       aes(x = fg, y = q50, ymin = q025, ymax = q975)) +
  geom_pointrange() +
  facet_wrap(~ parameter, scales = 'free', labeller = labeller(parameter=param_names)) +
  theme_classic() +
  theme(strip.background = element_blank(), panel.border = element_rect(color = 'black', fill = 'transparent')) +
  ggtitle('Pareto density, power law times exponential production', '1995')

ggplot(param_ci %>% filter(year == 1995, dens_model == 'weibull', prod_model == 'powerlawexp'), 
       aes(x = fg, y = q50, ymin = q025, ymax = q975)) +
  geom_pointrange() +
  facet_wrap(~ parameter, scales = 'free', labeller = labeller(parameter=param_names)) +
  theme_classic() +
  theme(strip.background = element_blank(), panel.border = element_rect(color = 'black', fill = 'transparent')) +
  ggtitle('Weibull density, power law times exponential production', '1995')
