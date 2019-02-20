# Functions to plot the ratio slopes themselves, and the fitted ratio trends over the observed ratio bin points.


gdrive_path <- '~/google_drive'
github_path <- '~/Documents/GitHub'

# Plot ratio slopes -------------------------------------------------------

slopes_light <- read.csv(file.path(gdrive_path, 'ForestLight/data/data_piecewisefits/lightpiecewise/lightpiecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
slopes_diam <- read.csv(file.path(gdrive_path, 'ForestLight/data/data_piecewisefits/piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)

library(tidyverse)

ratio_slopes_diam <- slopes_diam %>%
  filter(dens_model == 3, prod_model == 2, fg %in% c('fg1','fg2','fg3','fg4'), variable %in% c('density','total_production')) %>%
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
  filter(dens_model == 'ln', prod_model == 2, fg %in% c('fg1','fg2','fg3','fg4'), variable %in% c('density','total_production')) %>%
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


# Plot fitted ratios ------------------------------------------------------

# Note: these quantiles are not quite correct and will have to be manually recalculated.

fitted_diam <- read.csv(file.path(gdrive_path, 'ForestLight/data/data_piecewisefits/piecewise_ci_by_fg.csv'), stringsAsFactors = FALSE)
fitted_light <- read.csv(file.path(gdrive_path, 'ForestLight/data/data_piecewisefits/lightpiecewise/lightpiecewise_ci_by_fg.csv'), stringsAsFactors = FALSE)
observed_ratio_diam <- read.csv(file.path(gdrive_path, 'ForestLight/data/data_piecewisefits/observed_ratio_diam.csv'), stringsAsFactors = FALSE)
observed_ratio_light <- read.csv(file.path(gdrive_path, 'ForestLight/data/data_piecewisefits/lightpiecewise/observed_ratio_light.csv'), stringsAsFactors = FALSE)

# The below 2 data frames are from the "total workflow" (Do not run)
# observed_ratio_light <- cbind(light_per_area_bins_allclassified[,c(1,4,5)], fastslow_stats_bylight %>% filter(year==1995) %>% select(-bin, -year), breeder_stats_bylight %>% filter(year == 1995) %>% select(-bin, -year)) %>%
#   mutate(breeder_production_ratio = if_else(breeder_production_ratio == 0, as.numeric(NA), 1/breeder_production_ratio),
#          breeder_density_ratio = if_else(breeder_density_ratio == 0, as.numeric(NA), 1/breeder_density_ratio)) %>%
#   cbind(data.frame(n_fast = densitybin_byyear_bylight[[1]][[2]]$bin_count,
#                    n_slow = densitybin_byyear_bylight[[3]][[2]]$bin_count,
#                    n_pioneer = densitybin_byyear_bylight[[2]][[2]]$bin_count,
#                    n_breeder = densitybin_byyear_bylight[[4]][[2]]$bin_count))
# observed_ratio_diam <- cbind(dbhbin_allclassified[,c(1,4,5)], fastslow_stats_bydiam %>% filter(year==1995) %>% select(-bin, -year), breeder_stats_bydiam %>% filter(year == 1995) %>% select(-bin, -year)) %>%
#   mutate(breeder_production_ratio = if_else(breeder_production_ratio == 0, as.numeric(NA), 1/breeder_production_ratio),
#          breeder_density_ratio = if_else(breeder_density_ratio == 0, as.numeric(NA), 1/breeder_density_ratio)) %>%
#   cbind(data.frame(n_fast = densitybin_byyear_bydiam[[1]][[2]]$bin_count,
#                    n_slow = densitybin_byyear_bydiam[[3]][[2]]$bin_count,
#                    n_pioneer = densitybin_byyear_bydiam[[2]][[2]]$bin_count,
#                    n_breeder = densitybin_byyear_bydiam[[4]][[2]]$bin_count))
# write.csv(observed_ratio_light, file = '~/google_drive/ForestLight/data/data_piecewisefits/lightpiecewise/observed_ratio_light.csv', row.names = FALSE)
# write.csv(observed_ratio_diam, file = '~/google_drive/ForestLight/data/data_piecewisefits/observed_ratio_diam.csv', row.names = FALSE)

fitted_ratio_diam <- fitted_diam %>%
  filter(fg %in% c('fg1','fg2','fg3','fg4'), variable %in% c('density','total_production')) %>%
  select(-year) %>%
  gather(quantile, value, -dens_model, -prod_model, -fg, -dbh, -variable) %>%
  group_by(dens_model, prod_model, dbh, variable, quantile) %>%
  spread(fg, value) %>%
  mutate(pioneertobreeder = fg2/fg4, fasttoslow = fg1/fg3)

quant_colors <- scale_color_manual(values = c(q025 = 'gray40', q50 = 'red', q975 = 'gray40'))  
model_labeller <- labeller(dens_model = c('1' = '1 piece density', '2' = '2 piece density', '3' = '3 piece density', 'w' = 'Weibull density', 'ln' = 'lognormal density'),
                           prod_model = c('1' = '1 piece production', '2' = '2 piece production'))
exl <- expression(paste('Light received per crown area (W m'^-2, ')', sep = ''))
exd <- 'Diameter (cm)'
exf <- 'Fast:Slow'
exb <- 'Pioneer:Breeder'

pdiamdensfast <- ggplot(fitted_ratio_diam %>% filter(quantile %in% c('q025', 'q50', 'q975'), variable %in% 'density'), aes(x = dbh, y = fasttoslow)) +
  geom_line(aes(group = quantile, color = quantile)) +
  geom_point(data = observed_ratio_diam, aes(x = bin_midpoint, y = fastslow_density_ratio)) +
  scale_x_log10(name = exd) + 
  scale_y_log10(name = exf) +
  facet_grid(dens_model ~ prod_model, scales = 'free_y', labeller = model_labeller) +
  theme_bw() +
  quant_colors +
  ggtitle('Fast:Slow Density ratio', 'scaled by diameter')
pdiamdenspio <- ggplot(fitted_ratio_diam %>% filter(quantile %in% c('q025', 'q50', 'q975'), variable %in% 'density'), aes(x = dbh, y = pioneertobreeder)) +
  geom_line(aes(group = quantile, color = quantile)) +
  geom_point(data = observed_ratio_diam, aes(x = bin_midpoint, y = breeder_density_ratio)) +
  scale_x_log10(name = exd) + 
  scale_y_log10(name = exb) +
  facet_grid(dens_model ~ prod_model, scales = 'free_y', labeller = model_labeller) +
  theme_bw() +
  quant_colors +
  ggtitle('Pioneer:Breeder Density ratio', 'scaled by diameter')
pdiamprodfast <- ggplot(fitted_ratio_diam %>% filter(quantile %in% c('q025', 'q50', 'q975'), variable %in% 'total_production'), aes(x = dbh, y = fasttoslow)) +
  geom_line(aes(group = quantile, color = quantile)) +
  geom_point(data = observed_ratio_diam, aes(x = bin_midpoint, y = fastslow_production_ratio)) +
  scale_x_log10(name = exd) + 
  scale_y_log10(name = exf) +
  facet_grid(dens_model ~ prod_model, scales = 'free_y', labeller = model_labeller) +
  theme_bw() +
  quant_colors +
  ggtitle('Fast:Slow Production ratio', 'scaled by diameter')
pdiamprodpio <- ggplot(fitted_ratio_diam %>% filter(quantile %in% c('q025', 'q50', 'q975'), variable %in% 'total_production'), aes(x = dbh, y = pioneertobreeder)) +
  geom_line(aes(group = quantile, color = quantile)) +
  geom_point(data = observed_ratio_diam, aes(x = bin_midpoint, y = breeder_production_ratio)) +
  scale_x_log10(name = exd) + 
  scale_y_log10(name = exb) +
  facet_grid(dens_model ~ prod_model, scales = 'free_y', labeller = model_labeller) +
  theme_bw() +
  quant_colors +
  ggtitle('Pioneer:Breeder Production ratio', 'scaled by diameter')


fitted_ratio_light <- fitted_light %>%
  filter(fg %in% c('fg1','fg2','fg3','fg4'), variable %in% c('density','total_production')) %>%
  select(-year) %>%
  gather(quantile, value, -dens_model, -prod_model, -fg, -lightperarea, -variable) %>%
  group_by(dens_model, prod_model, lightperarea, variable, quantile) %>%
  spread(fg, value) %>%
  mutate(pioneertobreeder = fg2/fg4, fasttoslow = fg1/fg3)


plightdensfast <- ggplot(fitted_ratio_light %>% filter(quantile %in% c('q025', 'q50', 'q975'), variable %in% 'density'), aes(x = lightperarea, y = fasttoslow)) +
  geom_line(aes(group = quantile, color = quantile)) +
  geom_point(data = observed_ratio_light, aes(x = bin_midpoint, y = fastslow_density_ratio)) +
  scale_x_log10(name = exl) + 
  scale_y_log10(name = exf) +
  facet_grid(dens_model ~ prod_model, scales = 'free_y', labeller = model_labeller) +
  theme_bw() +
  quant_colors +
  ggtitle('Fast:Slow Density ratio', 'scaled by light per area')
plightdenspio <- ggplot(fitted_ratio_light %>% filter(quantile %in% c('q025', 'q50', 'q975'), variable %in% 'density'), aes(x = lightperarea, y = pioneertobreeder)) +
  geom_line(aes(group = quantile, color = quantile)) +
  geom_point(data = observed_ratio_light, aes(x = bin_midpoint, y = breeder_density_ratio)) +
  scale_x_log10(name = exl) + 
  scale_y_log10(name = exb) +
  facet_grid(dens_model ~ prod_model, scales = 'free_y', labeller = model_labeller) +
  theme_bw() +
  quant_colors +
  ggtitle('Pioneer:Breeder Density ratio', 'scaled by light per area')
plightprodfast <- ggplot(fitted_ratio_light %>% filter(quantile %in% c('q025', 'q50', 'q975'), variable %in% 'total_production'), aes(x = lightperarea, y = fasttoslow)) +
  geom_line(aes(group = quantile, color = quantile)) +
  geom_point(data = observed_ratio_light, aes(x = bin_midpoint, y = fastslow_production_ratio)) +
  scale_x_log10(name = exl) + 
  scale_y_log10(name = exf) +
  facet_grid(dens_model ~ prod_model, scales = 'free_y', labeller = model_labeller) +
  theme_bw() +
  quant_colors +
  ggtitle('Fast:Slow Production ratio', 'scaled by light per area')
plightprodpio <- ggplot(fitted_ratio_light %>% filter(quantile %in% c('q025', 'q50', 'q975'), variable %in% 'total_production'), aes(x = lightperarea, y = pioneertobreeder)) +
  geom_line(aes(group = quantile, color = quantile)) +
  geom_point(data = observed_ratio_light, aes(x = bin_midpoint, y = breeder_production_ratio)) +
  scale_x_log10(name = exl) + 
  scale_y_log10(name = exb) +
  facet_grid(dens_model ~ prod_model, scales = 'free_y', labeller = model_labeller) +
  theme_bw() +
  quant_colors +
  ggtitle('Pioneer:Breeder Production ratio', 'scaled by light per area')

fpfig <- '~/google_drive/ForestLight/figs/lightpowerlaws_feb2019/ratio_fitted_observed'

pdf(file.path(fpfig, 'ratios_scaled_by_diameter.pdf'), height = 7, width = 9)
  pdiamdensfast
  pdiamprodfast
  pdiamdenspio
  pdiamprodpio
dev.off()

pdf(file.path(fpfig, 'ratios_scaled_by_lightperarea.pdf'), height = 11, width = 9)
  plightdensfast
  plightprodfast
  plightdenspio
  plightprodpio
dev.off()