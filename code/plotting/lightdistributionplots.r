# Make initial plots of power laws

load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')
source('code/allfunctions27july.r')

library(dplyr)
library(ggplot2)

fakebin_across_years <- function(dat_values, dat_classes, edges, mean_type = 'geometric', n_census = 5) {
  qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  # add some padding just in case
  mins <- edges$bin_min
  mins[1] <- 0
  maxes <- edges$bin_max
  maxes[length(maxes)] <- Inf
  
  binstats <- t(sapply(1:length(mins), function(x) {
    indivs <- dat_values[dat_classes >= mins[x] & dat_classes < maxes[x]]
    if (mean_type == 'geometric') {
      mean_n <- exp(mean(log(indivs)))
      sd_n <- exp(sd(log(indivs)))
      ci_width <- qnorm(0.975) * sd(log(indivs)) / sqrt(length(indivs))
      ci_min <- exp(mean(log(indivs)) - ci_width)
      ci_max <- exp(mean(log(indivs)) + ci_width)
      
    } else {
      mean_n <- mean(indivs)
      sd_n <- sd(indivs)
      ci_width <- qnorm(0.975) * sd(indivs) / sqrt(length(indivs))
      ci_min <- mean_n - ci_width
      ci_max <- mean_n + ci_width
    }
    c(mean = mean_n, 
      sd = sd_n,
      quantile(indivs, probs = qprobs),
      ci_min = ci_min,
      ci_max = ci_max)
  }))
  dimnames(binstats)[[2]] <- c('mean', 'sd', 'q025', 'q25', 'median', 'q75', 'q975', 'ci_min', 'ci_max')
  data.frame(bin_midpoint = edges$bin_midpoint,
             bin_min = edges$bin_min,
             bin_max = edges$bin_max,
             mean_n_individuals = edges$bin_count / n_census,
             binstats)
}

# Create bins.
area_core <- 42.84

lightrange <- with(alltree_light_95, range(light_received/crownarea)) # 1.09 to 412

lightbinbounds <- exp(seq(log(lightrange[1]), log(lightrange[2]), length.out = 21))

alltree_light_95 <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), as.character(NA)))

lightabundbins_all <- with(alltree_light_95, logbin(x = light_received/crownarea, n = 20))
lightproductionbins_all <- with(alltree_light_95, logbin_setedges(x = light_received/crownarea, y = production, edges = lightabundbins_all))
lightindivprodbins_all <- with(alltree_light_95, fakebin_across_years(dat_values = production, dat_classes = light_received/crownarea, edges = lightabundbins_all, mean_type = 'geometric', n_census = 1))

lightabundbins_fg <- alltree_light_95 %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$light_received/.$crownarea, edges = lightabundbins_all))

lightproductionbins_fg <- alltree_light_95 %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$light_received/.$crownarea, y = .$production, edges = lightabundbins_all))

lightindivprodbins_fg <- alltree_light_95 %>%
  group_by(fg) %>%
  do(fakebin_across_years(dat_values = .$production, dat_classes = .$light_received/.$crownarea, edges = lightabundbins_all, mean_type = 'geometric', n_census = 1))

lightabundbins_fg <- rbind(data.frame(fg = 'all', lightabundbins_all, stringsAsFactors = FALSE), as.data.frame(lightabundbins_fg)) %>%
  mutate(bin_value = bin_value / area_core)
lightproductionbins_fg <- rbind(data.frame(fg = 'all', lightproductionbins_all, stringsAsFactors = FALSE), as.data.frame(lightproductionbins_fg)) %>%
  mutate(bin_value = bin_value / area_core)
lightindivprodbins_fg <- rbind(data.frame(fg = 'all', lightindivprodbins_all, stringsAsFactors = FALSE), as.data.frame(lightindivprodbins_fg))



p1 <- ggplot(lightabundbins_fg %>% filter(!is.na(fg)), aes(x = bin_midpoint, y = bin_value)) +
  geom_point() +
  facet_wrap(~ fg) +
  scale_x_log10(name = 'Light received per unit crown area') +
  scale_y_log10() +
  theme_bw() +
  ggtitle('Density')
  
p1b <- ggplot(lightabundbins_fg %>% filter(!is.na(fg), !fg %in% 'all'), aes(x = bin_midpoint, y = bin_value, color = fg)) +
  geom_point() +
  scale_x_log10(name = 'Light received per unit crown area') +
  scale_y_log10() +
  theme_bw() +
  ggtitle('Density')

p2 <- ggplot(alltree_light_95 %>% filter(!is.na(fg)), aes(x = light_received/crownarea, y = production)) +
  geom_hex(alpha = 0.4) +
  geom_hex(data = alltree_light_95 %>% mutate(fg = 'all'), alpha = 0.4) +
  geom_pointrange(data = lightindivprodbins_fg %>% filter(!is.na(fg)), aes(x = bin_midpoint, y = median, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = 'Light received per unit crown area') +
  scale_y_log10() +
  theme_bw() +
  ggtitle('Production (individual)') +
  scale_fill_gradient(trans = 'log', low = 'skyblue', high = 'darkblue')

p2b <- ggplot(lightindivprodbins_fg %>% filter(!is.na(fg), !fg %in% 'all'), aes(x = bin_midpoint, y = median, ymin = q25, ymax = q75, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Light received per unit crown area') +
  scale_y_log10() +
  theme_bw() +
  ggtitle('Production (individual)')

p3 <- ggplot(lightproductionbins_fg %>% filter(!is.na(fg)), aes(x = bin_midpoint, y = bin_value)) +
  geom_point() +
  facet_wrap(~ fg) +
  scale_x_log10(name = 'Light received per unit crown area') +
  scale_y_log10() +
  theme_bw() +
  ggtitle('Production (total)')

p3b <- ggplot(lightproductionbins_fg %>% filter(!is.na(fg), !fg %in% 'all'), aes(x = bin_midpoint, y = bin_value, color = fg)) +
  geom_point() +
  scale_x_log10(name = 'Light received per unit crown area') +
  scale_y_log10() +
  theme_bw() +
  ggtitle('Production (total)')

fpfig <- '~/google_drive/ForestLight/figs/lightpowerlaws_feb2019'
ggsave(file.path(fpfig, 'densitybylight_separate.png'), p1, height = 5, width = 9, dpi = 300)
ggsave(file.path(fpfig, 'densitybylight_together.png'), p1b, height = 5, width = 5, dpi = 300)
ggsave(file.path(fpfig, 'productionindividualbylight_separate.png'), p2, height = 5, width = 9, dpi = 300)
ggsave(file.path(fpfig, 'productionindividuallight_together.png'), p2b, height = 5, width = 5, dpi = 300)
ggsave(file.path(fpfig, 'productiontotalbylight_separate.png'), p3, height = 5, width = 9, dpi = 300)
ggsave(file.path(fpfig, 'productiontotalbylight_together.png'), p3b, height = 5, width = 5, dpi = 300)


# Add fitted lines --------------------------------------------------------

light_ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/lightpiecewise/lightpiecewise_ci_by_fg.csv', stringsAsFactors = FALSE) 
area_core <- 42.84

light_ci_df$fg[light_ci_df$fg == 'alltree'] <- 'all'

light_pred_dens <- light_ci_df %>%
  filter(variable == 'density') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

light_fitted_indivprod <- light_ci_df %>%
  filter(variable == 'production_fitted') %>%
  select(-variable)

light_fitted_totalprod <- light_ci_df %>%
  filter(variable == 'total_production_fitted') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

light_pred_indivprod <- light_ci_df %>%
  filter(variable == 'production') %>%
  select(-variable)

light_pred_totalprod <- light_ci_df %>%
  filter(variable == 'total_production') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 


p1fits <- ggplot(lightabundbins_fg %>% filter(!is.na(fg), bin_count > 10)) +
  geom_ribbon(data = light_pred_dens %>% filter(prod_model == 1, !dens_model %in% '1', !fg %in% 'unclassified'), aes(x = lightperarea, ymin = q025, ymax = q975, group = dens_model), fill = 'gray80') +
  geom_line(data = light_pred_dens %>% filter(prod_model == 1, !dens_model %in% '1', !fg %in% 'unclassified'), aes(x = lightperarea, y = q50, group = dens_model, color = dens_model)) +
  geom_point(aes(x = bin_midpoint, y = bin_value)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = 'Light received per unit crown area (Wm-2)') +
  scale_y_log10(name = 'Density (trees ha-1 per watt)', limits = c(1e-5, 1e2)) +
  theme_bw() +
  ggtitle('Density')

p2fits <- ggplot() +
  geom_ribbon(data = light_fitted_indivprod %>% filter(dens_model == '1', !fg %in% 'unclassified'), aes(x = lightperarea, ymin = q025, ymax = q975, group = factor(prod_model)), fill = 'gray80') +
  geom_line(data = light_fitted_indivprod %>% filter(dens_model == '1', !fg %in% 'unclassified'), aes(x = lightperarea, y = q50, group = factor(prod_model), color = factor(prod_model))) +
  geom_pointrange(data = lightindivprodbins_fg %>% filter(!is.na(fg), mean_n_individuals > 10), aes(x = bin_midpoint, y = median, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = 'Light received per unit crown area (Wm-2)') +
  scale_y_log10(name = 'Individual production (kg y-1)') +
  theme_bw() +
  ggtitle('Production (individual)') 

p3fits <- ggplot(lightproductionbins_fg %>% filter(!is.na(fg), bin_count > 10)) +
  geom_ribbon(data = light_fitted_totalprod %>% filter(!dens_model %in% '1', prod_model == 2, !fg %in% 'unclassified') %>% mutate(combo = paste(dens_model,prod_model,sep='x')), aes(x = lightperarea, ymin = q025, ymax = q975, group = combo), fill = 'gray80') +
  geom_line(data = light_fitted_totalprod %>% filter(!dens_model %in% '1', prod_model == 2, !fg %in% 'unclassified') %>% mutate(combo = paste(dens_model,prod_model,sep='x')), aes(x = lightperarea, y = q50, group = combo, color = combo)) +
  geom_point(aes(x = bin_midpoint, y = bin_value)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = 'Light received per unit crown area (Wm-2)') +
  scale_y_log10(name = 'Total production (kg y-1 ha-1 per watt)', limits = c(1e-4, 5e1)) +
  theme_bw() +
  ggtitle('Production (total)')

fpfig <- '~/google_drive/ForestLight/figs/lightpowerlaws_feb2019'

ggsave(file.path(fpfig, 'withfits_densitybylight_separate.png'), p1fits, height = 5, width = 7.5, dpi = 300)
ggsave(file.path(fpfig, 'withfits_productionindividualbylight_separate.png'), p2fits, height = 5, width = 7.5, dpi = 300)
ggsave(file.path(fpfig, 'withfits_productiontotalbylight_separate.png'), p3fits, height = 5, width = 7.5, dpi = 300)
