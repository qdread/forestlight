# Plot light received by crown area and volume, with individual and binned production
# Volume added 23 Mar 2019
# Growth/volume vs. light/volume added 25 Mar 2019

gdrive_path <- '~/google_drive/ForestLight'
github_path <- '~/Documents/GitHub/forestlight'

gdrive_path <- '/Users/jgradym/Google Drive/ForestLight'
github_path <- '/Users/jgradym/Documents/GitHub/forestlight'

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

library(dplyr)
library(ggplot2)
library(forestscaling)

# Create bins
# -----------
area_core <- 42.84
num_bins <- 20

log_midpoints <- function(a) exp(log(a)[-length(a)] + diff(log(a))/2)

lightperarearange <- with(alltree_light_95, range(light_received/crownarea)) # 1.04 to 412
lightperarearange[2] <- lightperarearange[2] + 1

lightperareabinbounds <- exp(seq(log(lightperarearange[1]), log(lightperarearange[2]), length.out = num_bins + 1))

lightperareabinedges <- data.frame(bin_min = lightperareabinbounds[-(num_bins+1)], bin_max = lightperareabinbounds[-1], bin_midpoint = log_midpoints(lightperareabinbounds))

lightpervolumerange <- with(alltree_light_95, range(light_received/crownvolume)) # 2.30 to 316
lightpervolumerange[2] <- lightpervolumerange[2] + 1

lightpervolumebinbounds <- exp(seq(log(lightpervolumerange[1]), log(lightpervolumerange[2]), length.out = num_bins + 1))

lightpervolumebinedges <- data.frame(bin_min = lightpervolumebinbounds[-(num_bins+1)], bin_max = lightpervolumebinbounds[-1], bin_midpoint = log_midpoints(lightpervolumebinbounds))


alltree_light_95 <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), as.character(NA)))

# Area bins
# ---------
lightabundbins_all <- with(alltree_light_95, logbin(x = light_received/crownarea, n = 20))
lightproductionbins_all <- with(alltree_light_95, logbin_setedges(x = light_received/crownarea, y = production, edges = lightabundbins_all))
lightindivprodbins_all <- alltree_light_95 %>%
  mutate(lightperarea_bin = cut(light_received/crownarea, breaks = c(lightperareabinedges$bin_min[1], lightperareabinedges$bin_max), labels = lightperareabinedges$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(lightperarea_bin) %>%
  do(c(n = nrow(.), quantile(.$production, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))

lightabundbins_fg <- alltree_light_95 %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$light_received/.$crownarea, edges = lightabundbins_all))

lightproductionbins_fg <- alltree_light_95 %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$light_received/.$crownarea, y = .$production, edges = lightabundbins_all))

lightindivprodbins_fg <- alltree_light_95 %>%
  mutate(lightperarea_bin = cut(light_received/crownarea, breaks = c(lightperareabinedges$bin_min[1], lightperareabinedges$bin_max), labels = lightperareabinedges$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, lightperarea_bin) %>%
  do(c(n = nrow(.), quantile(.$production, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))


lightabundbins_fg <- rbind(data.frame(fg = 'all', lightabundbins_all, stringsAsFactors = FALSE), as.data.frame(lightabundbins_fg)) %>%
  mutate(bin_value = bin_value / area_core)
lightproductionbins_fg <- rbind(data.frame(fg = 'all', lightproductionbins_all, stringsAsFactors = FALSE), as.data.frame(lightproductionbins_fg)) %>%
  mutate(bin_value = bin_value / area_core)
lightindivprodbins_fg <- data.frame(fg = 'all', lightindivprodbins_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightindivprodbins_fg)) %>%
  mutate(lightperarea_bin = as.numeric(as.character(lightperarea_bin))) %>%
  rename(bin_midpoint = lightperarea_bin)
           

# Volume bins
# -----------
lightpervolabundbins_all <- with(alltree_light_95, logbin(x = light_received/crownvolume, n = 20))
lightpervolproductionbins_all <- with(alltree_light_95, logbin_setedges(x = light_received/crownvolume, y = production, edges = lightpervolabundbins_all))
lightpervolindivprodbins_all <- alltree_light_95 %>%
  mutate(lightpervol_bin = cut(light_received/crownvolume, breaks = c(lightpervolumebinedges$bin_min[1], lightpervolumebinedges$bin_max), labels = lightpervolumebinedges$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(lightpervol_bin) %>%
  do(c(n = nrow(.), quantile(.$production, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))

lightpervolabundbins_fg <- alltree_light_95 %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$light_received/.$crownvolume, edges = lightpervolabundbins_all))

lightpervolproductionbins_fg <- alltree_light_95 %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$light_received/.$crownvolume, y = .$production, edges = lightpervolabundbins_all))

lightpervolindivprodbins_fg <- alltree_light_95 %>%
  mutate(lightpervol_bin = cut(light_received/crownvolume, breaks = c(lightpervolumebinedges$bin_min[1], lightpervolumebinedges$bin_max), labels = lightpervolumebinedges$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, lightpervol_bin) %>%
  do(c(n = nrow(.), quantile(.$production, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))


lightpervolabundbins_fg <- rbind(data.frame(fg = 'all', lightpervolabundbins_all, stringsAsFactors = FALSE), as.data.frame(lightpervolabundbins_fg)) %>%
  mutate(bin_value = bin_value / area_core)
lightpervolproductionbins_fg <- rbind(data.frame(fg = 'all', lightpervolproductionbins_all, stringsAsFactors = FALSE), as.data.frame(lightpervolproductionbins_fg)) %>%
  mutate(bin_value = bin_value / area_core)
lightpervolindivprodbins_fg <- data.frame(fg = 'all', lightpervolindivprodbins_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightpervolindivprodbins_fg)) %>%
  mutate(lightpervol_bin = as.numeric(as.character(lightpervol_bin))) %>%
  rename(bin_midpoint = lightpervol_bin)

# Added 25 Mar: production *per unit volume* as binned quantity
lightpervolprodpervolbins_all <- alltree_light_95 %>%
  mutate(lightpervol_bin = cut(light_received/crownvolume, breaks = c(lightpervolumebinedges$bin_min[1], lightpervolumebinedges$bin_max), labels = lightpervolumebinedges$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(lightpervol_bin) %>%
  do(c(n = nrow(.), quantile(.$production/.$crownvolume, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))
lightpervolprodpervolbins_fg <- alltree_light_95 %>%
  mutate(lightpervol_bin = cut(light_received/crownvolume, breaks = c(lightpervolumebinedges$bin_min[1], lightpervolumebinedges$bin_max), labels = lightpervolumebinedges$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, lightpervol_bin) %>%
  do(c(n = nrow(.), quantile(.$production/.$crownvolume, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))
lightpervolprodpervolbins_fg <- data.frame(fg = 'all', lightpervolprodpervolbins_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightpervolprodpervolbins_fg)) %>%
  mutate(lightpervol_bin = as.numeric(as.character(lightpervol_bin))) %>%
  rename(bin_midpoint = lightpervol_bin)

# Area plots --------------------------------------------------------------

# Expressions for axis labels
exlpa <- expression(paste('Light received per unit crown area (W m'^-2,')', sep = ''))
exlpv <- expression(paste('Light received per unit crown volume (W m'^-3,')', sep = ''))
exd <- expression(paste('Density (trees ha'^-1,')', sep = ''))
exindp <- expression(paste('Growth (kg y'^-1, ')', sep = ''))
extotp <- expression(paste('Growth (kg y'^-1, ' ha'^-1,  ')', sep = ''))


p1 <- ggplot(lightabundbins_fg %>% filter(!is.na(fg)), aes(x = bin_midpoint, y = bin_value)) +
  geom_point() +
  facet_wrap(~ fg) +
  scale_x_log10(name = exlpa) +
  scale_y_log10(name = exd) +
  theme_bw() +
  ggtitle('Density')
p1  

p1b <- ggplot(lightabundbins_fg %>% filter(!is.na(fg), !fg %in% 'all'), aes(x = bin_midpoint, y = bin_value, color = fg)) +
  geom_point() +
  scale_x_log10(name = exlpa) +
  scale_y_log10(exd) +
  theme_bw() +
  ggtitle('Density')
p1b

p2 <- ggplot(alltree_light_95 %>% filter(!is.na(fg)), aes(x = light_received/crownarea, y = production)) +
  geom_hex(alpha = 0.4) +
  geom_hex(data = alltree_light_95 %>% mutate(fg = 'all'), alpha = 0.4) +
  geom_pointrange(data = lightindivprodbins_fg %>% filter(!is.na(fg)), aes(x = bin_midpoint, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = exlpa) +
  scale_y_log10(exindp) +
  theme_bw() +
  ggtitle('Production (individual)') +
  scale_fill_gradient(trans = 'log', low = 'skyblue', high = 'darkblue', breaks = c(1, 10, 100))
p2

p2b <- ggplot(lightindivprodbins_fg %>% filter(!is.na(fg), !fg %in% 'all'), aes(x = bin_midpoint, y = q50, ymin = q25, ymax = q75, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = exlpa) +
  scale_y_log10(exindp) +
  theme_bw() +
  ggtitle('Production (individual)')
p2b

p3 <- ggplot(lightproductionbins_fg %>% filter(!is.na(fg)), aes(x = bin_midpoint, y = bin_value)) +
  geom_point() +
  facet_wrap(~ fg) +
  scale_x_log10(name = exlpa) +
  scale_y_log10(extotp) +
  theme_bw() +
  ggtitle('Production (total)')
p3

p3b <- ggplot(lightproductionbins_fg %>% filter(!is.na(fg), !fg %in% 'all'), aes(x = bin_midpoint, y = bin_value, color = fg)) +
  geom_point() +
  scale_x_log10(name = exlpa) +
  scale_y_log10(extotp) +
  theme_bw() +
  ggtitle('Production (total)')
p3b

fpfig <- file.path(gdrive_path, '/figs/lightpowerlaws_feb2019')
ggsave(file.path(fpfig, 'densitybylight_separate.png'), p1, height = 5, width = 9, dpi = 300)
ggsave(file.path(fpfig, 'densitybylight_together.png'), p1b, height = 5, width = 5, dpi = 300)
ggsave(file.path(fpfig, 'productionindividualbylight_separate.png'), p2, height = 5, width = 9, dpi = 300)
ggsave(file.path(fpfig, 'productionindividuallight_together.png'), p2b, height = 5, width = 5, dpi = 300)
ggsave(file.path(fpfig, 'productiontotalbylight_separate.png'), p3, height = 5, width = 9, dpi = 300)
ggsave(file.path(fpfig, 'productiontotalbylight_together.png'), p3b, height = 5, width = 5, dpi = 300)


# Area plots with fitted lines --------------------------------------------------------

light_ci_df <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/lightpiecewise/lightpiecewise_ci_by_fg.csv'), stringsAsFactors = FALSE) 
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
  scale_x_log10(name = exlpa) +
  scale_y_log10(name = exd, limits = c(1e-5, 1e2)) +
  theme_bw() +
  ggtitle('Density')
p1fits
p2fits <- ggplot() +
  geom_ribbon(data = light_fitted_indivprod %>% filter(dens_model == '1', !fg %in% 'unclassified'), aes(x = lightperarea, ymin = q025, ymax = q975, group = factor(prod_model)), fill = 'gray80') +
  geom_line(data = light_fitted_indivprod %>% filter(dens_model == '1', !fg %in% 'unclassified'), aes(x = lightperarea, y = q50, group = factor(prod_model), color = factor(prod_model))) +
  geom_pointrange(data = lightindivprodbins_fg %>% filter(!is.na(fg), mean_n_individuals > 10), aes(x = bin_midpoint, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = exlpa) +
  scale_y_log10(name = exindp) +
  theme_bw() +
  ggtitle('Production (individual)') 
p2fits

p3fits <- ggplot(lightproductionbins_fg %>% filter(!is.na(fg), bin_count > 10)) +
  geom_ribbon(data = light_fitted_totalprod %>% filter(!dens_model %in% '1', prod_model == 2, !fg %in% 'unclassified') %>% mutate(combo = paste(dens_model,prod_model,sep='x')), aes(x = lightperarea, ymin = q025, ymax = q975, group = combo), fill = 'gray80') +
  geom_line(data = light_fitted_totalprod %>% filter(!dens_model %in% '1', prod_model == 2, !fg %in% 'unclassified') %>% mutate(combo = paste(dens_model,prod_model,sep='x')), aes(x = lightperarea, y = q50, group = combo, color = combo)) +
  geom_point(aes(x = bin_midpoint, y = bin_value)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = exlpa) +
  scale_y_log10(name = extotp, limits = c(1e-4, 5e1)) +
  theme_bw() +
  ggtitle('Production (total)')
p3fits 


ggsave(file.path(fpfig, 'withfits_densitybylight_separate.png'), p1fits, height = 5, width = 7.5, dpi = 300)
ggsave(file.path(fpfig, 'withfits_productionindividualbylight_separate.png'), p2fits, height = 5, width = 7.5, dpi = 300)
ggsave(file.path(fpfig, 'withfits_productiontotalbylight_separate.png'), p3fits, height = 5, width = 7.5, dpi = 300)



# Volume plots ------------------------------------------------------------

p1vol <- ggplot(lightpervolabundbins_fg %>% filter(!is.na(fg), bin_value > 0), aes(x = bin_midpoint, y = bin_value)) +
  geom_point() +
  facet_wrap(~ fg) +
  scale_x_log10(name = exlpv) +
  scale_y_log10(name = exd) +
  theme_bw() +
  ggtitle('Density')
p1vol

p1bvol <- ggplot(lightpervolabundbins_fg %>% filter(!is.na(fg), !fg %in% 'all', bin_value > 0), aes(x = bin_midpoint, y = bin_value, color = fg)) +
  geom_point() +
  scale_x_log10(name = exlpv) +
  scale_y_log10(name = exd) +
  theme_bw() +
  ggtitle('Density')
p1bvol

p2vol <- ggplot(alltree_light_95 %>% filter(!is.na(fg)), aes(x = light_received/crownvolume, y = production)) +
  geom_hex(alpha = 0.4) +
  geom_hex(data = alltree_light_95 %>% mutate(fg = 'all'), alpha = 0.4) +
  geom_pointrange(data = lightpervolindivprodbins_fg %>% filter(!is.na(fg)), aes(x = bin_midpoint, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = exlpv) +
  scale_y_log10(exindp) +
  theme_bw() +
  ggtitle('Production (individual)') +
  scale_fill_gradient(trans = 'log', low = 'skyblue', high = 'darkblue', breaks = c(1,10,100))
p2vol

p2bvol <- ggplot(lightpervolindivprodbins_fg %>% filter(!is.na(fg), !fg %in% 'all'), aes(x = bin_midpoint, y = q50, ymin = q25, ymax = q75, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = exlpv) +
  scale_y_log10(exindp) +
  theme_bw() +
  ggtitle('Production (individual)')
p2bvol

p3vol <- ggplot(lightpervolproductionbins_fg %>% filter(!is.na(fg), bin_value > 0), aes(x = bin_midpoint, y = bin_value)) +
  geom_point() +
  facet_wrap(~ fg) +
  scale_x_log10(name = exlpv) +
  scale_y_log10(name = extotp) +
  theme_bw() +
  ggtitle('Production (total)')
p3vol

p3bvol <- ggplot(lightpervolproductionbins_fg %>% filter(!is.na(fg), !fg %in% 'all', bin_value > 0), aes(x = bin_midpoint, y = bin_value, color = fg)) +
  geom_point() +
  scale_x_log10(name = exlpv) +
  scale_y_log10(name = extotp) +
  theme_bw() +
  ggtitle('Production (total)')
p3bvol

# Growth per volume vs light per volume plots (indiv) ---------------------

exppv <- expression(paste('Growth per unit crown volume (kg y'^-1, ' m'^-3,')', sep = ''))

# Plot with medians and central 50% interval, faceted
ggplot(lightpervolprodpervolbins_fg %>% filter(!is.na(fg), !fg %in% 'all'), aes(x = bin_midpoint, y = q50, ymin = q25, ymax = q75)) +
  geom_pointrange() +
  facet_wrap(~ fg) +
  geom_abline(intercept=-3, slope = 1, color ="darkgray",linetype="dashed", size=1.5)+
  scale_x_log10(name = exlpv) +
  scale_y_log10(name = exppv) +
  theme_bw() 

# Plot with medians and central 50% interval, all together
# Dodge the positions of the bars.
ggplot(lightpervolprodpervolbins_fg %>% filter(!is.na(fg), !fg %in% 'all'), aes(x = bin_midpoint, y = q50, ymin = q25, ymax = q75, color = fg)) +
  geom_pointrange(position = position_dodge(width = 0.1)) +
  scale_x_log10(name = exlpv) +
  scale_y_log10(name = exppv) +
  theme_bw() 

# Plot with medians & interval, over hexagons, faceted
ggplot(lightpervolprodpervolbins_fg %>% filter(!is.na(fg), !fg %in% 'all')) +
  geom_hex(data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = light_received/crownvolume, y = production/crownvolume)) +
  geom_pointrange(aes(x = bin_midpoint, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = exlpv) +
  scale_y_log10(name = exppv) +
  scale_fill_gradient(trans = 'log', low = 'skyblue', high = 'darkblue', breaks = c(1,10,100,1000)) +
  theme_bw() 


# Plot with medians & interval, over point cloud, faceted
ggplot(lightpervolprodpervolbins_fg %>% filter(!is.na(fg), !fg %in% 'all')) +
  geom_point(data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = light_received/crownvolume, y = production/crownvolume), color = 'skyblue', alpha = 0.05) +
  geom_pointrange(aes(x = bin_midpoint, y = q50, ymin = q25, ymax = q75)) +
  geom_abline(intercept=-3, slope = 1, color ="darkgray",linetype="dashed", size=1.5)+
  facet_wrap(~ fg) +
  scale_x_log10(name = exlpv) +
  scale_y_log10(name = exppv) +
  theme_bw() 

str(alltree_light_95)
# Very simple model
loglogregressions <- alltree_light_95 %>%
  #filter(light_received/crownvolume > 10)
  group_by(fg) %>%
  do(model = lm(log10(production/crownvolume) ~ log10(light_received/crownvolume), data = .))

lapply(loglogregressions$model, summary)

lm1 <- lm(log10(production/crownvolume) ~ log10(light_received/crownvolume), data = alltree_light_95)
summary(lm1)
