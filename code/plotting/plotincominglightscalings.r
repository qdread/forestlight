# Plots of "raw light scaling" aka incoming energy scaling
# Fitted and predicted values, as well as comparing slopes to the production scaling
# Use 1995 data

library(tidyverse)

### 
# Fitted slopes
rawlightslopes <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/rawlightpiecewise/rawlightpiecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE)
rawlightslopes$fg <- factor(rawlightslopes$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))
rawlightslopes$variable <- factor(rawlightslopes$variable, labels = c("Density", "Individual Incoming Energy", "Total Incoming Energy"))
colors <- c("sienna4", "yellowgreen", "springgreen4")

# Using 3 segment density and 2 segment production
ggplot(rawlightslopes %>% filter(dens_model == 3, prod_model == 2, !fg %in% 'Unclassified'), 
       aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg,labeller = label_value) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = "springgreen4", size = 0.3) +
  geom_hline(yintercept = -2, linetype = 'dashed', col = "sienna4", size = 0.3) +
  geom_hline(yintercept = 2, linetype = 'dashed', col = "yellowgreen", size = 0.3) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_ribbon(alpha = 0.4, size = 0.3) +
  geom_line() +
  scale_x_log10(name = 'Diameter (cm)', expand = c(0,0)) +
  theme_bw() + theme(axis.text = element_text(color = "black"))+
  theme(strip.background = element_blank(), panel.grid = element_blank(), legend.position = 'bottom') +
  labs(y = 'Fitted Slope') +  #coord_fixed(ratio = .1)+
  ggtitle('Fitted slopes \n 3 segment density model and 2 segment incoming energy model')+
  theme(plot.title = element_text(hjust=0.5)) #+ theme(legend.title=element_blank())

slopes <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE)
slopes$fg <- factor(slopes$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))
slopes$variable <- factor(slopes$variable, labels = c("Density", "Individual Growth", "Total Growth"))
colors <- c("sienna4", "yellowgreen", "springgreen4")

# Using 3 segment density and 2 segment production
ggplot(slopes %>% filter(dens_model == 3, prod_model == 1, !fg %in% 'Unclassified'), 
       aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg,labeller = label_value) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = "springgreen4", size = 0.3) +
  geom_hline(yintercept = -2, linetype = 'dashed', col = "sienna4", size = 0.3) +
  geom_hline(yintercept = 2, linetype = 'dashed', col = "yellowgreen", size = 0.3) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_ribbon(alpha = 0.4, size = 0.3) +
  geom_line() +
  scale_x_log10(name = 'Diameter (cm)', expand = c(0,0)) +
  theme_bw() + theme(axis.text = element_text(color = "black"))+
  theme(strip.background = element_blank(), panel.grid = element_blank(), legend.position = 'bottom') +
  labs(y = 'Fitted Slope') +  #coord_fixed(ratio = .1)+
  ggtitle('Fitted slopes \n 3 segment density model and 2 segment production model')+
  theme(plot.title = element_text(hjust=0.5)) #+ theme(legend.title=element_blank())

# Get the slope + ci for the middle portion
prod_slopes_10cm <- slopes %>% filter(variable == 'Total Growth', dbh == unique(slopes$dbh)[38], dens_model == 3, prod_model == 2)
light_slopes_10cm <- rawlightslopes %>% filter(variable == 'Total Incoming Energy', dbh == unique(rawlightslopes$dbh)[38], dens_model == 3, prod_model == 2)

ggplot(rbind(prod_slopes_10cm, light_slopes_10cm) %>% filter(!fg %in% c('All','Unclassified')), aes(x = fg, y = q50, ymin = q025, ymax = q975, group = variable, color = variable)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_errorbar(position = position_dodge(width = 0.1), width = 0.2) +
  geom_point(position = position_dodge(width = 0.1)) +
  theme_classic() +
  ggtitle('Symmetry of total growth and total incoming energy patterns', 'Slope of scaling relationship at 10 cm dbh')

###
# Observed values with fitted lines superimposed (update to Figure 4)

# Load fitted values
rawlight_ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/rawlightpiecewise/rawlightpiecewise_ci_by_fg.csv', stringsAsFactors = FALSE) 
area_core <- 42.84

rawlight_ci_df$fg[rawlight_ci_df$fg == 'alltree'] <- 'all'

rawlight_pred_dens <- rawlight_ci_df %>%
  filter(variable == 'density') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

rawlight_fitted_indivprod <- rawlight_ci_df %>%
  filter(variable == 'production_fitted') %>%
  select(-variable)

rawlight_fitted_totalprod <- rawlight_ci_df %>%
  filter(variable == 'total_production_fitted') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

rawlight_pred_indivprod <- rawlight_ci_df %>%
  filter(variable == 'production') %>%
  select(-variable)

rawlight_pred_totalprod <- rawlight_ci_df %>%
  filter(variable == 'total_production') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

# Load observed values
rawlight_obs_totalprod <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/rawlightpiecewise/observed_light_bin.csv', stringsAsFactors = FALSE) %>% mutate(bin_value = bin_value/area_core)


# The below df is from the "total workflow" (Do not run)
# rawlight_observed_bin <- cbind(fg = 'all', lightreceivedbin_alltree_byyear[[2]]) %>%
#   rbind(map2_dfr(lightreceivedbin_fg_byyear, c('fg1','fg2','fg3','fg4','fg5','unclassified'), ~ cbind(fg = .y, .x[[2]])))
# write.csv(rawlight_observed_bin, file = '~/google_drive/ForestLight/data/data_piecewisefits/rawlightpiecewise/observed_light_bin.csv', row.names = FALSE)

# Modify fitted values so that they do not go above the data.
bin_maxes <- rawlight_obs_totalprod %>% group_by(fg) %>% summarize(max = max(bin_midpoint[bin_count > 10]))
rawlight_fitted_totalprod <- rawlight_fitted_totalprod %>%
  left_join(bin_maxes) %>%
  filter(dbh <= max)

prawlightfits <- ggplot(rawlight_obs_totalprod %>% filter(!is.na(fg), !fg %in% 'unclassified', bin_count > 10)) +
  geom_ribbon(data = rawlight_fitted_totalprod %>% filter(prod_model == 2, !fg %in% 'unclassified'), aes(x = dbh, ymin = q025, ymax = q975), fill = 'gray80') +
  geom_line(data = rawlight_fitted_totalprod %>% filter(prod_model == 2, !fg %in% 'unclassified'), aes(x = dbh, y = q50)) +
  geom_point(aes(x = bin_midpoint, y = bin_value)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = expression(paste('Total incoming light (W ha'^-1,')', sep=''))) +
  theme_bw() 

mycols <- c('black', RColorBrewer::brewer.pal(5,'Set1'))

prawlightfits_onefig <- ggplot(rawlight_obs_totalprod %>% filter(!is.na(fg), !fg %in% 'unclassified', bin_count > 10)) +
  geom_ribbon(data = rawlight_fitted_totalprod %>% filter(prod_model == 2, !fg %in% 'unclassified'), aes(x = dbh, ymin = q025, ymax = q975, fill = fg), alpha = 0.4) +
  geom_line(data = rawlight_fitted_totalprod %>% filter(prod_model == 2, !fg %in% 'unclassified'), aes(x = dbh, y = q50, color = fg)) +
  geom_point(aes(x = bin_midpoint, y = bin_value, color = fg)) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = expression(paste('Total incoming light (W ha'^-1,')', sep='')), limits = c(1e2,1e5)) +
  theme_bw() +
  scale_color_manual(values=mycols) + scale_fill_manual(values=mycols)
  

fpfig <- '~/google_drive/ForestLight/figs/lightpowerlaws_feb2019'
ggsave(file.path(fpfig, 'totallightscaling_separateplots.png'), prawlightfits, height = 5, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'totallightscaling_oneplot.png'), prawlightfits_onefig, height = 5, width = 6, dpi = 300)
