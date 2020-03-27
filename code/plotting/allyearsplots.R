# Plots for appendix with data from all years
# QDR / Forestlight / 26 March 2020

# This will include plots of the slope parameters and the actual fits superimposed on the data.
# For the fits, just use the point estimate of the parameters to draw the plots.


# Load data ---------------------------------------------------------------

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users/jgradym/Google Drive/ForestLight'))

library(tidyverse)
library(forestscaling)

# # Load model fits
# load('~/Dropbox/Q/projects/forestlight/fitsallyears.RData')

# Load extracted parameters and slopes
prodpars_byfg_df <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/allyearfits_prodpars.csv'))
denspars_byfg_df <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/allyearfits_denspars.csv'))

allpars_wide <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/allyearfits_allslopes.csv'))

# We also need to load the observed values for the 5 census periods, so that we can plot the binned points.
obs_dens <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_dens.csv')) %>%
  filter(!fg %in% 'unclassified') %>%
  mutate(fg = if_else(fg == 'all', 'alltree', fg))
obs_indivprod <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_indivprod.csv')) %>%
  filter(!fg %in% 'unclassified') %>%
  mutate(fg = if_else(fg == 'all', 'alltree', fg))
obs_totalprod <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_totalprod.csv')) %>%
  filter(!fg %in% 'unclassified') %>%
  mutate(fg = if_else(fg == 'all', 'alltree', fg))


# Themes
theme_set(theme_bw() +
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.position = 'bottom'))

guild_colors <- c("black", "#99AF3C", "#1D5128", "#1e3160", "#6BA3BF", "#D6D6D6")
guild_fills <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "ivory")
fg_names <-  c('all', 'fast', 'tall', 'slow', 'short', 'medium')

fg_label <- as_labeller(setNames(fg_names, unique(obs_dens$fg)))

# Slopes across all years -------------------------------------------------



plotpars_long <- allpars_wide %>%
  pivot_longer(cols = -c(year, fg), names_to = 'parameter') %>%
  filter(parameter %in% c('indivprod_slope', 'dens_slope_mid', 'totalprod_slope_mid'), !fg %in% 'unclassified') %>%
  mutate(parameter = factor(parameter, labels = c('abundance', 'individual growth', 'total production')))

# Slopes for midsize trees, for all 5 censuses.
p_slopes <- ggplot(plotpars_long, aes(x = year, y = value, color = fg, fill = fg, shape = parameter)) +
  geom_hline(yintercept = c(-2, 0, 2), linetype = 'dashed', alpha = 0.5) +
  geom_line(show.legend = FALSE) +
  geom_point(color = 'black', size = 2) +
  scale_color_manual(values = guild_colors) +
  scale_shape_manual(values = c(21, 22, 23), guide = FALSE) +
  scale_fill_manual(values = guild_fills, name = 'Functional guild', labels = fg_names, guide = guide_legend(override.aes = list(shape = 21))) +
  scale_y_continuous(name = 'slope')


# Fitted and observed plots -----------------------------------------------

# Plot just the fitted values

prodpars_wide <- prodpars_byfg_df %>%
  pivot_wider(names_from = parameter)

denspars_wide <- denspars_byfg_df %>%
  pivot_wider(names_from = parameter)

# Evaluate production and density functions for each of the Fgs in each year at many points
dbh_pred <- logseq(1, 300, 101)

prod_fitted <- prodpars_wide %>%
  group_by(year, fg) %>%
  group_modify(~ data.frame(x = dbh_pred, y = powerlaw_log(dbh_pred, .$beta0, .$beta1))) %>%
  filter(!fg %in% 'unclassified')

# Density needs to be corrected by multiplying by number of individuals and dividing by area
dens_multipliers <- obs_dens %>% group_by(fg, year) %>%
  summarize(multiplier = sum(bin_count)/42.84)

dens_fitted <- denspars_wide %>%
  group_by(year, fg) %>%
  group_modify(~ data.frame(x = dbh_pred, y = pdf_3part(x = dbh_pred, xmin = 1, alpha_low = .$alpha_low, alpha_mid = .$alpha_mid, alpha_high = .$alpha_high, tau_low = .$tau_low, tau_high = .$tau_high))) %>%
  filter(!fg %in% 'unclassified') %>%
  left_join(dens_multipliers) %>%
  mutate(y = y * multiplier)

p_dens <- ggplot() +
  geom_point(data = obs_dens %>% filter(bin_count >= 20, bin_value > 0), aes(x = bin_midpoint, y = bin_value, group = factor(year), color = factor(year))) +
  geom_line(data = dens_fitted, aes(x = x, y = y, group = factor(year), color = factor(year))) +
  facet_wrap(~ fg, labeller = fg_label, nrow = 3) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = parse(text = 'Abundance~(ha^-1~cm^-1)'), limits = c(1e-3, 2500)) +
  scale_color_viridis_d(name = 'year', direction = -1)

p_prod <- ggplot() +
  geom_point(data = obs_indivprod %>% filter(mean_n_individuals >= 20), aes(x = bin_midpoint, y = median, group = factor(year), color = factor(year))) +
  geom_line(data = prod_fitted, aes(x = x, y = y, group = factor(year), color = factor(year))) +
  facet_wrap(~ fg, labeller = fg_label, nrow = 3) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = parse(text='Individual~Growth~(kg~y^-1)')) + 
  scale_color_viridis_d(name = 'year', direction = -1)


# Z-score plots -----------------------------------------------------------

# z scores by year
allpars_long <- allpars_wide %>%
  pivot_longer(cols = -c(year, fg), names_to = 'parameter') 

zscores <- allpars_long %>%
  group_by(fg, parameter) %>%
  mutate(z = (value - mean(value))/sd(value))

p_zscores <- ggplot(zscores, aes(x = factor(year), y = abs(z), fill = factor(year))) + 
  geom_boxplot(color = 'gray30', size = 1) +
  scale_fill_viridis_d(direction = -1) +
  theme(legend.position = 'none') +
  labs(x = 'year', y = 'parameter z-score')


# Save plots --------------------------------------------------------------

ggsave(file.path(gdrive_path, 'figs/plots_all_years/slopes_all_years.png'), p_slopes, height = 6, width = 5, dpi = 300)
ggsave(file.path(gdrive_path, 'figs/plots_all_years/density_all_years.png'), p_dens, height = 7, width = 5, dpi = 300)
ggsave(file.path(gdrive_path, 'figs/plots_all_years/production_all_years.png'), p_prod, height = 7, width = 5, dpi = 300)
ggsave(file.path(gdrive_path, 'figs/plots_all_years/zscores_all_years.png'), p_zscores, height = 5, width = 5, dpi = 300)
