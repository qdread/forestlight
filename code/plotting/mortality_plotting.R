# Plot mortality binned data (binned in the mortality_processing.r script)
# QDR / Forestlight / 03 Oct 2019

# Load data ---------------------------------------------------------------

library(tidyverse)

# Check out this slick trick so that we no longer need to comment out any paths - just run this and it sees if it's Quentin or not.
user <- Sys.info()['user']
gdrive_path <- ifelse(user == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))
github_path <- ifelse(user == 'qread', '~/Documents/GitHub/forestlight', file.path('/Users',user,'Documents/GitHub/forestlight'))

bin_mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_mortalitybins.csv')) # Load binned data
mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_mortalityindividuals.csv')) # Load raw data


# Make some simple plots --------------------------------------------------

# Color mapping
fg_labels <- c('Fast','LL Pioneer', 'Slow', 'SL Breeder', 'Medium')
guild_fills_nb <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")
guild_colors_nb <- c("#3B4403", "#02330A", "#031a49", "#02394F", "#595A5B")



# Show functional groups
ggplot(bin_mort %>% filter(variable == 'light_per_area', !fg %in% c('all','unclassified')) %>% mutate(fg = factor(fg, labels = fg_labels)), 
       aes(x = bin_midpoint, y = mortality, color = fg)) +
  geom_point(aes(size = lived + died)) +
  scale_x_log10(name = parse(text = 'Light~received~per~unit~crown~area~(W~m^-2)')) +
  scale_y_continuous(name = 'Mortality rate 1990-1995') +
  scale_size_continuous(name = 'Individuals', trans = 'log', breaks = 10^(0:4)) +
  theme_bw() +
  ggtitle('Mortality rate by light per area')

# Plot the logistic curves for each functional group
# Only use light as the predictor, not dbh
# Also do a very simple model fit (not even Bayesian!)

# Fit the model outside ggplot to make sure it's correct
mort_lightarea_fits <- mort %>% 
  filter(!fg %in% 'unclassified') %>%
  mutate(dead = as.numeric(!alive)) %>%
  group_by(fg) %>%
  do(fit = glm(dead ~ log10(light_received_byarea), family = 'binomial', data = .))


ggplot(mort %>% filter(!fg %in% c('unclassified'))  %>% mutate(fg = factor(fg, labels = fg_labels)),
       aes(x = light_received_byarea, y = as.numeric(!alive), group = fg, color = fg)) +
  geom_smooth(method = 'glm', method.args = list(family = "binomial")) +
  scale_x_log10(name = parse(text = 'Light~received~per~unit~crown~area~(W~m^-2)')) +
  scale_y_continuous(name = 'Mortality rate 1990-1995') +
  theme_bw()

# Combine
plightarea <- ggplot() +
  geom_smooth(data = mort %>% filter(!fg %in% c('unclassified')) %>% mutate(fg = factor(fg, labels = fg_labels)),
              aes(x = light_received_byarea, y = as.numeric(!alive), group = fg, color = fg),
              method = 'glm', method.args = list(family = "binomial"),
              alpha = 0.3, lwd = 1.5) +
  geom_point(data = bin_mort %>% filter(variable == 'light_per_area', !fg %in% c('all','unclassified'), (lived+died) > 7)  %>% mutate(fg = factor(fg, labels = fg_labels)),
             aes(x = bin_midpoint, y = mortality, fill = fg, size = lived + died),
             pch = 21) +
  scale_x_log10(name = parse(text = 'Light~received~per~unit~crown~area~(W~m^-2)')) +
  scale_y_continuous(name = 'Mortality rate 1990-1995') +
  scale_size_continuous(name = 'Individuals', trans = 'log', breaks = 10^(0:3)) +
  scale_color_manual(values = guild_fills_nb) +
  scale_fill_manual(values = guild_fills_nb) +
  theme_bw() +
  ggtitle('Mortality rate by light per area')



# Same plot by diameter and light per volume ------------------------------

pdiam <- ggplot() +
  geom_smooth(data = mort %>% filter(!fg %in% c('unclassified')) %>% mutate(fg = factor(fg, labels = fg_labels)),
              aes(x = dbh, y = as.numeric(!alive), group = fg, color = fg),
              method = 'glm', method.args = list(family = "binomial"),
              alpha = 0.3, lwd = 1.5) +
  geom_point(data = bin_mort %>% filter(variable == 'dbh', !fg %in% c('all','unclassified'), (lived+died) > 5)  %>% mutate(fg = factor(fg, labels = fg_labels)),
             aes(x = bin_midpoint, y = mortality, fill = fg, size = lived + died),
             pch = 21) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_continuous(name = 'Mortality rate 1990-1995') +
  scale_size_continuous(name = 'Individuals', trans = 'log', breaks = 10^(0:3)) +
  scale_color_manual(values = guild_fills_nb) +
  scale_fill_manual(values = guild_fills_nb) +
  theme_bw() +
  ggtitle('Mortality rate by size')

plightvolume <- ggplot() +
  geom_smooth(data = mort %>% filter(!fg %in% c('unclassified')) %>% mutate(fg = factor(fg, labels = fg_labels)),
              aes(x = light_received_byvolume, y = as.numeric(!alive), group = fg, color = fg),
              method = 'glm', method.args = list(family = "binomial"),
              alpha = 0.3, lwd = 1.5) +
  geom_point(data = bin_mort %>% filter(variable == 'light_per_volume', !fg %in% c('all','unclassified'), lived+died > 5)  %>% mutate(fg = factor(fg, labels = fg_labels)),
             aes(x = bin_midpoint, y = mortality, fill = fg, size = lived + died),
             pch = 21) +
  scale_x_log10(name = parse(text = 'Light~received~per~unit~crown~volume~(W~m^-3)')) +
  scale_y_continuous(name = 'Mortality rate 1990-1995') +
  scale_size_continuous(name = 'Individuals', trans = 'log', breaks = 10^(0:3)) +
  scale_color_manual(values = guild_fills_nb) +
  scale_fill_manual(values = guild_fills_nb) +
  theme_bw() +
  ggtitle('Mortality rate by light per volume')


# Write plots to pdf ------------------------------------------------------

pdf(file.path(gdrive_path, 'figs/mortality_bubblecharts.pdf'), height = 6, width = 7)
  plightarea; pdiam; plightvolume
dev.off()


# Diameter plot with logit scale ------------------------------------------

ysc <- scale_y_continuous(breaks = c(0.03, 0.3, .1), labels = c(0.03, 0.3, 0.1), limits = c(0.02, .5),
                   name = expression(paste("Mortality (5 yr"^-1,")")), trans = 'logit')

pdiamlogit <- ggplot() +
  # geom_smooth(data = mort %>% filter(!fg %in% c('unclassified')) %>% mutate(fg = factor(fg, labels = fg_labels)),
  #             aes(x = dbh, y = as.numeric(!alive), group = fg, color = fg),
  #             method = 'glm', method.args = list(family = "binomial"),
  #             alpha = 0.3, lwd = 1.5) +
  geom_point(data = bin_mort %>% filter(variable == 'dbh', !fg %in% c('all','unclassified'), (lived+died) > 10, died > 0)  %>% mutate(fg = factor(fg, labels = fg_labels)),
             aes(x = bin_midpoint, y = mortality, fill = fg, size = lived + died),
             pch = 21) +
  scale_x_log10(name = 'Diameter (cm)') +
  ysc +
  scale_size_continuous(name = 'Individuals', trans = 'log', breaks = 10^(0:3)) +
  scale_color_manual(values = guild_fills_nb) +
  scale_fill_manual(values = guild_fills_nb) +
  theme_bw() +
  ggtitle('Mortality rate by size')

ggsave(file.path(gdrive_path, 'figs/mortality_diameter_logity.pdf'), pdiamlogit, height = 5, width = 6)
