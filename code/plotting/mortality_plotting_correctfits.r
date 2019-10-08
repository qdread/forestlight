# Plot mortality binned data (binned in the mortality_processing.r script)
# this version uses the correct mixed model fit lines
# QDR / Forestlight / 08 Oct 2019

# Load data ---------------------------------------------------------------

library(tidyverse)

# Check out this slick trick so that we no longer need to comment out any paths - just run this and it sees if it's Quentin or not.
user <- Sys.info()['user']
gdrive_path <- ifelse(user == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))
github_path <- ifelse(user == 'qread', '~/Documents/GitHub/forestlight', file.path('/Users',user,'Documents/GitHub/forestlight'))

bin_mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting_aug2018/obs_mortalitybins.csv')) # Load binned data
mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting_aug2018/obs_mortalityindividuals.csv')) # Load raw data
fitted_mort <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/mortality_ci_by_fg.csv')) # Load fitted values

# Make some simple plots --------------------------------------------------

# Color mapping
fg_labels <- c('Fast','LL Pioneer', 'Slow', 'SL Breeder', 'Medium')
guild_fills_nb <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")
guild_colors_nb <- c("#3B4403", "#02330A", "#031a49", "#02394F", "#595A5B")

plightarea <- ggplot(fitted_mort %>% mutate(fg = factor(fg, labels = fg_labels))) +
  geom_ribbon(aes(x = light_per_area, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
  geom_line(aes(x = light_per_area, y = q50, group = fg, color = fg)) +
  geom_point(data = bin_mort %>% filter(variable == 'light_per_area', !fg %in% c('all','unclassified'), (lived+died) > 7)  %>% mutate(fg = factor(fg, labels = fg_labels)),
             aes(x = bin_midpoint, y = mortality, fill = fg, size = lived + died),
             pch = 21) +
  scale_x_log10(name = parse(text = 'Light~received~per~unit~crown~area~(W~m^-2)'), limits = c(1, 412), expand = c(0, 0)) +
  scale_y_continuous(name = 'Mortality rate 1990-1995') +
  scale_size_continuous(name = 'Individuals', trans = 'log', breaks = 10^(0:3)) +
  scale_color_manual(values = guild_fills_nb) +
  scale_fill_manual(values = guild_fills_nb) +
  theme_bw() +
  ggtitle('Mortality rate by light per area')
