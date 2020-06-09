# Exploratory script: richness of each FG, and ratios, vs. diameter
# QDR / forestlight / 9 June 2020

library(tidyverse)
library(forestscaling)

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users/jgradym/Google Drive/ForestLight'))
github_path <- ifelse(Sys.info()['user'] == 'qread', '~/documents/github/', file.path('/Users/jgradym/Documents/Github'))

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.R'))

# Load 1995 bin data (to make sure we're using the same bin breaks as other figs)
load(file.path(gdrive_path, 'data/data_binned/bin_object_singleyear.RData'))

# Create binned richness data ---------------------------------------------

# Use existing bin bounds
bin_bounds <- fastslow_stats_bydiam_byyear[1:20, c('bin_midpoint', 'bin_min', 'bin_max')]
bin_bounds_light <- fastslow_stats_bylight_byyear[1:20, c('bin_midpoint', 'bin_min', 'bin_max')]
# Get richness of each FG by each bin

bin_x_fg <- expand_grid(fg = c(paste0('fg', 1:5), 'unclassified'), bin = 1:20) %>%
  left_join(bin_bounds %>% mutate(bin = 1:20))

bin_x_fg_light <- expand_grid(fg = c(paste0('fg', 1:5), 'unclassified'), bin = 1:20) %>%
  left_join(bin_bounds_light %>% mutate(bin = 1:20))

# 1995 data (132,982 trees)
dat <- alltreedat[[3]] %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

bin_x_fg <- bin_x_fg %>%
  cbind(pmap_dfr(bin_x_fg, function(fg, bin_min, bin_max, ...) {
    data.frame(richness = length(unique(dat$sp[dat$fg %in% fg & dat$dbh_corr >= bin_min & dat$dbh_corr < bin_max])))
  }))

bin_x_fg_light <- bin_x_fg_light %>%
  cbind(pmap_dfr(bin_x_fg_light, function(fg, bin_min, bin_max, ...) {
    data.frame(richness = length(unique(dat$sp[dat$fg %in% fg & dat$dbh_corr >= bin_min & dat$dbh_corr < bin_max])))
  }))

str(bin_x_fg)
# Quick ggplot ------------------------------------------------------------

guild_colors <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")

## richness for each group

# Arith y scale
(p_rich <- ggplot(bin_x_fg %>% filter(!fg %in% 'unclassified' & richness > 0), aes(x = bin_midpoint, y = richness, color = fg)) + 
  geom_point(size = 4) + scale_x_log10(name = 'diameter') + theme_plant() +
  scale_color_manual(values = guild_colors) +
  theme(legend.position = 'right'))

# Log y scale
(p_rich_log <- p_rich + scale_y_log10(limits = c(0.6, 140)))

## ratios of richness

richness_wide <- bin_x_fg %>%
  pivot_wider(id_cols = c(bin, bin_midpoint, bin_min, bin_max), names_from = fg, values_from = richness) %>%
  mutate(richness_ratio_fastslow = fg1/fg3,
         richness_ratio_pioneerbreeder = fg2/fg4) %>%
  mutate_if(is.double, ~ if_else(is.finite(.), ., as.numeric(NA)))

# Richness ratio 
(p_ratio <- ggplot(richness_wide, aes(x = bin_midpoint)) +
  geom_point(size = 4, color = 'gray', aes(y = richness_ratio_fastslow)) +
  geom_point(size = 4, color = 'black', aes(y = richness_ratio_pioneerbreeder)) +
  annotate(geom = 'text', x = 1, y = 30, label = 'fast:slow', size = 6, color = 'gray', hjust = 0) +
  annotate(geom = 'text', x = 1, y = 25, label = 'pioneer:breeder', size = 6, color = 'black', hjust = 0) +
  scale_x_log10(name = 'diameter') + scale_y_log10(name = 'richness ratio', limit = c(0.5, 10)) + theme_plant())

(p_ratio <- ggplot(richness_wide, aes(x = bin_midpoint)) +
    geom_point(size = 4, color = 'gray', aes(y = richness_ratio_fastslow)) +
    geom_point(size = 4, color = 'black', aes(y = richness_ratio_pioneerbreeder)) +
    annotate(geom = 'text', x = 1, y = 30, label = 'fast:slow', size = 6, color = 'gray', hjust = 0) +
    annotate(geom = 'text', x = 1, y = 25, label = 'pioneer:breeder', size = 6, color = 'black', hjust = 0) +
    scale_x_log10(name = 'diameter') + scale_y_log10(name = 'richness ratio', limit = c(0.5, 2)) + theme_plant())

#----- by light------------------------

(p_rich_light <- ggplot(bin_x_fg_light %>% filter(!fg %in% 'unclassified' & richness > 0), 
                        aes(x = bin_midpoint, y = richness, color = fg)) + 
   geom_point(size = 4) + scale_x_log10(name = 'diameter') + theme_plant() +
   scale_color_manual(values = guild_colors) +
   scale_y_log10() + #limits = c(0.6, 140)
   theme(legend.position = 'right'))


richness_wide_light <- bin_x_fg_light %>%
  pivot_wider(id_cols = c(bin, bin_midpoint, bin_min, bin_max), names_from = fg, values_from = richness) %>%
  mutate(richness_ratio_fastslow = fg1/fg3,
         richness_ratio_pioneerbreeder = fg2/fg4) %>%
  mutate_if(is.double, ~ if_else(is.finite(.), ., as.numeric(NA)))


(p_ratio_light <- ggplot(richness_wide_light, aes(x = bin_midpoint)) +
    geom_point(size = 4, color = 'gray', aes(y = richness_ratio_fastslow)) +
    geom_point(size = 4, color = 'black', aes(y = richness_ratio_pioneerbreeder)) +
    annotate(geom = 'text', x = 1, y = 30, label = 'fast:slow', size = 6, color = 'gray', hjust = 0) +
    annotate(geom = 'text', x = 1, y = 25, label = 'pioneer:breeder', size = 6, color = 'black', hjust = 0) +
    scale_x_log10(name = 'light') + scale_y_log10(name = 'richness ratio', limit = c(0.5, 10)) + theme_plant())

(p_ratio_light <- ggplot(richness_wide_light, aes(x = bin_midpoint)) +
    geom_point(size = 4, color = 'gray', aes(y = richness_ratio_fastslow)) +
    geom_point(size = 4, color = 'black', aes(y = richness_ratio_pioneerbreeder)) +
    annotate(geom = 'text', x = 1, y = 30, label = 'fast:slow', size = 6, color = 'gray', hjust = 0) +
    annotate(geom = 'text', x = 1, y = 25, label = 'pioneer:breeder', size = 6, color = 'black', hjust = 0) +
    scale_x_log10(name = 'light') + scale_y_log10(name = 'richness ratio') + theme_plant())


totallightbins_fg[1:10,]
bin_x_fg_light[1:10,] 

obs_totalprod[90:100,] 
bin_x_fg[90:100,] 
