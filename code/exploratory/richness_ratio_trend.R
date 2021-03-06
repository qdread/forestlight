# Exploratory script: richness of each FG, and ratios, vs. diameter
# QDR / forestlight / 9 June 2020

library(tidyverse)
library(forestscaling)

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users/jgradym/Google Drive/ForestLight'))
github_path <- ifelse(Sys.info()['user'] == 'qread', '~/documents/github/MSU_repos', file.path('/Users/jgradym/Documents/Github'))

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

# 1995 data with light values (113,651 trees)
dat_light <- dat %>%
  filter(!is.na(light_received_byarea))

bin_x_fg <- bin_x_fg %>%
  cbind(pmap_dfr(bin_x_fg, function(fg, bin_min, bin_max, ...) {
    sp_ids <- as.character(dat$sp[dat$fg %in% fg & dat$dbh_corr >= bin_min & dat$dbh_corr < bin_max])
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min))

bin_x_fg_light <- bin_x_fg_light %>%
  cbind(pmap_dfr(bin_x_fg_light, function(fg, bin_min, bin_max, ...) {
    sp_ids <- as.character(dat_light$sp[dat_light$fg %in% fg & dat_light$light_received_byarea >= bin_min & dat_light$light_received_byarea < bin_max])
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min))

# Quick ggplot ------------------------------------------------------------

guild_colors <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")

## richness for each group

# Arith y scale
(p_rich <- ggplot(bin_x_fg %>% 
                   # arrange(desc(fg)) %>% # what !
                    filter(!fg %in% 'unclassified' & richness > 0), aes(x = bin_midpoint, y = richness_by_bin_width, color = fg)) + 
  geom_point(size = 4) + scale_x_log10(name = 'diameter') + theme_plant() +
  scale_color_manual(values = guild_colors) +
  theme(legend.position = 'right'))

# Log y scale
(p_rich_log <- p_rich + scale_y_log10()) #limits = c(0.6, 140)))

## ratios of richness

richness_wide <- bin_x_fg %>%
  pivot_wider(id_cols = c(bin, bin_midpoint, bin_min, bin_max), names_from = fg, values_from = c(richness, n_individuals)) %>%
  mutate(richness_ratio_fastslow = richness_fg1/richness_fg3,
         richness_ratio_pioneerbreeder = richness_fg2/richness_fg4,
         lowest_n_fastslow = pmin(n_individuals_fg1, n_individuals_fg3),
         lowest_rich_fastslow  = pmin(richness_fg1, richness_fg3),
         lowest_n_pioneerbreeder = pmin(n_individuals_fg2, n_individuals_fg4),
         lowest_rich_pioneerbreeder = pmin(richness_fg2, richness_fg4)) %>%
  mutate_if(is.double, ~ if_else(is.finite(.), ., as.numeric(NA)))

# Richness ratio 
(p_ratio <- ggplot(richness_wide, aes(x = bin_midpoint)) +
  geom_point(size = 4, color = 'gray', aes(y = richness_ratio_fastslow)) +
  geom_point(size = 4, color = 'black', aes(y = richness_ratio_pioneerbreeder)) +
  annotate(geom = 'text', x = 1, y = 30, label = 'fast:slow', size = 6, color = 'gray', hjust = 0) +
  annotate(geom = 'text', x = 1, y = 25, label = 'pioneer:breeder', size = 6, color = 'black', hjust = 0) +
  scale_x_log10(name = 'diameter') + scale_y_log10(name = 'richness ratio', limit = c(0.5, 10)) + theme_plant())

dfs <- richness_wide %>% filter(lowest_n_fastslow >= 20)
lm1 <- lm(log(richness_ratio_fastslow) ~ log(bin_midpoint), data= dfs)
summary(lm1) # 0.23 (0.18, 0.28), r2 = 0.88
confint(lm1) 
dst <- richness_wide %>% filter(lowest_n_pioneerbreeder >= 20)
lm2 <- lm(log(richness_ratio_pioneerbreeder) ~ log(bin_midpoint), data=dst)
summary(lm2) #0.80 (0.68, 0.91), r2 = 0.97
confint(lm2)

(p_ratio <- ggplot(richness_wide, aes(x = bin_midpoint)) +
    geom_point(size = 4, color = 'gray', aes(y = richness_ratio_fastslow)) +
    geom_point(size = 4, color = 'black', aes(y = richness_ratio_pioneerbreeder)) +
    annotate(geom = 'text', x = 1, y = 30, label = 'fast:slow', size = 6, color = 'gray', hjust = 0) +
    annotate(geom = 'text', x = 1, y = 25, label = 'pioneer:breeder', size = 6, color = 'black', hjust = 0) +
    scale_x_log10(name = 'diameter') + scale_y_log10(name = 'richness ratio', limit = c(0.5, 2)) + theme_plant())

#----- by light------------------------

(p_rich_light <- ggplot(bin_x_fg_light %>% 
                          arrange(desc(fg)) %>% # ok
                          filter(!fg %in% 'unclassified' & richness > 0), 
                        aes(x = bin_midpoint, y = richness, color = fg)) + 
   geom_point(size = 4) + scale_x_log10(name = 'light per unit crown area') + theme_plant() +
   scale_color_manual(values = guild_colors) +
   scale_y_log10() + #limits = c(0.6, 140)
   theme(legend.position = 'right'))


richness_wide_light <- bin_x_fg_light %>%
  pivot_wider(id_cols = c(bin, bin_midpoint, bin_min, bin_max), names_from = fg, values_from = c(richness, n_individuals)) %>%
  mutate(richness_ratio_fastslow = richness_fg1/richness_fg3,
         richness_ratio_pioneerbreeder = richness_fg2/richness_fg4,
         lowest_n_fastslow = pmin(n_individuals_fg1, n_individuals_fg3),
         lowest_rich_fastslow  = pmin(richness_fg1, richness_fg3),
         lowest_n_pioneerbreeder = pmin(n_individuals_fg2, n_individuals_fg4),
         lowest_rich_pioneerbreeder = pmin(richness_fg2, richness_fg4)) %>%
  mutate_if(is.double, ~ if_else(is.finite(.), ., as.numeric(NA)))

lfs <- richness_wide_light %>% filter(lowest_n_fastslow >= 20)
lm3 <- lm(log(richness_ratio_fastslow) ~ log(bin_midpoint), data= lfs)
summary(lm3) # 0.25 (0.17 - 0.32) r2 = 0.96
confint(lm3)
lst <- richness_wide_light %>% filter(lowest_n_pioneerbreeder >= 20)
lm4 <- lm(log(richness_ratio_pioneerbreeder) ~ log(bin_midpoint), data=lst)
summary(lm4) #0.80, r2 = 0.96
confint(lm4) #0.69, 0.91


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
