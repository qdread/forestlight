# Bin and plot the following:
# 1. light per volume (uncorrected) ~ diameter
# 2. light per volume, corrected with several different k values ~ diameter


# Load data ---------------------------------------------------------------

library(tidyverse)

user <- Sys.info()['user']
gdrive_path <- ifelse(user == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))
github_path <- ifelse(user == 'qread', '~/Documents/GitHub/forestlight', file.path('/Users',user,'Documents/GitHub/forestlight'))

source(file.path(github_path, 'code/allfunctions27july.r'))
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

alltree_light_95 <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), as.character(NA)))

dbhbin1995 <- with(alltree_light_95, logbin(x = dbh_corr, n = 20))

# Correct incident light capture percentage, using crown depth in meters, for a range of k values
# Crown depth varies from ~1-38

pct_light_captured <- function(depth, k) 1 - exp(-k * depth)


# Plot analytically with a series of k ------------------------------------

depth_values <- seq(0, 40, length.out = 101)
k_values <- seq(0.3, 0.7, by = 0.05)

cross_values <- cross2(depth_values, k_values)

plot_dat <- data.frame(depth = map_dbl(cross_values, 1),
                       k = map_dbl(cross_values, 2),
                       pct_captured = map_dbl(cross_values, ~ pct_light_captured(depth = .[[1]], k = .[[2]])))

ggplot(plot_dat, aes(x = depth, y = pct_captured, group = k, color = k)) +
  geom_line(size = 1) +
  theme_minimal() +
  scale_color_viridis_c(breaks = seq(0.3, 0.7, by = 0.1)) +
  scale_y_continuous(name = 'proportion light captured') +
  scale_x_continuous(name = 'crown depth (m)', limits = c(0, 20))


# Calculate percent light capture for all trees given different k's -------

k_values <- c(seq(0.3, 0.7, by = 0.1), Inf)

treedat95 <- alltree_light_95 %>% 
  select(dbh_corr, light, crownarea, crownvolume, crowndepth, light_received_byarea, light_received_byvolume)

# Correction factors for each k
corr_factors <- map(k_values, ~ pct_light_captured(depth = treedat95$crowndepth, k = .))

corrected_light_by_volume <- map(corr_factors, ~ . * treedat95$light_received_byvolume)


# Bin the data for visualization ------------------------------------------

# Combine corrected values into a data frame
treedat_withcorr <- map2_dfr(k_values, corrected_light_by_volume, ~ data.frame(dbh_corr = treedat95$dbh_corr, 
                                                                               k = .x,
                                                                               corrected_light_received_by_volume = .y))

# Do bins for each of the groups
bin_breaks <- c(dbhbin1995$bin_min[1], dbhbin1995$bin_max[-length(dbhbin1995$bin_max)], dbhbin1995$bin_max[length(dbhbin1995$bin_max)] + 1)

correctedlightvol_fakebin <- treedat_withcorr %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = bin_breaks, labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(k, dbh_bin) %>%
  do(quantile(.$corrected_light_received_by_volume, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975'))) %>%
  ungroup %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))



# Create plots ------------------------------------------------------------

ggplot(correctedlightvol_fakebin %>% mutate(k = factor(k)), aes(x = dbh_bin, y = q50, ymin = q025, ymax = q975, group = k, color = k, fill = k)) +
  geom_line(position = position_dodge(width = 0.04)) +
  geom_point(position = position_dodge(width = 0.04)) +
  geom_errorbar(position = position_dodge(width = 0.04)) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = parse(text = 'Light~captured~per~volume~(W~m^-3)')) +
  theme_minimal() +
  ggtitle('Light capture per volume for different values of k', 'k = Inf is the uncorrected value we used before')


# Slopes for the different k values ---------------------------------------

treedat_withcorr %>%
  group_by(k) %>%
  summarize(scaling_coefficient = lm(log10(corrected_light_received_by_volume) ~ log10(dbh_corr))$coeff[2])
