# Plots of diameter and height increment by diameter for all FGs.


# Load raw data -----------------------------------------------------------


gdrive_path <- '~/google_drive/ForestLight'
github_path <- '~/Documents/GitHub/forestlight'

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))
source(file.path(github_path, 'code/allfunctions27july.r'))
load(file.path(gdrive_path, 'data/BCI_raw/bcidata/bciqcrun.R'))


# Calculate growth rates --------------------------------------------------


# We need to "annualize" the 5-year increments in height and in diameter to 1-year rates.

annual_increment <- function(x_old, x_new, census_interval = 5, new_interval = 1){
  rate <-  (x_new / x_old)^(1/census_interval) - 1
  x_oldplus1year <- x_old * (1 + rate)^new_interval
  return(x_oldplus1year - x_old)
}

annual_height_growth <- annual_increment(bci.full3$height_corr, bci.full4$height_corr)
annual_diameter_growth <- annual_increment(bci.full3$dbh_corr, bci.full4$dbh_corr)
annual_biomass_growth <- annual_increment(bci.full3$agb_corr, bci.full4$agb_corr)

annual_growth_rates <- data.frame(treeID = bci.full4$treeID,
                                  height_growth = annual_height_growth,
                                  diameter_growth = annual_diameter_growth,
                                  biomass_growth = annual_biomass_growth)

library(tidyverse)

# Convert the units to proper units
# height stays in m, diameter in cm, biomass in kg
growthdat95 <- alltreedat[[3]] %>% 
  left_join(annual_growth_rates) %>%
  mutate(diameter_growth = diameter_growth / 10,
         biomass_growth = biomass_growth * 1000)


# Get bin values for plotting ---------------------------------------------

# Calculate the DBH bin bounds (same as the ones used for all other dbh analysis)
n <- 20
dbh_range <- log(range(unlist(map(alltreedat, 'dbh_corr'))))
log_bin_edges <- seq(dbh_range[1], dbh_range[2], length.out = n + 1)
bin_edges <- exp(log_bin_edges) # get edges of bins
bin_midpoints <- exp(log_bin_edges[-length(log_bin_edges)] + diff(log_bin_edges)/2)

# Get quantiles from each bin
growthdat95_binned <- growthdat95 %>%
  select(fg, dbh_corr, height_growth, diameter_growth, biomass_growth) %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = bin_edges, include.lowest = TRUE)) %>%
  gather(variable, growth, -fg, -dbh_corr, -dbh_bin) %>%
  group_by(dbh_bin, variable) %>%
  summarize(median = quantile(growth, 0.5),
            minimum = quantile(growth, 0.025),
            maximum = quantile(growth, 0.975)) %>%
  mutate(bin_midpoint = bin_midpoints[as.numeric(dbh_bin)])

# Create plots ------------------------------------------------------------

# Plots with medians, 2.5%, and 97.5%
# Using the same 20 bins as all other bin plots

pbin1 <- ggplot(growthdat95_binned %>% filter(variable == 'biomass_growth'), aes(x = bin_midpoint, y = median, ymin = minimum, ymax = maximum)) +
  geom_errorbar(width = 0) +
  geom_point(pch = 21) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = 'Biomass growth (kg/y)') +
  theme_bw()

pbin2 <- ggplot(growthdat95_binned %>% filter(variable == 'diameter_growth'), aes(x = bin_midpoint, y = median, ymin = minimum, ymax = maximum)) +
  geom_errorbar(width = 0) +
  geom_point(pch = 21) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = 'Diameter growth (cm/y)') +
  theme_bw()

pbin3 <- ggplot(growthdat95_binned %>% filter(variable == 'height_growth'), aes(x = bin_midpoint, y = median, ymin = minimum, ymax = maximum)) +
  geom_errorbar(width = 0) +
  geom_point(pch = 21) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = 'Height growth (m/y)') +
  theme_bw()

# Hexagon plots
# Using same hex color scale as other hex plots
hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000), 
                                             labels = c(1,10,100,1000), limits=c(1,5000))

phex1 <- ggplot(growthdat95, aes(x = dbh_corr, y = biomass_growth)) +
  geom_hex() +
  hex_scale_log_colors +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = 'Biomass growth (kg/y)') +
  theme_bw()

phex2 <- ggplot(growthdat95, aes(x = dbh_corr, y = diameter_growth)) +
  geom_hex() +
  hex_scale_log_colors +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = 'Diameter growth (cm/y)', limits = c(0.01, 4)) +
  theme_bw()

phex3 <- ggplot(growthdat95, aes(x = dbh_corr, y = height_growth)) +
  geom_hex() +
  hex_scale_log_colors +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = 'Height growth (m/y)', limits = c(0.001, 2)) +
  theme_bw()

# Write to pdf for now
pdf(file.path(gdrive_path, 'figs/height_diameter_biomass_increments.pdf'), height = 6, width = 6)
pbin1; pbin2; pbin3; phex1; phex2; phex3
dev.off()


# Plots with bins, by FG --------------------------------------------------

# Get quantiles from each bin for each fg
fg_growthdat95_binned <- growthdat95 %>%
  select(fg, dbh_corr, height_growth, diameter_growth, biomass_growth) %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = bin_edges, include.lowest = TRUE),
         fg = factor(fg)) %>%
  gather(variable, growth, -fg, -dbh_corr, -dbh_bin) %>%
  group_by(fg, dbh_bin, variable) %>%
  summarize(n = n(),
            median = quantile(growth, 0.5),
            minimum = quantile(growth, 0.025),
            maximum = quantile(growth, 0.975)) %>%
  filter(n >= 3) %>%
  mutate(bin_midpoint = bin_midpoints[as.numeric(dbh_bin)])

cols <- c("#BFE046", "#27408b", "#267038", "#87Cefa", "ivory")
dodge_width <- 0.03

fgpbin1 <- ggplot(fg_growthdat95_binned %>% filter(variable == 'biomass_growth', !is.na(fg)), aes(x = bin_midpoint, y = median, ymin = minimum, ymax = maximum, group = fg, fill = fg)) +
  geom_errorbar(width = 0, color = 'black', position = position_dodge(width = dodge_width)) +
  geom_point(pch = 21, position = position_dodge(width = dodge_width)) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = 'Biomass growth (kg/y)') +
  scale_fill_manual(values = cols) +
  theme_bw()

fgpbin2 <- ggplot(fg_growthdat95_binned %>% filter(variable == 'diameter_growth', !is.na(fg)), aes(x = bin_midpoint, y = median, ymin = minimum, ymax = maximum, group = fg, fill = fg)) +
  geom_errorbar(width = 0, color = 'black', position = position_dodge(width = dodge_width)) +
  geom_point(pch = 21, position = position_dodge(width = dodge_width)) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = 'Diameter growth (cm/y)') +
  scale_fill_manual(values = cols) +
  theme_bw()

fgpbin3 <- ggplot(fg_growthdat95_binned %>% filter(variable == 'height_growth', !is.na(fg)), aes(x = bin_midpoint, y = median, ymin = minimum, ymax = maximum, group = fg, fill = fg)) +
  geom_errorbar(width = 0, color = 'black', position = position_dodge(width = dodge_width)) +
  geom_point(pch = 21, position = position_dodge(width = dodge_width)) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = 'Height growth (m/y)') +
  scale_fill_manual(values = cols) +
  theme_bw()
