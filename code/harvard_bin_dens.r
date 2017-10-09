# Harvard census
hfstems <- read.csv('C:/Users/Q/google_drive/ForestLight/HF/hf253-04-stems-2014.csv', stringsAsFactors = FALSE)

length(unique(hfstems$tree.id))

# Source code for log binning function
source('code/allfunctions27july.r')

# Set number of bins       
numbins <- 20

dbhbin_hf <- logbin(x = hfstems$dbh, y = NULL, n = numbins)

area_core <- 35

library(dplyr)
library(cowplot)

# Figure 3b
dbhbin_hf %>%
  mutate(bin_value = bin_value/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_value)) +
  geom_point(size = 2) +
  geom_abline(intercept = 4, slope = -2) +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  panel_border(colour = 'black')

# Export Harvard Forest
dbhbin_hf$value_per_ha <- dbhbin_hf$bin_value/35
write.csv(dbhbin_hf, 'C:/Users/Q/google_drive/ForestLight/data/data_04oct/harvard_diameter_bin.csv', row.names = FALSE)
