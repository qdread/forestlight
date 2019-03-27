gdrive_path <- '~/google_drive'
github_path <- '~/Documents/GitHub/forestlight'

load(file.path(gdrive_path, 'ForestLight/data/rawdataobj_alternativecluster.r'))
source(file.path(github_path, 'code/allfunctions27july.r'))

library(dplyr)
library(ggplot2)

alltree_light_95 <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), as.character(NA)))

# Create bins
# -----------
area_core <- 42.84
num_bins <- 20

allyeardbh <- unlist(lapply(alltreedat[2:6], '[', , 'dbh_corr'))
dbhbin_all <- logbin(x = allyeardbh, y = NULL, n = numbins)

totallightpervolumebins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = light_received_byvolume, edges = dbhbin_all))
totallightpervolumebins_fg <- alltree_light_95 %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$light_received_byvolume, edges = dbhbin_all))

totallightpervolumebins_fg <- rbind(data.frame(fg = 'all', totallightpervolumebins_all, stringsAsFactors = FALSE), as.data.frame(totallightpervolumebins_fg)) %>%
  mutate(bin_value = bin_value / area_core)

# Make plot
# ---------

mycols <- c('black', RColorBrewer::brewer.pal(5,'Set1'))

plightbyvol_onefig <- ggplot(totallightpervolumebins_fg %>% filter(bin_count > 0)) +
  geom_point(aes(x = bin_midpoint, y = bin_value, color = fg)) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = expression(paste('Total incoming light (relative units)', sep=''))) +
  theme_bw() +
  scale_color_manual(values=mycols) + scale_fill_manual(values=mycols)
