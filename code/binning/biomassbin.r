# Bins of aboveground biomass

biomassbin_alltree_1995 <- with(alltreedat[[3]], logbin_setedges(x = dbh_corr, y = agb_corr, edges = dbhbin_all))
biomassbin_byfg_1995 <- alltreedat[[3]] %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$agb_corr, edges = dbhbin_all))

biomassbin_1995 <- rbind(data.frame(fg = 'all', biomassbin_alltree_1995), ungroup(biomassbin_byfg_1995))

write.csv(biomassbin_1995, '~/google_drive/ForestLight/data/biomassbin_1995.csv', row.names = FALSE)

#### make quick plot

biomassbin_1995 <- read.csv('~/google_drive/ForestLight/data/biomassbin_1995.csv', stringsAsFactors = FALSE)

ggplot(biomassbin_1995, aes(x = bin_midpoint, y = bin_value, color = fg)) +
  geom_point() +
  theme_bw() +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = c('black', RColorBrewer::brewer.pal(6, 'Set1')))
