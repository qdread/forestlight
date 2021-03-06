# Bins of aboveground biomass

biomassbin_alltree_1995 <- with(alltreedat[[3]], logbin_setedges(x = dbh_corr, y = agb_corr, edges = dbhbin_all))
biomassbin_byfg_1995 <- alltreedat[[3]] %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$agb_corr, edges = dbhbin_all))

biomassbin_1995 <- rbind(data.frame(fg = 'all', biomassbin_alltree_1995), ungroup(biomassbin_byfg_1995))

write.csv(biomassbin_1995, '~/google_drive/ForestLight/data/biomassbin_1995.csv', row.names = FALSE)

#### make quick plot

#Q Dog
gdrive_path <- '~/google_drive/ForestLight'
github_path <- '~/Documents/GitHub/forestlight'

#Grady
gdrive_path <- '/Users/jgradym/Google Drive/ForestLight'
github_path <- '/Users/jgradym/Documents/GitHub/forestlight'

#Grady_2
gdrive_path <- '/Users/johngrady/Google Drive/ForestLight'
github_path <- '/Users/johngrady/Documents/GitHub/forestlight'

biomassbin_1995 <- read.csv(file.path(gdrive_path,'data/biomassbin_1995.csv'), stringsAsFactors = FALSE)


fast_sum <- sum(biomassbin_1995$bin_value[biomassbin_1995$fg == "fg1"])
pioneer_sum <- sum(biomassbin_1995$bin_value[biomassbin_1995$fg == "fg2"])
slow_sum <- sum(biomassbin_1995$bin_value[biomassbin_1995$fg == "fg3"])
ratio_slow_fast <- slow_sum/fast_sum
ggplot(biomassbin_1995, aes(x = bin_midpoint, y = bin_value, color = fg)) +
  geom_point() +
  theme_bw() +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = c('black', RColorBrewer::brewer.pal(6, 'Set1')))

sums <- aggregate(biomassbin_1995$bin_value, by = list(Category = biomassbin_1995$fg), FUN = sum)
ggplot(data = sums, aes(x = Category, y = x)) + 
  geom_bar(stat = "identity") + theme_plant
  
