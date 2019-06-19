# Plot figure 5 with fits
# 19 June 2019

gdrive_path <- '~/google_drive/ForestLight'

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))
pred_area <- read.csv(file.path(gdrive_path, 'data/fig5areafit.csv'))
pred_volume <- read.csv(file.path(gdrive_path, 'data/fig5volumefit.csv'))

regdata <- with(alltree_light_95, data.frame(dbh = dbh_corr, light_area = light_received/crownarea, light_volume = light_received/crownvolume))

library(ggplot2)

# Basic plots for area and volume
# No formatting applied.

ggplot(regdata, aes(x=dbh, y=light_area)) + geom_point(alpha = 0.05) + scale_x_log10() + scale_y_log10() +
  geom_line(data = pred_area, aes(x=dbh, y=fitted), size = 1, color = 'red') +
  geom_line(data = pred_area, aes(x=dbh, y=cimin), linetype = 'dotted', color = 'red') +
  geom_line(data = pred_area, aes(x=dbh, y=cimax), linetype = 'dotted', color = 'red') 

ggplot(regdata, aes(x=dbh, y=light_volume)) + geom_point(alpha = 0.05) + scale_x_log10() + scale_y_log10() +
  geom_line(data = pred_volume, aes(x=dbh, y=fitted), size = 1, color = 'red') +
  geom_line(data = pred_volume, aes(x=dbh, y=cimin), linetype = 'dotted', color = 'red') +
  geom_line(data = pred_volume, aes(x=dbh, y=cimax), linetype = 'dotted', color = 'red') 
