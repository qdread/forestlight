# Get slopes of light per crown area / diameter, and light per crown volume / diameter

library(ggplot2)
library(brms)

# Load data (changing file path if necessary)
load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')

# Do log linear regression of light received / crown area versus dbh

regdata <- with(alltree_light_95, data.frame(dbh = dbh_corr, light_area = light_received/crownarea, light_volume = light_received/crownvolume))

# Plot for sanity check
#ggplot(regdata, aes(x=dbh, y=light_volume)) + geom_point(alpha = 0.05) + scale_x_log10() + scale_y_log10()

# fit model.
options(mc.cores = 3)

reg_area <- brm(log10(light_area) ~ log10(dbh), data = regdata, chains = 3, iter = 2000, warmup = 1000)
reg_volume <- brm(log10(light_volume) ~ log10(dbh), data = regdata, chains = 3, iter = 2000, warmup = 1000)

summ_area <- summary(reg_area)
summ_volume <- summary(reg_volume)

fixef(reg_area)
fixef(reg_volume)

bayes_R2(reg_area)
bayes_R2(reg_volume)

dbh_rad_slope <- 0.658
dbh_vol_slope <- 2.02

dbh_vol_slope - 2*dbh_rad_slope # Should be exactly equal to the difference between the two regression slopes. Works out.

# Save the fitted values, with intervals, to a CSV to load and plot.
dbh_pred <- exp(seq(log(1.1), log(315), length.out = 101))
pred_area <- fitted(reg_area, newdata = data.frame(dbh = dbh_pred), summary = TRUE)
pred_area <- data.frame(dbh = dbh_pred, fitted = 10^pred_area[,'Estimate'], cimin = 10^pred_area[,'Q2.5'], cimax = 10^pred_area[,'Q97.5'])

ggplot(regdata, aes(x=dbh, y=light_area)) + geom_point(alpha = 0.05) + scale_x_log10() + scale_y_log10() +
  geom_line(data = pred_area, aes(x=dbh, y=fitted), size = 1, color = 'red') +
  geom_line(data = pred_area, aes(x=dbh, y=cimin), linetype = 'dotted', color = 'red') +
  geom_line(data = pred_area, aes(x=dbh, y=cimax), linetype = 'dotted', color = 'red') 


pred_volume <- fitted(reg_volume, newdata = data.frame(dbh = dbh_pred), summary = TRUE)
pred_volume <- data.frame(dbh = dbh_pred, fitted = 10^pred_volume[,'Estimate'], cimin = 10^pred_volume[,'Q2.5'], cimax = 10^pred_volume[,'Q97.5'])

ggplot(regdata, aes(x=dbh, y=light_volume)) + geom_point(alpha = 0.05) + scale_x_log10() + scale_y_log10() +
  geom_line(data = pred_volume, aes(x=dbh, y=fitted), size = 1, color = 'red') +
  geom_line(data = pred_volume, aes(x=dbh, y=cimin), linetype = 'dotted', color = 'red') +
  geom_line(data = pred_volume, aes(x=dbh, y=cimax), linetype = 'dotted', color = 'red') 

write.csv(pred_area, '~/google_drive/ForestLight/data/fig5areafit.csv', row.names = FALSE)
write.csv(pred_volume, '~/google_drive/ForestLight/data/fig5volumefit.csv', row.names = FALSE)
