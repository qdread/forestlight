# Get slopes of light per crown area / diameter, and light per crown volume / diameter

library(tidyverse)
library(brms)

# Load data
gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users/jgradym/Google_Drive/ForestLight'))
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

# Do log linear regression of light received / crown area versus dbh
# Updated 20 April 2020: Use the corrected crown volume.

regdata <- alltree_light_95 %>%
  select(dbh_corr, light_received_byarea, light_received_byvolume) %>%
  setNames(c('dbh', 'light_area', 'light_volume'))

# Plot for sanity check
#ggplot(regdata, aes(x=dbh, y=light_volume)) + geom_point(alpha = 0.05) + scale_x_log10() + scale_y_log10()

# fit model.
options(mc.cores = 3)
set.seed(111)

reg_area <- brm(log10(light_area) ~ log10(dbh), data = regdata, chains = 3, iter = 2000, warmup = 1000, save_warmup = FALSE)
reg_volume <- brm(log10(light_volume) ~ log10(dbh), data = regdata, chains = 3, iter = 2000, warmup = 1000, save_warmup = FALSE)

save(reg_area, reg_volume, file = '~/Documents/temp/fig1brmsfits.RData')

summ_area <- summary(reg_area)
summ_volume <- summary(reg_volume)

fixef(reg_area)
fixef(reg_volume)

brms::bayes_R2(reg_area)
brms::bayes_R2(reg_volume)

# extract samples to calc Bayes R2
samples_area <- map_dfr(reg_area$fit@sim$samples, ~ data.frame(intercept = .$b_Intercept, slope = .$b_log10dbh))

# Use median parameter estimate to get Bayesian R2.
area_log_fitted <- cbind(1, log10(regdata$dbh)) %*% fixef(reg_area)[,'Estimate'] # Fitted values
area_residuals <- log10(regdata$light_area) - area_log_fitted
area_var_fitted <- var(area_log_fitted)
area_var_resid <- var(area_residuals)
area_rsq <- area_var_fitted / (area_var_fitted + area_var_resid)

volume_log_fitted <- cbind(1, log10(regdata$dbh)) %*% fixef(reg_volume)[,'Estimate'] # Fitted values
volume_residuals <- log10(regdata$light_volume) - volume_log_fitted
volume_var_fitted <- var(volume_log_fitted)
volume_var_resid <- var(volume_residuals)
volume_rsq <- volume_var_fitted / (volume_var_fitted + volume_var_resid)

# Area 0.680, Volume 0.262


# save summary statistics
coef_table <-rbind(
  bind_rows(data.frame(regression = 'light per unit crown area versus dbh', parameter = c('intercept', 'slope'), fixef(reg_area)),
            data.frame(regression = 'light per unit crown area versus dbh', parameter = c('r-squared'), Estimate = area_rsq)),
  bind_rows(data.frame(regression = 'light per unit crown volume versus dbh', parameter = c('intercept', 'slope'), fixef(reg_volume)),
            data.frame(regression = 'light per unit crown volume versus dbh', parameter = c('r-squared'), Estimate = volume_rsq)))

write_csv(coef_table, file.path(gdrive_path, 'data/clean_summary_tables/fig1_light_by_size_parameters.csv'))

###############################################

# Plot the results
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

#write.csv(pred_area, '~/google_drive/ForestLight/data/fig5areafit.csv', row.names = FALSE)
#write.csv(pred_volume, '~/google_drive/ForestLight/data/fig5volumefit.csv', row.names = FALSE)
