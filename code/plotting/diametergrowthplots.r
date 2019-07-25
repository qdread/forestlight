# Individual diameter growth plots
#Q Dog
gdrive_path <- '~/google_drive/ForestLight'
github_path <- '~/Documents/GitHub/forestlight'

# Plot of parameters with 95% credible intervals
params <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/diamgrowth_piecewise_paramci_by_fg.csv'), stringsAsFactors = FALSE)

params %>% 
  filter(model == 2, parameter == 'beta1_high', !fg %in% 'unclassified') %>%
  ggplot(aes(x = fg, y = q50, ymin = q025, ymax = q975)) +
  geom_errorbar(width = 0.2) +
  geom_point() + 
  scale_y_continuous(limits =c(0, 0.7), expand = c(0,0)) +
  theme_bw()

# Plot with binned values and quantiles, along with fits

# bin edges
binmins <- unique(obs_dens$bin_min)
binmaxes <- unique(obs_dens$bin_max)
binbreaks <- c(binmins, binmaxes[length(binmaxes)])
binmidpoints <- exp(log(binbreaks)[-length(binbreaks)] + diff(log(binbreaks))/2)

# Load raw data
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

# We need to annualize dbh increment to a 1 year difference since the raw dbh increment is 5 year difference.

annual_increment <- function(dbh_old, dbh_new, census_interval = 5, new_interval = 1){
  rate <-  (dbh_new / dbh_old)^(1/census_interval) - 1
  dbh_oldplus1year <- dbh_old * (1 + rate)^new_interval
  return(dbh_oldplus1year - dbh_old)
}

dat <- alltreedat[[3]] %>%
  select(sp, fg, dbh_corr, dbhlastcensus) %>%
  mutate(dbh_increment = annual_increment(dbhlastcensus, dbh_corr),
         fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

# Create bins
binall <- dat %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = binbreaks, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  summarize(q025 = quantile(dbh_increment, 0.025),
            q25 = quantile(dbh_increment, 0.25),
            q50 = quantile(dbh_increment, 0.5),
            q75 = quantile(dbh_increment, 0.75),
            q975 = quantile(dbh_increment, 0.975),
            mean = exp(mean(log(dbh_increment))),
            sd = exp(sd(log(dbh_increment))),
            ci_width = qnorm(0.975) * sd(log(dbh_increment)) / sqrt(length(dbh_increment)),
            ci_min = exp(mean(log(dbh_increment)) - ci_width),
            ci_max = exp(mean(log(dbh_increment)) + ci_width),
            mean_n_individuals = length(dbh_increment))

binfg <- dat %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = binbreaks, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  summarize(q025 = quantile(dbh_increment, 0.025),
            q25 = quantile(dbh_increment, 0.25),
            q50 = quantile(dbh_increment, 0.5),
            q75 = quantile(dbh_increment, 0.75),
            q975 = quantile(dbh_increment, 0.975),
            mean = exp(mean(log(dbh_increment))),
            sd = exp(sd(log(dbh_increment))),
            ci_width = qnorm(0.975) * sd(log(dbh_increment)) / sqrt(length(dbh_increment)),
            ci_min = exp(mean(log(dbh_increment)) - ci_width),
            ci_max = exp(mean(log(dbh_increment)) + ci_width),
            mean_n_individuals = length(dbh_increment))

   

bindbhincrement <- rbind(cbind(fg = 'all', binall) %>% ungroup, binfg %>% ungroup) %>%
  mutate(bin_midpoint = binmidpoints[as.numeric(dbh_bin)],
         year = 1995
         )

ggplot(bindbhincrement %>% filter(!fg %in% 'unclassified'), aes(x = dbh, y = q50, ymin = q25, ymax = q75, color = fg)) +
  geom_pointrange() +
  theme_bw()

# Add fitted values
fittedvals <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/diamgrowth_piecewise_ci_by_fg.csv'), stringsAsFactors = FALSE) %>%
  filter(variable == 'diameter_growth_fitted') %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))


guild_fills2 <- guild_fills; guild_fills2[6] <- 'red'
ggplot(bindbhincrement %>% filter(!fg %in% 'unclassified'), aes(x = bin_midpoint, y = q50, ymin = q25, ymax = q75, color = fg)) +
  geom_pointrange() +
  geom_ribbon(data = fittedvals %>% filter(prod_model == 2, !fg %in% 'unclassified'), aes(x = dbh, ymin=q025, ymax=q975), alpha = 0.15) +
  geom_line(data = fittedvals %>% filter(prod_model == 2, !fg %in% 'unclassified'), aes(x = dbh), size = 1) +
  theme_bw() +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values=guild_fills2)

# Write files for observed and fitted.
fp <- file.path(gdrive_path, 'data/data_piecewisefits')
fp_obs <- file.path(gdrive_path, 'data/data_forplotting_aug2018')

write.csv(bindbhincrement, file = file.path(fp_obs, 'obs_indivdiamgrowth.csv'), row.names = FALSE)
write.csv(fittedvals, file = file.path(fp, 'fitted_indivdiamgrowth.csv'), row.names = FALSE)
