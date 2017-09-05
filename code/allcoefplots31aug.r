# Coefficient plots and tables, 31 Aug
# Results of new analysis

# List of coefficient plots/tables to generate:

# 1. density scaling by diameter, with and without cutoff. AICs for each, as well as confidence intervals (allyears, 90light, and 95light)
# 2. density scaling by light, with and without cutoff (AICs, CIs) (90, 95)
# 3. individual production by diameter, with and without light received used as a covariate; AIC and CI for each (90, 95)
# 4. individual production by light received (90 and 95)
# 5. individual light received by diameter
# 6. binned production by diameter (all years, 90, and 95)
# 7. binned light received by diameter (90 and 95)

# 1. Get density scaling by diameter coefficients, with and without cutoff, for each year.

library(dplyr)

guilds <- c('all','shade','intermediate','gap','unclassified','shade_quantile','intermediate_quantile','gap_quantile')
allyear_dens_coeff <- list()

for (n in 1:length(allyears_names)) {
  i <- allyears_names[n]
  dat_coeffs <- get(paste0(i, '_paretofits'))
  dat_ci <- get(paste0(i, '_paretobootci'))
  dat_coeffs <- lapply(dat_coeffs, function(x) c(x$fit_pareto@coef, x$fit_cutoff@coef))
  dat_ci <- lapply(dat_ci, t)
  dat_ci <- do.call('rbind', dat_ci)
  
  dat <- cbind(data.frame(year = rep(years,each=3),  model = c('pareto','pareto_cutoff','pareto_cutoff'), coefficient = c('alpha','alpha','L')), value = unlist(dat_coeffs), ci_min = dat_ci[,1], ci_max = dat_ci[,2])
  
  dat <- cbind(guild = guilds[n], dat)
  allyear_dens_coeff[[n]] <- dat
}

allyear_dens_coeff <- do.call('rbind', allyear_dens_coeff)

# Create plot.

allyear_dens_coeff %>%
  filter(guild %in% c('shade','intermediate','gap'), coefficient == 'alpha') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max, colour = model, group = model)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling slopes', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'coefficients_densityscaling_evenshadegroups.png'), height = 6, width = 9, dpi = 300)

allyear_dens_coeff %>%
  filter(guild %in% c('shade','intermediate','gap'), coefficient == 'L') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling cutoffs', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'coefficients_densityscalingcutoffs_evenshadegroups.png'), height = 6, width = 9, dpi = 300)

allyear_dens_coeff %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile'), coefficient == 'alpha') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max, colour = model, group = model)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling slopes', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'coefficients_densityscaling_quantileshadegroups.png'), height = 6, width = 9, dpi = 300)

allyear_dens_coeff %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile'), coefficient == 'L') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling cutoffs', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'coefficients_densityscalingcutoffs_quantileshadegroups.png'), height = 6, width = 9, dpi = 300)

# AIC scores for the coefficients in number 1.
allyear_deltaaic <- lapply(allyears_names, function(x) {
  dat <- get(paste0(x, '_paretofits'))
  lapply(dat, function(z) with(z, aic_pareto - aic_cutoff))
})

allyear_deltaaic <- data.frame(guild = rep(guilds,each=6), year=years, deltaaicc = unlist(allyear_deltaaic))

allyear_deltaaic %>%
  filter(guild %in% c('shade','intermediate','gap')) %>%
  ggplot(aes(x = guild, y = deltaaicc)) +
  facet_wrap(~ year) +
  geom_point() + 
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('AIC improvement from including cutoff', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'deltaaicc_evenshadegroups.png'), height=6, width=9, dpi=300)

allyear_deltaaic %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile')) %>%
  ggplot(aes(x = guild, y = deltaaicc)) +
  facet_wrap(~ year) +
  geom_point() + 
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('AIC improvement from including cutoff', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'deltaaicc_quantileshadegroups.png'), height=6, width=9, dpi=300)

# 2. Density by light received coefficients

dens_light_coeff1990 <- list()

for (n in 1:length(names1990)) {
  i <- names1990[n]
  dat_coeffs <- get(paste0(i, '_lightparetofits'))
  dat_ci <- get(paste0(i, '_lightparetobootci'))
  dat_coeffs <- c(dat_coeffs$fit_pareto@coef, dat_coeffs$fit_cutoff@coef)
  dat_ci <- t(dat_ci)

  dat <- cbind(data.frame(year = 1990,  model = c('pareto','pareto_cutoff','pareto_cutoff'), coefficient = c('alpha','alpha','L')), value = dat_coeffs, ci_min = dat_ci[,1], ci_max = dat_ci[,2])
  
  dat <- cbind(guild = guilds[n], dat)
  dens_light_coeff1990[[n]] <- dat
}

dens_light_coeff1995 <- list()

for (n in 1:length(names1995)) {
  i <- names1995[n]
  dat_coeffs <- get(paste0(i, '_lightparetofits'))
  dat_ci <- get(paste0(i, '_lightparetobootci'))
  dat_coeffs <- c(dat_coeffs$fit_pareto@coef, dat_coeffs$fit_cutoff@coef)
  dat_ci <- t(dat_ci)
  
  dat <- cbind(data.frame(year = 1995,  model = c('pareto','pareto_cutoff','pareto_cutoff'), coefficient = c('alpha','alpha','L')), value = dat_coeffs, ci_min = dat_ci[,1], ci_max = dat_ci[,2])
  
  dat <- cbind(guild = guilds[n], dat)
  dens_light_coeff1995[[n]] <- dat
}

dens_light_coeff <- do.call('rbind', c(dens_light_coeff1990, dens_light_coeff1995))

dens_light_coeff %>%
  filter(guild %in% c('shade','intermediate','gap'), coefficient == 'alpha') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max, colour = model, group = model)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling slopes by light received', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'coefficients_densitybylightscaling_evenshadegroups.png'), height = 6, width = 8, dpi = 300)

dens_light_coeff %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile'), coefficient == 'alpha') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max, colour = model, group = model)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling slopes by light received', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'coefficients_densitybylightscaling_quantileshadegroups.png'), height = 6, width = 9, dpi = 300)


dens_light_coeff %>%
  filter(guild %in% c('shade','intermediate','gap'), coefficient == 'L') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling cutoffs by light received', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'coefficients_densitybylightscalingcutoffs_evenshadegroups.png'), height = 6, width = 8, dpi = 300)

dens_light_coeff %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile'), coefficient == 'L') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling cutoffs by light received', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'coefficients_densitybylightscalingcutoffs_quantileshadegroups.png'), height = 6, width = 8, dpi = 300)

# Delta aic for the plots in 2.
delta1990 <- lapply(names1990, function(x) with(get(paste0(x,'_lightparetofits')), aic_pareto - aic_cutoff))
delta1995 <- lapply(names1990, function(x) with(get(paste0(x,'_lightparetofits')), aic_pareto - aic_cutoff))

deltaaiclight <- data.frame(guild = guilds, year = rep(c(1990,1995),each=8), deltaaicc = unlist(c(delta1990,delta1995)))

deltaaiclight %>%
  filter(guild %in% c('shade','intermediate','gap')) %>%
  ggplot(aes(x = guild, y = deltaaicc)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'blue') +
  facet_wrap(~ year) +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('AICc improvement from adding cutoff to light-based scaling', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'deltaaic_lightscaling_evenshadegroups.png'), height=5, width=8, dpi=300)
  
deltaaiclight %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile')) %>%
  ggplot(aes(x = guild, y = deltaaicc)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'blue') +
  facet_wrap(~ year) +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('AICc improvement from adding cutoff to light-based scaling', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'deltaaic_lightscaling_quantileshadegroups.png'), height=5, width=8, dpi=300)

# 3. Individual production by diameter coefficients

prod_slopes_1990 <- lapply(names1990, function(x) get(paste0(x, '_dens_slopes')))
for (i in 1:length(guilds)) prod_slopes_1990[[i]] <- data.frame(year = 1990, guild = guilds[i], prod_slopes_1990[[i]])
prod_slopes_1995 <- lapply(names1995, function(x) get(paste0(x, '_dens_slopes')))
for (i in 1:length(guilds)) prod_slopes_1995[[i]] <- data.frame(year = 1995, guild = guilds[i], prod_slopes_1995[[i]])

prod_slopes <- do.call('rbind', c(prod_slopes_1990, prod_slopes_1995))

prod_slopes %>%
  filter(guild %in% c('shade','intermediate','gap')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax, group = withindex, colour = withindex)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Production slopes with and without light', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'coefficients_individualproductionscaling_evenshadegroups.png'), height = 6, width = 8, dpi = 300)

prod_slopes %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax, group = withindex, colour = withindex)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Production slopes with and without light', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'coefficients_individualproductionscaling_quantileshadegroups.png'), height = 6, width = 8, dpi = 300)

# 4. Individual production by light received

prod_light_slopes_1990 <- lapply(names1990, function(x) {
  dat <- get(x)
  lm_dat <- lm(I(log10(production)) ~ I(log10(light_received)), data = dat)
  c(slope = as.numeric(lm_dat$coef[2]), cimin = confint(lm_dat)[2,1], cimax = confint(lm_dat)[2,2])
})

prod_light_slopes_1990 <- data.frame(year = 1990, guild = guilds, do.call('rbind', prod_light_slopes_1990))

prod_light_slopes_1995 <- lapply(names1995, function(x) {
  dat <- get(x)
  lm_dat <- lm(I(log10(production)) ~ I(log10(light_received)), data = dat)
  c(slope = as.numeric(lm_dat$coef[2]), cimin = confint(lm_dat)[2,1], cimax = confint(lm_dat)[2,2])
})

prod_light_slopes_1995 <- data.frame(year = 1995, guild = guilds, do.call('rbind', prod_light_slopes_1995))

prod_light_slopes <- rbind(prod_light_slopes_1990, prod_light_slopes_1995)

prod_light_slopes %>%
  filter(guild %in% c('shade','intermediate','gap')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Production slopes as a function of light received', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'coefficients_individualproductionbylightscaling_evenshadegroups.png'), height = 6, width = 8, dpi = 300)


prod_light_slopes %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Production slopes as a function of light received', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'coefficients_individualproductionbylightscaling_quantileshadegroups.png'), height = 6, width = 8, dpi = 300)


# 5. Individual light received by diameter

light_diam_slopes_1990 <- lapply(names1990, function(x) {
  dat <- get(x)
  lm_dat <- lm(I(log10(light_received)) ~ I(log10(dbh_corr)), data = dat)
  c(slope = as.numeric(lm_dat$coef[2]), cimin = confint(lm_dat)[2,1], cimax = confint(lm_dat)[2,2])
})
light_diam_slopes_1990 <- data.frame(year = 1990, guild = guilds, do.call('rbind', light_diam_slopes_1990))

light_diam_slopes_1995 <- lapply(names1995, function(x) {
  dat <- get(x)
  lm_dat <- lm(I(log10(light_received)) ~ I(log10(dbh_corr)), data = dat)
  c(slope = as.numeric(lm_dat$coef[2]), cimin = confint(lm_dat)[2,1], cimax = confint(lm_dat)[2,2])
})
light_diam_slopes_1995 <- data.frame(year = 1995, guild = guilds, do.call('rbind', light_diam_slopes_1995))


light_diam_slopes <- rbind(light_diam_slopes_1990, light_diam_slopes_1995)

light_diam_slopes %>%
  filter(guild %in% c('shade','intermediate','gap')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Light received slopes as a function of diameter', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'coefficients_lightbydiameter_evenshadegroups.png'), height = 6, width = 8, dpi = 300)

light_diam_slopes %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Light received slopes as a function of diameter', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'coefficients_lightbydiameter_quantileshadegroups.png'), height = 6, width = 8, dpi = 300)

# 6. Binned production by diameter

prod_bin_slopes <- lapply(allyears_names, function(x) {
  lm_list <- get(paste0(x, 'prod_lm'))
  lapply(lm_list, function(z) with(z, c(slope = as.numeric(z$coef[2]), cimin = confint(z)[2,1], cimax = confint(z)[2,2])))
})
for (i in 1:length(prod_bin_slopes)) prod_bin_slopes[[i]] <- data.frame(guild = guilds[i], do.call('rbind', prod_bin_slopes[[i]]))
prod_bin_slopes <- cbind(year = years, do.call('rbind',prod_bin_slopes))

prod_bin_slopes %>%
  filter(guild %in% c('shade','intermediate','gap')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Total (binned) production slopes', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'coefficients_totalproductionscaling_evenshadegroups.png'), height = 6, width = 9, dpi = 300)

prod_bin_slopes %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Total (binned) production slopes', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'coefficients_totalproductionscaling_quantileshadegroups.png'), height = 6, width = 9, dpi = 300)

# 7. Binned light received by diameter

light_bin_slopes_1990 <- lapply(names1990, function(x) {
  lm_dat <- get(paste0(x, 'light_lm'))
  c(slope = as.numeric(lm_dat$coef[2]), cimin = confint(lm_dat)[2,1], cimax = confint(lm_dat)[2,2])
})
light_bin_slopes_1990 <- data.frame(year = 1990, guild = guilds, do.call('rbind', light_bin_slopes_1990))

light_bin_slopes_1995 <- lapply(names1995, function(x) {
  lm_dat <- get(paste0(x, 'light_lm'))
  c(slope = as.numeric(lm_dat$coef[2]), cimin = confint(lm_dat)[2,1], cimax = confint(lm_dat)[2,2])
})
light_bin_slopes_1995 <- data.frame(year = 1995, guild = guilds, do.call('rbind', light_bin_slopes_1995))

light_bin_slopes <- rbind(light_bin_slopes_1990, light_bin_slopes_1995)

light_bin_slopes %>%
  filter(guild %in% c('shade','intermediate','gap')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Total (binned) light received slopes', 'Shade-tolerance groups based on even division of axis')

ggsave(file.path(fpfig, 'coefficients_totallightscaling_evenshadegroups.png'), height = 6, width = 8, dpi = 300)

light_bin_slopes %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Total (binned) light received slopes', 'Shade-tolerance groups based on quantiles')

ggsave(file.path(fpfig, 'coefficients_totallightscaling_quantileshadegroups.png'), height = 6, width = 8, dpi = 300)
