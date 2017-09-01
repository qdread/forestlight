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

allyear_dens_coeff %>%
  filter(guild %in% c('shade','intermediate','gap'), coefficient == 'L') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling cutoffs', 'Shade-tolerance groups based on even division of axis')

allyear_dens_coeff %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile'), coefficient == 'alpha') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max, colour = model, group = model)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling slopes', 'Shade-tolerance groups based on quantiles')

allyear_dens_coeff %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile'), coefficient == 'L') %>%
  ggplot(aes(x = guild, y = value, ymin = ci_min, ymax = ci_max)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  ggtitle('Density scaling cutoffs', 'Shade-tolerance groups based on quantiles')

# 2. Density by light received coefficients

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
  theme(strip.background = element_blank())

prod_slopes %>%
  filter(guild %in% c('shade_quantile','intermediate_quantile','gap_quantile')) %>%
  ggplot(aes(x = guild, y = slope, ymin = cimin, ymax = cimax, group = withindex, colour = withindex)) +
  facet_wrap(~ year) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank())

# 4. Individual production by light received

# 5. Individual light received by diameter

# 6. Binned production by diameter

# 7. Binned light received by diameter
