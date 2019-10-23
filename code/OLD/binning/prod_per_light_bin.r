# Plot production per light received as a function of tree diameter (BCI)
# Total and by shade/gap

# All points
alltreedat[[3]] %>%
  ggplot(aes(x = dbh_corr, y = production/light_received)) +
  geom_hex() +
  scale_fill_gradient(low='gray95', high='black') +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Production per unit light') +
  panel_border(colour = 'black')


fakebin_across_years <- function(dat_values, dat_classes, edges) {
  qprobs <- c(0.025, 0.5, 0.975)
  # add some padding just in case
  mins <- edges$bin_min
  mins[1] <- 0
  maxes <- edges$bin_max
  maxes[length(maxes)] <- Inf
  
  binstats <- t(sapply(1:length(mins), function(x) {
    indivs <- dat_values[dat_classes >= mins[x] & dat_classes < maxes[x]]
    c(mean = mean(indivs), quantile(indivs, probs = qprobs))
  }))
  dimnames(binstats)[[2]] <- c('mean', 'q025', 'median', 'q975')
  data.frame(bin_midpoint = edges$bin_midpoint,
             bin_min = edges$bin_min,
             bin_max = edges$bin_max,
             binstats)
}


# Bin
alltreelight <- alltreedat[[3]] %>% 
  filter(!is.na(light_received)) %>%
  mutate(energy_received = light_received * 86400 * 365.25)
prodperlight_alltree <- with(alltreelight, 
  fakebin_across_years(dat_values = 1e12 * production/energy_received, dat_classes = dbh_corr, edges = dbhbin_alltree))
  

prodperlight_alltree %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Production per unit light (mg MJ-1)') +
  panel_border(colour = 'black')
