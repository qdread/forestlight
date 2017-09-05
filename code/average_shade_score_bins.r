# Weighted average of pca score plotted versus light received per unit crown area.

ggplot(alltree_light_90, aes(x=log10(light_received/crownarea), y=pca_scores)) + geom_hex() # Incorrect.

# Create log bin, then somehow weight by the pca score of each individual tree.

numbins <- 20

lightperareabin_1990 <- with(alltree_light_90, logbin(x=light_received/crownarea, y=NULL, n=numbins))

# For each tree in 1990, determine which of the 20 light per area bins it belongs to.
whichbin1990 <- sapply(with(alltree_light_90, light_received/crownarea), function(x) which(x >= lightperareabin_1990$bin_min & x < lightperareabin_1990$bin_max))

shadescores_1990 <- alltree_light_90 %>%
  mutate(light_area_bin = whichbin1990) %>%
  group_by(light_area_bin) %>%
  summarize(shade_score = mean(na.omit(pca_scores)),
            shade_ci = qnorm(0.975) * sd(na.omit(pca_scores))/sum(!is.na(pca_scores)))

lightperareabin_1990 <- cbind(lightperareabin_1990, shadescores_1990)
fit_1990 <- lm(shade_score ~ I(log10(bin_midpoint)), data = lightperareabin_1990)

label_df <- data.frame(bin_midpoint = 50, shade_score = 0.6, lab = paste0('slope = ', round(fit_1990$coefficients[2],2), ', r2 = ', round(summary(fit_1990)$r.sq, 2)))

p1990 <- ggplot(lightperareabin_1990, aes(x = bin_midpoint, y = shade_score)) +
  geom_pointrange(aes(ymin = shade_score - shade_ci, ymax = shade_score + shade_ci)) +
  stat_smooth(method = 'lm', se = FALSE) +
  geom_text(data=label_df, aes(label=lab)) +
  scale_x_log10(name = expression(paste('Light received (W m'^-2,')'))) +
  labs(y = 'Average shade tolerance of bin') +
  ggtitle('1990') +
  panel_border(colour='black')


lightperareabin_1995 <- with(alltree_light_95, logbin(x=light_received/crownarea, y=NULL, n=numbins))

# For each tree in 1995, determine which of the 20 light per area bins it belongs to.
whichbin1995 <- sapply(with(alltree_light_95, light_received/crownarea), function(x) which(x >= lightperareabin_1995$bin_min & x < lightperareabin_1995$bin_max))

shadescores_1995 <- alltree_light_95 %>%
  mutate(light_area_bin = whichbin1995) %>%
  group_by(light_area_bin) %>%
  summarize(shade_score = mean(na.omit(pca_scores)),
            shade_ci = qnorm(0.975) * sd(na.omit(pca_scores))/sum(!is.na(pca_scores)))

lightperareabin_1995 <- cbind(lightperareabin_1995, shadescores_1995)
fit_1995 <- lm(shade_score ~ I(log10(bin_midpoint)), data = lightperareabin_1995)

label_df <- data.frame(bin_midpoint = 50, shade_score = 0.6, lab = paste0('slope = ', round(fit_1995$coefficients[2],2), ', r2 = ', round(summary(fit_1995)$r.sq, 2)))

p1995 <- ggplot(lightperareabin_1995, aes(x = bin_midpoint, y = shade_score)) +
  geom_pointrange(aes(ymin = shade_score - shade_ci, ymax = shade_score + shade_ci)) +
  stat_smooth(method = 'lm', se = FALSE) +
  geom_text(data=label_df, aes(label=lab)) +
  scale_x_log10(name = expression(paste('Light received (W m'^-2,')'))) +
  labs(y = 'Average shade tolerance of bin') +
  ggtitle('1995') +
  panel_border(colour='black')

grid2 <- plot_grid(p1990,p1995)
ggsave(file.path(fpfig, 'shadescore_by_lightperarea.png'), height=5, width=8, dpi=300)
