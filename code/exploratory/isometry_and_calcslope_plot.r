library(reshape2)
melt_pars <- melt(param_ci, id.vars=1:3)
cast_pars <- dcast(melt_pars, fg+year~parameter+variable)

p_mean_panels <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_display)) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975), alpha = 0.5, fill = 'red') +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50), color = 'red') +
  geom_errorbar(aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max)) +
  geom_point(aes(x = bin_midpoint, y = mean)) +
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = x_max_q50 * 0.5, xend = x_max_q50 * 2, y = (4/9) * k_q50 * G_q50 * (x_max_q50 * 0.5), yend = (4/9) * k_q50 * G_q50 * (x_max_q50 * 2)), color = 'blue', size = 1) +
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = x_max_q025 * 0.5, xend = x_max_q025 * 2, y = (4/9) * k_q025 * G_q025 * (x_max_q025 * 0.5), yend = (4/9) * k_q025 * G_q025 * (x_max_q025 * 2)), color = 'blue', size = 1, linetype = 'dashed') +
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = x_max_q975 * 0.5, xend = x_max_q975 * 2, y = (4/9) * k_q975 * G_q975 * (x_max_q975 * 0.5), yend = (4/9) * k_q975 * G_q975 * (x_max_q975 * 2)), color = 'blue', size = 1, linetype = 'dashed') +
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = x_max_q50 * 0.5, xend = x_max_q50 * 2, y = (y_max_q50 * 0.5), yend = (y_max_q50 * 2)), color = 'green', size = 1) +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA))


# Version from 27 Apr. Manually find x and y locations for the slope segment to be plotted.

segment_location <- pred_light_5groups %>%
  group_by(fg, year) %>%
  summarize(xmax = sum(light_area[which.max(diff(q50)):(1+which.max(diff(q50)))])/2,
            ymax = sum(q50[which.max(diff(q50)):(1+which.max(diff(q50)))])/2)

cast_pars <- left_join(cast_pars, segment_location)

p_mean_panels <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_display)) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975), alpha = 0.5, fill = 'red') +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50), color = 'red') +
  geom_errorbar(aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max)) +
  geom_point(aes(x = bin_midpoint, y = mean)) +
  geom_abline(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(slope=log_slope_q50, intercept=-2)) +
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = xmax * 0.5, xend = xmax * 2, y = ymax * 0.5 , yend = ymax * 2) , color = 'blue', size = 1) +
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = xmax * 0.5, xend = xmax * 2, y = ymax * 0.5, yend = ymax * 2), color = 'green', size = 1) +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA))
