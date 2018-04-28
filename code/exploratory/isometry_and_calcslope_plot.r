library(reshape2)
melt_pars <- melt(param_ci, id.vars=1:3)
cast_pars <- dcast(melt_pars, fg+year~parameter+variable)

# Version from 27 Apr. Manually find x and y locations for the slope segment to be plotted.

segment_location <- pred_light_5groups %>%
  group_by(fg, year) %>%
  summarize(xmax = sum(light_area[which.max(diff(q50)):(1+which.max(diff(q50)))])/2,
            ymax = sum(q50[which.max(diff(q50)):(1+which.max(diff(q50)))])/2)

cast_pars <- left_join(cast_pars, segment_location)

p_mean_panels <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_display)) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50), color = 'red') +
  geom_errorbar(aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max)) +
  geom_point(aes(x = bin_midpoint, y = mean)) +
  
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = xmax * 0.5, xend = xmax * 2, y = ymax * 0.5, yend = ymax * 2), color = 'green', size = 1) +
  scale_x_log10(name = title_x) + 
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = xmax * 0.5, xend = xmax * 2, y = ymax * 0.5^log_slope_q50 , yend = ymax * 2^log_slope_q50) , color = 'blue', size = 1) +
  scale_y_log10(name = title_y) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA))

ggsave('C:/Users/Q/google_drive/ForestLight/figs/5slopes.pdf', height = 8, width = 14)
