area_core <- 42.84

predbin_all <- predbin_all %>%
  rename(dbh = bin_midpoint) %>% 
  mutate(fg = if_else(fg=='alltree','all',fg)) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core))

plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg3'),
               model_fit_density = 'weibull', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 500),
               y_breaks = c(0.1, 1, 10, 100),
               preddat = predbin_all)
plot_totalprod(year_to_plot = 1995,
                fg_names = c('fg3'),
                model_fit_density = 'weibull', 
                model_fit_production = 'powerlaw',
                y_limits = c(0.01, 500),
                y_breaks = c(0.1, 1, 10, 100)
                )
