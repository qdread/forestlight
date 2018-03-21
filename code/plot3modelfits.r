# Check diagnostics and plot.

library(rstan)
load('C:/Users/Q/Dropbox/projects/forestlight/stanoutput/localsubsamplefit5000_20Mar_alltree.RData')

summary(fit_ppow_all)[[1]]
summary(fit_pexp_all)[[1]]
summary(fit_pbert_all)[[1]]
summary(fit_wpow_all)[[1]]
summary(fit_wexp_all)[[1]]
summary(fit_wbert_all)[[1]]


library(bayesplot)

mcmc_trace(as.array(fit_ppow_all))
mcmc_trace(as.array(fit_pexp_all))
mcmc_trace(as.array(fit_pbert_all))
mcmc_trace(as.array(fit_wpow_all))
mcmc_trace(as.array(fit_wexp_all))
mcmc_trace(as.array(fit_wbert_all))

# Plot raw data from 1995 with the fits on top.
library(cowplot)


p_prod <- ggplot() +
  geom_hex(data = alltreedat[[3]], aes(x = dbh_corr, y = production)) +
  scale_fill_continuous(low = 'gray80', high = 'black') +
  scale_x_log10() + scale_y_log10()

library(dplyr)

dat_prod <- ci_df %>%
  filter(dens_model == 'weibull', variable == 'production', fg == 'alltree')
  

p_prod +
  geom_line(data = dat_prod, aes(x = dbh, y = q50, color = prod_model)) +
  geom_line(data = dat_prod, aes(x = dbh, y = q025, color = prod_model), linetype = 'dotted') +
  geom_line(data = dat_prod, aes(x = dbh, y = q975, color = prod_model), linetype = 'dotted')

ggsave('C:/Users/Q/google_drive/ForestLight/figs/credible_interval_plots/3modelindivproduction.png', height=6, width=7.5, dpi=400)

# Plot total production data from 1995 with the fits on top

p_totalprod <- ggplot() +
  geom_point(data = all_totalprod_1995 %>% mutate(bin_value = bin_value/42.84), aes(x = bin_midpoint, y = bin_value)) +
  scale_x_log10() + scale_y_log10()

dat_totprod <- ci_df %>%
  filter(dens_model == 'weibull', variable == 'total_production', fg == 'alltree') %>%
  mutate_at(vars(starts_with('q')), funs(./42.84))

p_totalprod +
  geom_line(data = dat_totprod, aes(x = dbh, y = q50, color = prod_model)) +
  geom_line(data = dat_totprod, aes(x = dbh, y = q025, color = prod_model), linetype = 'dotted') +
  geom_line(data = dat_totprod, aes(x = dbh, y = q975, color = prod_model), linetype = 'dotted')

ggsave('C:/Users/Q/google_drive/ForestLight/figs/credible_interval_plots/3modeltotalproduction.png', height=6, width=7.5, dpi=400)
