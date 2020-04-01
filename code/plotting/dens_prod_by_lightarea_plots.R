# Plots of density and production scaled by light per area (intermediate plots used to construct Figure 6A)
# QDR / ForestLight / 26 Mar 2020


# Load data ---------------------------------------------------------------

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users/jgradym/Google Drive/ForestLight'))

library(tidyverse)
library(rstan)
library(forestscaling)

#### Extract model output to get the fitted values, slopes, etc.
load('~/Dropbox/Q/projects/forestlight/fits_bylight_forratio.RData')

# source the extra extraction functions that aren't in the package
source('~/Documents/GitHub/forestscalingworkflow/R_functions/model_output_extraction_functions.r')
source('~/Documents/GitHub/forestlight/stan/get_ratio_slopes_fromfit.R')


# Get the statistics on the ratio trends.
la_pred <- logseq(1,412,101)
x_min <- 7
parnames_prod <- c('beta0', 'beta1')

# Predicted values for each fit.
dens_pred_fg_lightarea <- map(densfit_alltrees, function(fit) {
  pars_fg <- extract(fit, c('alpha')) %>% bind_cols
  pmap(pars_fg, pdf_pareto, x = la_pred, xmin = x_min)
})
prod_pred_fg_lightarea <- map(prodfit_alltrees, function(fit) {
  pars_fg <- extract(fit, parnames_prod) %>% bind_cols
  pmap(pars_fg, powerlaw_log, x = la_pred)
})

# Get total production including correction factors.

corr_factor <- function(y, y_fit, n_pars) {
  y_fit <- do.call(cbind, y_fit)
  # Get residuals by subtracting log y from linear predictor
  resids <- -1 * sweep(log(y_fit), 1, log(y))
  
  # Sum of squared residuals
  ssq_resid <- apply(resids^2, 2, sum)
  # Standard error of estimates
  sse <- (ssq_resid / (length(y) - n_pars))^0.5
  # Correction factors
  exp((sse^2)/2)
}

# We need all fitted values for production, not just the 101 values, to calculate correction factor.
# Recreate stan data
# Load data
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r')) # doesn't include imputed values

get_stan_data <- function(dat, x_min) with(dat, list(N = nrow(dat), x = dat$light_received_byarea, y = dat$production, x_min = x_min))

stan_data_list <- alltree_light_95 %>%
  filter(!recruit) %>%
  filter(!fg %in% 5, !is.na(fg), light_received_byarea >= x_min) %>%
  mutate(fg = paste0('fg', fg)) %>%
  group_by(fg) %>%
  group_map(~ get_stan_data(., x_min))

prod_pred_all_lightarea <- map2(stan_data_list, prodfit_alltrees, function(dat, fit) {
  pars_fg <- extract(fit, parnames_prod) %>% bind_cols
  pmap(pars_fg, powerlaw_log, x = dat$x)
})

prod_cf_fg_lightarea <- map2(stan_data_list, prod_pred_all_lightarea, ~ corr_factor(y = .x$y, y_fit = .y, n_pars = 2))

totalprod_pred_fg_lightarea <- map2(dens_pred_fg_lightarea, prod_pred_fg_lightarea, ~ do.call(cbind, .x) * do.call(cbind, .y)) %>%
  map2(prod_cf_fg_lightarea, ~ sweep(.x, 2, .y, `*`))


# Get credible intervals --------------------------------------------------

# Multiply density and total production times number of individuals, and divide by area
area_core <- 42.84
fg_names <- c('fast', 'tall', 'slow', 'short')

dens_pred_fg_lightarea_quantiles <- map2(dens_pred_fg_lightarea, map(stan_data_list, 'N'), function(dat, N) {
  do.call(cbind, dat) %>%
    sweep(2, N/area_core, `*`) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    t %>% as.data.frame %>% setNames(c('q025', 'q50', 'q975')) %>%
    mutate(light_area = la_pred)
})

dens_pred_dat <- map2_dfr(dens_pred_fg_lightarea_quantiles, fg_names,
                      ~ data.frame(fg = .y, .x)) %>%
  mutate(fg = factor(fg, levels = fg_names))

prod_pred_fg_lightarea_quantiles <- map(prod_pred_fg_lightarea, function(dat) {
  do.call(cbind, dat) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    t %>% as.data.frame %>% setNames(c('q025', 'q50', 'q975')) %>%
    mutate(light_area = la_pred)
}) 

prod_pred_dat <- map2_dfr(prod_pred_fg_lightarea_quantiles, fg_names,
                          ~ data.frame(fg = .y, .x)) %>%
  mutate(fg = factor(fg, levels = fg_names))

totalprod_pred_fg_lightarea_quantiles <- map2(totalprod_pred_fg_lightarea, map(stan_data_list, 'N'), function(dat, N) {
  dat %>%
    sweep(2, N/area_core, `*`) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    t %>% as.data.frame %>% setNames(c('q025', 'q50', 'q975')) %>%
    mutate(light_area = la_pred)
})

totalprod_pred_dat <- map2_dfr(totalprod_pred_fg_lightarea_quantiles, fg_names,
                          ~ data.frame(fg = .y, .x)) %>%
  mutate(fg = factor(fg, levels = fg_names))


# Plots -------------------------------------------------------------------

# Need to manually create observed data here: density, indiv. production, and total production *scaled by light per area*

data_to_bin <- alltree_light_95 %>%
  filter(fg %in% 1:4) %>%
  mutate(fg = factor(fg, labels = fg_names)) %>%
  select(fg, light_received_byarea, production)

# Determine bin edges by binning all
binedgedat <- with(data_to_bin, logbin(light_received_byarea, n = 20))

obs_dens <- data_to_bin %>%
  group_by(fg) %>%
  group_modify(~ logbin_setedges(x = .$light_received_byarea, edges = binedgedat))

obs_totalprod <- data_to_bin %>%
  group_by(fg) %>%
  group_modify(~ logbin_setedges(x = .$light_received_byarea, y = .$production, edges = binedgedat))

obs_indivprod <- data_to_bin %>%
  group_by(fg) %>%
  group_modify(~ cloudbin_across_years(dat_values = .$production, dat_classes = .$light_received_byarea, edges = binedgedat, n_census = 1))

# Themes
alpha_level <- 0.5
guild_colors <- c("black", "#99AF3C", "#1D5128", "#1e3160", "#6BA3BF", "#D6D6D6")
guild_fills <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "ivory")

theme_set(theme_bw() +
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.position = 'bottom'))

fill_scale <- scale_fill_manual(values = guild_fills[2:5], name = 'Functional guild', labels = fg_names, guide = guide_legend(override.aes = list(shape = 21)))
color_scale <- scale_color_manual(values = guild_fills[2:5], name = 'Functional guild', labels = fg_names, guide = FALSE)

# Density
p_dens <- ggplot() +
  geom_ribbon(data = dens_pred_dat %>% filter(light_area >= 7), aes(x = light_area, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = alpha_level) +
  geom_line(data = dens_pred_dat %>% filter(light_area >= 7), aes(x = light_area, y = q50, group = fg, color = fg)) +
  geom_point(data = obs_dens %>% filter(bin_count >= 20, bin_value > 0), aes(x = bin_midpoint, y = bin_value/area_core, group = fg, fill = fg), shape = 21, color = 'black', size = 2, show.legend = FALSE) +
  scale_x_log10(name = parse(text = 'Light~per~crown~area~(W~m^-2)'), limits=c(7,400)) +
  scale_y_log10(name = parse(text = 'Abundance~(ha^-1~cm^-1)')) +
  fill_scale +
  color_scale

# Production
p_indivprod <- ggplot() +
  geom_ribbon(data = prod_pred_dat %>% filter(light_area >= 7), aes(x = light_area, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = alpha_level) +
  geom_line(data = prod_pred_dat %>% filter(light_area >= 7), aes(x = light_area, y = q50, group = fg, color = fg)) +
  geom_point(data = obs_indivprod %>% filter(mean_n_individuals >= 20), aes(x = bin_midpoint, y = median, group = fg, fill = fg), , shape = 21, color = 'black', size = 2, show.legend = FALSE) +
  scale_x_log10(name = parse(text = 'Light~per~crown~area~(W~m^-2)'), limits=c(7,400)) +
  scale_y_log10(name = parse(text = 'Individual~growth~(kg~y^-1)')) +
  fill_scale +
  color_scale

# Total production
p_totalprod <- ggplot() +
  geom_ribbon(data = totalprod_pred_dat %>% filter(light_area >= 7), aes(x = light_area, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = alpha_level) +
  geom_line(data = totalprod_pred_dat %>% filter(light_area >= 7), aes(x = light_area, y = q50, group = fg, color = fg)) +
  geom_point(data = obs_totalprod %>% filter(bin_count >= 20, bin_value > 0), aes(x = bin_midpoint, y = bin_value/area_core, group = fg, fill = fg), shape = 21, color = 'black', size = 2, show.legend = FALSE) +
  scale_x_log10(name = parse(text = 'Light~per~crown~area~(W~m^-2)'), limits = c(7,400)) +
  scale_y_log10(name = parse(text = 'Total~production~(kg~y^-1~ha^-1~cm^-1)')) +
  fill_scale +
  color_scale

# Save plots
fp_fig <- file.path(gdrive_path, 'figs/light_area_plots_mar2020')

ggsave(file.path(fp_fig, 'density_scaled_lightperarea.png'), p_dens, height = 5, width = 5, dpi = 300)
ggsave(file.path(fp_fig, 'production_scaled_lightperarea.png'), p_indivprod, height = 5, width = 5, dpi = 300)
ggsave(file.path(fp_fig, 'totalproduction_scaled_lightperarea.png'), p_totalprod, height = 5, width = 5, dpi = 300)


# Table of exponents ------------------------------------------------------

dens_table <- densfit_summaries %>%
  map2_dfr(fg_names, ~ data.frame(fg = .y, parameter = row.names(.x), as.data.frame(.x))) %>%
  filter(!parameter %in% 'lp__') %>%
  setNames(c('fg','parameter','mean','se_mean','sd',paste0('q',c('025','25','50','75','975')), 'n_eff', 'Rhat')) %>%
  mutate_at(vars(mean, starts_with('q')), ~ - (. + 1)) %>%
  mutate(parameter = 'density slope')

prod_table <- prodfit_summaries %>%
  map2_dfr(fg_names, ~ data.frame(fg = .y, parameter = c('growth intercept', 'growth slope', 'growth st.dev.', 'lp__'), as.data.frame(.x))) %>%
  filter(!parameter %in% 'lp__') %>%
  setNames(c('fg','parameter','mean','se_mean','sd',paste0('q',c('025','25','50','75','975')), 'n_eff', 'Rhat')) %>%
  arrange(parameter)

# Manually get total production slope quantiles by taking beta1 - alpha - 1 for each iteration, then take mean and quantiles

totalprod_slopes <- map2(densfit_alltrees, prodfit_alltrees, function(dfit, pfit) {
  extract(pfit, 'beta1')[[1]] - extract(dfit, 'alpha')[[1]] - 1
})

totalprod_table <- totalprod_slopes %>%
  map2_dfr(fg_names, ~ data.frame(fg = .y, 
                                  parameter = 'total production slope',
                                  mean = mean(.x),
                                  se_mean = sd(.x)/sqrt(length(.x)),
                                  sd = sd(.x),
                                  t(quantile(.x, probs = c(.025, .25, .5, .75, .975))),
                                  n_eff = NA,
                                  Rhat = NA)) %>%
  setNames(c('fg','parameter','mean','se_mean','sd',paste0('q',c('025','25','50','75','975')), 'n_eff', 'Rhat')) 

final_table <- bind_rows(dens_table, prod_table, totalprod_table) %>%
  select(-q25, -q75, -n_eff, -Rhat, -sd)

write_csv(final_table, file.path(gdrive_path, 'data/clean_summary_tables/slopes_lightperarea.csv'))
