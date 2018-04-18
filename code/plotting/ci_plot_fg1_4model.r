# Plot CI for functional group 1, all combinations, 1995
source('stan/extract_ci_stan.r')

densitybin_5census <- read.csv('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/densitybin_5census.csv', stringsAsFactors = FALSE)
dbh_pred_bins <- densitybin_5census[densitybin_5census$fg == 'all', 'bin_midpoint']

# Minimum x values for Pareto.
min_n <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/stanoutput/min_n.csv', stringsAsFactors = FALSE)

# Get confidence intervals
xmins <- min_n$xmin[min_n$year == 1995]
ns <- min_n$n[min_n$year == 1995]
ci_ppow <- dens_prod_ci(fit_ppow, dbh_pred = dbh_pred_bins, dens_form = 'pareto', prod_form = 'powerlaw', x_min = xmins[1], n_indiv = ns[2])
ci_pexp <- dens_prod_ci(fit_pexp, dbh_pred = dbh_pred_bins, dens_form = 'pareto', prod_form = 'powerlawexp', x_min = xmins[1], n_indiv = ns[2])
ci_wpow <- dens_prod_ci(fit_wpow, dbh_pred = dbh_pred_bins, dens_form = 'weibull', prod_form = 'powerlaw', x_min = NULL, n_indiv = ns[2])
ci_wexp <- dens_prod_ci(fit_wexp, dbh_pred = dbh_pred_bins, dens_form = 'weibull', prod_form = 'powerlawexp', x_min = NULL, n_indiv = ns[2])

ci_df <- rbind(data.frame(dens_model = 'pareto', prod_model = 'powerlaw', ci_ppow),
               data.frame(dens_model = 'pareto', prod_model = 'powerlawexp', ci_pexp),
               data.frame(dens_model = 'weibull', prod_model = 'powerlaw', ci_wpow),
               data.frame(dens_model = 'weibull', prod_model = 'powerlawexp', ci_wexp))

# Plot fits only.
library(cowplot)

area_core <- 42.84

ggplot(dplyr::filter(ci_df, variable == 'density'),
       aes(x = dbh, color = interaction(dens_model, prod_model), group = interaction(dens_model, prod_model))) +
  scale_x_log10() + scale_y_log10() +
  geom_line(aes(y = q50/area_core)) +
  geom_line(aes(y = q025/area_core), linetype = 'dotted') +
  geom_line(aes(y = q975/area_core), linetype = 'dotted')

ggplot(dplyr::filter(ci_df, variable == 'production'),
       aes(x = dbh, color = interaction(dens_model, prod_model), group = interaction(dens_model, prod_model))) +
  scale_x_log10() + scale_y_log10() +
  geom_line(aes(y = q50)) +
  geom_line(aes(y = q025), linetype = 'dotted') +
  geom_line(aes(y = q975), linetype = 'dotted')

ggplot(dplyr::filter(ci_df, variable == 'total_production'),
       aes(x = dbh, color = interaction(dens_model, prod_model), group = interaction(dens_model, prod_model))) +
  scale_x_log10() + scale_y_log10() +
  geom_line(aes(y = q50/area_core)) +
  geom_line(aes(y = q025/area_core), linetype = 'dotted') +
  geom_line(aes(y = q975/area_core), linetype = 'dotted')

# Plot data with the fits.
library(dplyr)

load('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/rawdataobj_22jan.r')
group_names <- c('all','all_classified','fg1','fg2','fg3','fg4','fg5','unclassified')
source('code/allfunctions27july.r')


fakebin_across_years <- function(dat_values, dat_classes, edges, mean_type = 'geometric', n_census = 5) {
  qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  # add some padding just in case
  mins <- edges$bin_min
  mins[1] <- 0
  maxes <- edges$bin_max
  maxes[length(maxes)] <- Inf
  
  binstats <- t(sapply(1:length(mins), function(x) {
    indivs <- dat_values[dat_classes >= mins[x] & dat_classes < maxes[x]]
    if (mean_type == 'geometric') {
      mean_n <- exp(mean(log(indivs)))
      sd_n <- exp(sd(log(indivs)))
      ci_width <- qnorm(0.975) * sd(log(indivs)) / sqrt(length(indivs))
      ci_min <- exp(mean(log(indivs)) - ci_width)
      ci_max <- exp(mean(log(indivs)) + ci_width)
      
    } else {
      mean_n <- mean(indivs)
      sd_n <- sd(indivs)
      ci_width <- qnorm(0.975) * sd(indivs) / sqrt(length(indivs))
      ci_min <- mean_n - ci_width
      ci_max <- mean_n + ci_width
    }
    c(mean = mean_n, 
      sd = sd_n,
      quantile(indivs, probs = qprobs),
      ci_min = ci_min,
      ci_max = ci_max)
  }))
  dimnames(binstats)[[2]] <- c('mean', 'sd', 'q025', 'q25', 'median', 'q75', 'q975', 'ci_min', 'ci_max')
  data.frame(bin_midpoint = edges$bin_midpoint,
             bin_min = edges$bin_min,
             bin_max = edges$bin_max,
             mean_n_individuals = edges$bin_count / n_census,
             binstats)
}


# Set number of bins       
numbins <- 20

# Bin all trees including unclassified
allyeardbh <- unlist(lapply(alltreedat[2:6], '[', , 'dbh_corr'))
dbhbin_all <- logbin(x = allyeardbh, y = NULL, n = numbins)
dbhbin_all_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_all))
# Make a version of alltreedat without the unclassified trees
alltreedat_classified <- lapply(alltreedat, function(x) subset(x, !is.na(fg)))
# Bin classified trees. (log binning of density)
allyeardbh_classified <- unlist(lapply(alltreedat_classified[2:6], '[', , 'dbh_corr'))
dbhbin_allclassified <- logbin(x = allyeardbh_classified, y = NULL, n = numbins)

dbhbin_fg_byyear <- list()

for (i in 1:6) {
  dbhbin_fg_byyear[[i]] <- lapply(fgdat[[i]][2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_allclassified))
}

fg1_dens_1995 <- dbhbin_fg_byyear[[1]][[2]]

fg1_prod_1995 <- fakebin_across_years(dat_values = fgdat[[1]][[3]]$production, dat_classes = fgdat[[1]][[3]]$dbh_corr, edges = fg1_dens_1995, n_census = 1)

totalprodbin_fg_byyear <- list()

for (i in 1:6) {
  totalprodbin_fg_byyear[[i]] <- lapply(fgdat[[i]][2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_allclassified))
}

fg1_totalprod_1995 <- totalprodbin_fg_byyear[[1]][[2]]

###

fpfig <- 'C:/Users/Q/google_drive/ForestLight/figs/credible_interval_plots'
global_diam_xlimits <- c(1, 316) # Maximum and minimum for x axis if it's diameter.

p_fg1_indivprod <- fg1_prod_1995 %>%
  filter(!is.na(mean), mean > 0) %>%
  ggplot() +
  geom_point(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975), color = 'black') +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
  panel_border(colour = 'black')

dat_prod <- dplyr::filter(ci_df, variable == 'production', dens_model == 'pareto')

p_fg1_indivprod +
  geom_ribbon(data = dat_prod, aes(ymin = q025, ymax = q975, x = dbh, fill = prod_model, group = prod_model), alpha = 0.5) +
  geom_line(data = dat_prod, aes(y = q50, x = dbh, color = prod_model, group = prod_model)) +
  theme(legend.position = 'bottom') + ggtitle('Individual production')
ggsave(file.path(fpfig, 'fg1_1995_production.png'), height=5, width=5, dpi=300)

p_fg1_dens <- fg1_dens_1995 %>%
  filter(!is.na(bin_value), bin_value > 0) %>%
  mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
  mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
  ggplot() +
  geom_point(aes(x = bin_midpoint, y = bin_yvalue), color = 'black') +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  panel_border(colour = 'black')

dat_dens <- dplyr::filter(ci_df, variable == 'density', prod_model == 'powerlaw')

p_fg1_dens +
  geom_ribbon(data = dat_dens, aes(ymin = q025/area_core, ymax = q975/area_core, x = dbh, fill = dens_model, group = dens_model), alpha = 0.5) +
  geom_line(data = dat_dens, aes(y = q50/area_core, x = dbh, color = dens_model, group = dens_model)) +
  theme(legend.position = 'bottom') + ggtitle('Density')
ggsave(file.path(fpfig, 'fg1_1995_density.png'), height = 5, width = 5, dpi = 300)

p_fg1_totalprod <- fg1_totalprod_1995 %>%
  filter(!is.na(bin_value), bin_value > 0) %>%
  mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
  mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
  ggplot() +
  geom_point(aes(x = bin_midpoint, y = bin_yvalue)) +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total production (kg ha'^-1,' y'^-1,')')), limits = c(.01,500)) +
  panel_border(colour = 'black')

dat_totalprod <- dplyr::filter(ci_df, variable == 'total_production')

p_fg1_totalprod +
  geom_ribbon(data = dat_totalprod, aes(ymin = q025/area_core, ymax = q975/area_core, x = dbh, fill = interaction(dens_model, prod_model), group = interaction(dens_model, prod_model)), alpha = 0.5) +
  geom_line(data = dat_totalprod, aes(y = q50/area_core, x = dbh, color = interaction(dens_model, prod_model), group = interaction(dens_model, prod_model))) +
  scale_fill_discrete(name = 'Functional forms', labels = c('D Pareto\nP Powerlaw', 'D Weibull\nP Powerlaw', 'D Pareto\nP Powerlaw*Exponential', 'D Weibull\nP Powerlaw*Exponential')) +
  theme(legend.position = 'bottom', legend.text = element_text(size=8)) + ggtitle('Total production') + guides(colour = FALSE)
ggsave(file.path(fpfig, 'fg1_1995_totalproduction.png'), height=5, width=6.5, dpi=300)
