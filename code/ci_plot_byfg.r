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

totalprodbin_fg_byyear <- list()

for (i in 1:6) {
  totalprodbin_fg_byyear[[i]] <- lapply(fgdat[[i]][2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_allclassified))
}

dens_plot_fg <- function(dat_dens_points, dat_dens_fits, fgname, year) {
  fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")
  years <- c(1990, 1995, 2000, 2005, 2010)
  p_all_dens <- dat_dens_points[[which(fgname == fg_names)]][[which(year == years)]] %>%
    filter(!is.na(bin_value), bin_value > 0) %>%
    mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
    mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
    ggplot() +
    geom_point(aes(x = bin_midpoint, y = bin_yvalue), color = 'black') +
    scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
    scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
    panel_border(colour = 'black')
  
  dat_dens <- dplyr::filter(dat_dens_fits, fg == fgname, variable == 'density', prod_model == 'powerlaw')
  
  p_all_dens +
    geom_ribbon(data = dat_dens, aes(ymin = q025/area_core, ymax = q975/area_core, x = dbh, fill = dens_model, group = dens_model), alpha = 0.5) +
    geom_line(data = dat_dens, aes(y = q50/area_core, x = dbh, color = dens_model, group = dens_model)) +
    theme(legend.position = 'bottom') + ggtitle('Density', 'contingent on production fit as power law')
  
}

# Density functional form is better for some groups than others.
dens_plot_fg(dbhbin_fg_byyear, ci_df, "fg1", 1995)
dens_plot_fg(dbhbin_fg_byyear, ci_df, "fg2", 1995)
dens_plot_fg(dbhbin_fg_byyear, ci_df, "fg3", 1995)
dens_plot_fg(dbhbin_fg_byyear, ci_df, "fg4", 1995)
dens_plot_fg(dbhbin_fg_byyear, ci_df, "fg5", 1995)

prod_plot_fg <- function(dat_prod_fits, dat_dens_points, fgname, year) {
  fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")
  years <- c(1990, 1995, 2000, 2005, 2010)
  fg_prod <- fakebin_across_years(dat_values = fgdat[[which(fgname == fg_names)]][[which(year == years) + 1]]$production, dat_classes = fgdat[[which(fgname == fg_names)]][[which(year == years) + 1]]$dbh_corr, edges = dat_dens_points[[which(fgname == fg_names)]][[which(year == years)]], n_census = 1)
  p_all_indivprod <- fg_prod %>%
    filter(!is.na(mean), mean > 0) %>%
    ggplot() +
    geom_point(aes(x = bin_midpoint, y = mean), color = 'black') +
    scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
    scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
    panel_border(colour = 'black')
  
  dat_prod <- dplyr::filter(dat_prod_fits, fg == 'alltree', variable == 'production', dens_model == 'pareto')
  
  p_all_indivprod +
    geom_ribbon(data = dat_prod, aes(ymin = q025, ymax = q975, x = dbh, fill = prod_model, group = prod_model), alpha = 0.5) +
    geom_line(data = dat_prod, aes(y = q50, x = dbh, color = prod_model, group = prod_model)) +
    theme(legend.position = 'bottom') + ggtitle('Individual production' , 'contingent on density fit as Pareto')
  
}

prod_plot_fg(ci_df, dbhbin_fg_byyear, "fg1", 1995)
prod_plot_fg(ci_df, dbhbin_fg_byyear, "fg2", 1995)
prod_plot_fg(ci_df, dbhbin_fg_byyear, "fg3", 1995)
prod_plot_fg(ci_df, dbhbin_fg_byyear, "fg4", 1995)
prod_plot_fg(ci_df, dbhbin_fg_byyear, "fg5", 1995)


totalprod_plot_fg <- function(dat_totalprod_points, dat_totalprod_fits, fgname, year) {
  fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")
  years <- c(1990, 1995, 2000, 2005, 2010)
  p_all_totalprod <- dat_totalprod_points[[which(fgname == fg_names)]][[which(year == years)]] %>%
    filter(!is.na(bin_value), bin_value > 0) %>%
    mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
    mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
    ggplot() +
    geom_point(aes(x = bin_midpoint, y = bin_yvalue)) +
    scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
    scale_y_log10(name = expression(paste('Total production (kg ha'^-1,' y'^-1,')')), limits = c(.1,10000)) +
    panel_border(colour = 'black')
  
  dat_totalprod <- dplyr::filter(dat_totalprod_fits, fg == fgname, variable == 'total_production')
  
  p_all_totalprod +
    geom_ribbon(data = dat_totalprod, aes(ymin = q025/area_core, ymax = q975/area_core, x = dbh, fill = interaction(dens_model, prod_model), group = interaction(dens_model, prod_model)), alpha = 0.5) +
    geom_line(data = dat_totalprod, aes(y = q50/area_core, x = dbh, color = interaction(dens_model, prod_model), group = interaction(dens_model, prod_model))) +
    scale_fill_discrete(name = 'Functional forms', labels = c('D Pareto\nP Powerlaw', 'D Weibull\nP Powerlaw', 'D Pareto\nP Powerlaw*Exponential', 'D Weibull\nP Powerlaw*Exponential')) +
    theme(legend.position = 'bottom', legend.text = element_text(size=8)) + ggtitle('Total production') + guides(colour = FALSE)
  
}

totalprod_plot_fg(totalprodbin_fg_byyear, ci_df, "fg1", 1995)
totalprod_plot_fg(totalprodbin_fg_byyear, ci_df, "fg2", 1995)
totalprod_plot_fg(totalprodbin_fg_byyear, ci_df, "fg3", 1995)
totalprod_plot_fg(totalprodbin_fg_byyear, ci_df, "fg4", 1995)
totalprod_plot_fg(totalprodbin_fg_byyear, ci_df, "fg5", 1995)
