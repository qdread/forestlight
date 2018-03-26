# Credible intervals for all fits for all years.
# Edit 26 March. Get credible intervals for the parameters as well.

fp <- '~/forestlight/stanoutput'
fnames <- dir(fp, pattern = 'fit_')

source('~/forestlight/stancode/extract_ci_stan.r')

library(purrr)
library(tidyr)
library(dplyr)
library(rstan)

z <- data.frame(filename = fnames, stringsAsFactors = FALSE) %>%
  separate(filename, into = c('fit','model_name','fg','year','r'), remove = FALSE) %>%
  mutate(dens_model = if_else(model_name %in% c('ppow', 'pexp'), 'pareto', 'weibull'),
         prod_model = if_else(model_name %in% c('ppow','wpow'), 'powerlaw', 'powerlawexp')) %>%
  select(filename, year, dens_model, prod_model, fg)


min_n <- read.csv('~/forestlight/stancode/min_n.csv', stringsAsFactors = FALSE)
dbh_pred <- exp(seq(log(1.2), log(315), length.out = 50))

all_cis <- pmap(z, function(filename, year, fg, dens_model, prod_model) {
  x_min_i <- min_n$xmin[min_n$year == year & min_n$fg == fg]
  n_i <- min_n$n[min_n$year == year & min_n$fg == fg]
  load(file.path(fp, filename))
  ci <- dens_prod_ci(fit, dbh_pred, dens_form = dens_model, prod_form = prod_model, x_min = x_min_i, n_indiv = n_i)
  data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg, ci)
})

ci_df <- do.call(rbind, all_cis)
write.csv(ci_df, '~/forestlight/ci_by_fg.csv', row.names = FALSE)

# Credible intervals for parameters
pareto_par <- c('alpha')
weibull_par <- c('m', 'n')
powerlaw_par <- c('beta0', 'beta1')
powerlawexp_par <- c('beta0', 'beta1', 'a', 'b', 'c')

param_cis <- pmap(z, function(filename, year, fg, dens_model, prod_model) {
  load(file.path(fp, filename))
  if (dens_model == 'pareto') get_pars <- pareto_par else get_pars <- weibull_par
  if (prod_model == 'powerlaw') get_pars <- c(get_pars, powerlaw_par) else get_pars <- c(get_pars, powerlawexp_par)
  
  summ_fit <- summary(fit)
  cbind(data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg), summ_fit[[1]][get_pars, ])
})

param_cis <- map(param_cis, function(x) cbind(parameter = dimnames(x)[[1]], x))
param_cis <- do.call(rbind, param_cis)
names(param_cis)[9:13] <- c('q025', 'q25', 'q50', 'q75', 'q975')
write.csv(param_cis, '~/forestlight/paramci_by_fg.csv', row.names = FALSE)

# Draw plots (locally) ----------------------------------------------------

ci_df <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/ci_by_fg.csv', stringsAsFactors = FALSE)

# Load actual data points
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

all_dens_1995 <- dbhbin_all_byyear[[2]]

all_prod_1995 <- fakebin_across_years(dat_values = alltreedat[[3]]$production, dat_classes = alltreedat[[3]]$dbh_corr, edges = all_dens_1995, n_census = 1)

totalprodbin_alltree_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_all))

all_totalprod_1995 <- totalprodbin_alltree_byyear[[2]]

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

dens_plot_fg <- function(dat_dens_points, dat_dens_fits, fgname, year_to_plot) {
  fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")
  years <- c(1990, 1995, 2000, 2005, 2010)
  p_all_dens <- dat_dens_points[[which(fgname == fg_names)]][[which(year_to_plot == years)]] %>%
    filter(!is.na(bin_value), bin_value > 0) %>%
    mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
    mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
    ggplot() +
    geom_point(aes(x = bin_midpoint, y = bin_yvalue), color = 'black') +
    scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
    scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
    panel_border(colour = 'black')
  
  dat_dens <- dplyr::filter(dat_dens_fits, fg == fgname, variable == 'density', prod_model == 'powerlaw', year == year_to_plot)
  
  p_all_dens +
    geom_ribbon(data = dat_dens, aes(ymin = q025/area_core, ymax = q975/area_core, x = dbh, fill = dens_model, group = dens_model), alpha = 0.5) +
    geom_line(data = dat_dens, aes(y = q50/area_core, x = dbh, color = dens_model, group = dens_model)) +
    theme(legend.position = 'bottom') + ggtitle('Density')
  
}

library(cowplot)

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.
global_diam_xlimits <- c(1, 316) # Maximum and minimum for x axis if it's diameter.


# Density functional form is better for some groups than others.
dens_plot_fg(dbhbin_fg_byyear, ci_df, "fg1", 1995)
dens_plot_fg(dbhbin_fg_byyear, ci_df, "fg2", 1995)
dens_plot_fg(dbhbin_fg_byyear, ci_df, "fg3", 1995)
dens_plot_fg(dbhbin_fg_byyear, ci_df, "fg4", 1995)
dens_plot_fg(dbhbin_fg_byyear, ci_df, "fg5", 1995)


prod_plot_fg <- function(dat_prod_fits, dat_dens_points, fgname, year_to_plot) {
  fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")
  years <- c(1990, 1995, 2000, 2005, 2010)
  fg_prod <- fakebin_across_years(dat_values = fgdat[[which(fgname == fg_names)]][[which(year_to_plot == years) + 1]]$production, dat_classes = fgdat[[which(fgname == fg_names)]][[which(year_to_plot == years) + 1]]$dbh_corr, edges = dat_dens_points[[which(fgname == fg_names)]][[which(year_to_plot == years)]], n_census = 1)
  p_all_indivprod <- fg_prod %>%
    filter(!is.na(mean), mean > 0) %>%
    ggplot() +
    geom_point(aes(x = bin_midpoint, y = mean), color = 'black') +
    scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
    scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
    panel_border(colour = 'black')
  
  dat_prod <- dplyr::filter(dat_prod_fits, fg == 'alltree', variable == 'production', dens_model == 'pareto', year == year_to_plot)
  
  p_all_indivprod +
    geom_ribbon(data = dat_prod, aes(ymin = q025, ymax = q975, x = dbh, fill = prod_model, group = prod_model), alpha = 0.5) +
    geom_line(data = dat_prod, aes(y = q50, x = dbh, color = prod_model, group = prod_model)) +
    theme(legend.position = 'bottom') + ggtitle('Individual production')
  
}

prod_plot_fg(ci_df, dbhbin_fg_byyear, "fg1", 1995)
prod_plot_fg(ci_df, dbhbin_fg_byyear, "fg2", 1995)
prod_plot_fg(ci_df, dbhbin_fg_byyear, "fg3", 1995)
prod_plot_fg(ci_df, dbhbin_fg_byyear, "fg4", 1995)
prod_plot_fg(ci_df, dbhbin_fg_byyear, "fg5", 1995)

totalprod_plot_fg <- function(dat_totalprod_points, dat_totalprod_fits, fgname, year_to_plot) {
  fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")
  years <- c(1990, 1995, 2000, 2005, 2010)
  p_all_totalprod <- dat_totalprod_points[[which(fgname == fg_names)]][[which(year_to_plot == years)]] %>%
    filter(!is.na(bin_value), bin_value > 0) %>%
    mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
    mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
    ggplot() +
    geom_point(aes(x = bin_midpoint, y = bin_yvalue)) +
    scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
    scale_y_log10(name = expression(paste('Total production (kg ha'^-1,' y'^-1,')')), limits = c(.1,10000)) +
    panel_border(colour = 'black')
  
  dat_totalprod <- dplyr::filter(dat_totalprod_fits, fg == fgname, variable == 'total_production', year == year_to_plot)
  
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



# Write plots to files ----------------------------------------------------

fpfig <- 'C:/Users/Q/google_drive/ForestLight/figs/credible_interval_plots'

for (i in c('fg1','fg2','fg3','fg4','fg5')) ggsave(file.path(fpfig, paste0(i, '_1995_density.png')),
                                                   dens_plot_fg(dbhbin_fg_byyear, ci_df, i, 1995),
                                                   height = 5, width = 5, dpi = 300)

for (i in c('fg1','fg2','fg3','fg4','fg5')) ggsave(file.path(fpfig, paste0(i, '_1995_production.png')),
                                                   prod_plot_fg(ci_df, dbhbin_fg_byyear, i, 1995),
                                                   height = 5, width = 5, dpi = 300)

for (i in c('fg1','fg2','fg3','fg4','fg5')) ggsave(file.path(fpfig, paste0(i, '_1995_totalproduction.png')),
                                                   totalprod_plot_fg(totalprodbin_fg_byyear, ci_df, i, 1995),
                                                   height = 5, width = 6.5, dpi = 300)
