# Clean plotting code using "packaged" functions
# Plots all main and supplemental figures in forest scaling MS

########### ============== ##############
### WHICH MODEL FITS TO USE IN PLOTS ####
########### ============== ##############

# change these if needed.
DENS = 3
PROD = 1


# Set path to data on google drive
#devtools::install_github('qdread/forestscaling')

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users/jgradym/Google Drive/ForestLight'))

library(forestscaling) # Packaged all the functions and ggplot2 themes here!
library(tidyverse)
library(egg)
library(scales)
library(RColorBrewer)
library(gtable)
library(grid)
library(reshape2)
library(hexbin)

# Define color schemes and labels
guild_fills <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "ivory")
guild_fills <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray87")
guild_fills2 <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray")
guild_colors <- c("black", "#99AF3C", "#1D5128", "#1e3160", "#6BA3BF", "#D6D6D6")
guild_fills_nb <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "ivory")
guild_fills_nb <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray87")
guild_fills_nb0 <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")
guild_colors_nb0 <- c("#99AF3C", "#1D5128", "#1e3160", "#6BA3BF", "#D6D6D6")
guild_colors_nb <- c("#3B4403", "#02330A", "#031a49", "#02394F", "#595A5B")
guild_fills_nb <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")

fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('Fast','LL Pioneer', 'Slow', 'SL Breeder', 'Medium')

geom_size <- 4

# Some Plotting Code


################################################################################################
# ------------------------------ Fig 1: hand drawn schematics ---------------------------------
################################################################################################



################################################################################################
# ------------------------------ Fig 2: Life Histories ---------------------------------
################################################################################################

# ### Fig 2a, PCA
# Load Nadja's data (new functional groups 25 June)
# fg5 is the new column (we originally used fg from the older df)
fgbci <- read.table(file.path(gdrive_path, 'data/Ruger/fgroups_dynamics_new.txt'), stringsAsFactors = FALSE)
# Correct functional groups so that: 1 fast, 2 pioneer, 3 slow, 4 breeder, 5 intermediate
# Old 1,2,3,4,5 --> New 2,3,1,4,5
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5))

# Currently X1new is high for slow species and low for fast species
# Currently X2new is high for pioneer species and low for breeder species
# Correct these
fgbci$PC_slow_to_fast <- -fgbci$X1new
fgbci$PC_breeder_to_pioneer <- fgbci$X2new


#Classification description: 

# Quentin: I went through my email and found the explanation Nadja gave for how groups are defined. 
# So it's basically if you took a circle out of the middle of a pizza, then cut the remaining ring into four 90 degree slices. 
#I think the percentiles divided by 2 is a fine cutoff, since ultimately it is arbitrary, 
#and the way it's split seems fine. In any case all the analyses we did using the continuous variables show the same result, qualitatively.

# Nadja: I attach the new classification, it is column 'fg5'.
#The values to split between groups are -0.951 and 1.284 on the first axis 
#(these are 10th and 90th percentiles divided by 2), on axis 2 it's -0.976 and 1.097.
#Then, I divided species scores by the absolute value of these values and 
#included all species within a circle of radius 1 in the 'intermediate' group. 
#Otherwise, species were split at 45° angles between axes.
# Load data ----
lab_x <- expression(paste(italic('Slow'), 'to', italic('Fast')))
fgs <- ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, fill = factor(fg5))) +
  geom_point(shape = 21, size = geom_size, color = "black") + 
  labs(x = 'Slow to Fast', y = 'Breeders to Pioneers') +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_color_manual(values = guild_colors_nb, labels = fg_labels, name = 'functional group')+
  scale_fill_manual(values = guild_fills_nb) + theme_plant()
fgs

p2a  <- set_panel_size(fgs, width=unit(14.3,"cm"), height=unit(14.3,"cm"))
grid.newpage()
grid.draw(p2a)
pdf(file.path(gdrive_path, 'Figures/Fig_2/Fig_2a.pdf'))
grid.draw(p2a)
dev.off()

#########################################3

#########  Fig 2b
# Axis titles
title_x <- expression(paste('Light per Crown Area (W m'^-2,')',sep=''))
title_y <- expression(atop('Growth per Crown Area', paste('(kg yr'^-1, ' m'^-2,')', sep='')))
scale_y_log10(name =  expression(atop('Growth per Crown Area',
                                      paste('(kg y'^-1, ' m'^-2,')'))))
# Colors
colors <- c('black',"#BFE046","#267038" , "#27408b", "#87Cefa", "gray87" ) # #BFE046= light green,#267038=dark green,27408b=dark blue,"#87Cefa" = light blue,    
year_to_plot = 1995
fg_colors <- c(fg1 = "#BFE046", fg2 = "#27408b" , fg3 = "#267038", fg4 =  "#87Cefa", fg5 = "gray70",  alltree = 'black') #unclassified = 'brown',
fg_colors2 <- c(fg1 = "#BFE046", fg2 = "#27408b" , fg3 = "#267038", fg4 =  "#87Cefa", fg5 = "gray87",  alltree = 'black') #unclassified = 'brown',

fp <- file.path(gdrive_path, 'data/data_forplotting')

# New File path needed
obs_light_binned <- read.csv(file.path(fp, 'obs_light_binned.csv'), stringsAsFactors = FALSE)
obs_light_raw <- read.csv(file.path(fp, 'obs_light_raw.csv'), stringsAsFactors = FALSE)
pred_light <- read.csv(file.path(fp, 'pred_light.csv'), stringsAsFactors = FALSE)
param_ci <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/lightbyarea_paramci_by_fg.csv'), stringsAsFactors = FALSE)

# Get rid of the predicted points that are outside the limits of the observed data for each FG
obs_limits <- obs_light_binned %>%
  group_by(fg, year) %>%
  summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))

pred_light <- pred_light %>%
  left_join(obs_limits) %>%
  filter(light_area >= min_obs & light_area <= max_obs)

pred_light_5groups <- pred_light %>% filter(!fg %in% c('alltree','unclassified'))

melt_pars <- melt(param_ci, id.vars=1:3)
cast_pars <- dcast(melt_pars, fg+year~parameter+variable)

dodge_width <- 0.03
error_bar_width <- 0.04

# Do some additional computation to correct the error bar width for the number of groups in each bin
obs_light_binned_plotdata <- obs_light_binned %>% filter(year == year_to_plot, mean_n_individuals >= 20, !fg %in% c('alltree', 'unclassified')) %>%
  group_by(bin_midpoint, year) %>% 
  mutate(width = sum(c('fg1','fg2','fg3','fg4','fg5') %in% fg)) %>% 
  ungroup

p_mean_1panel <- ggplot(obs_light_binned_plotdata) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
              aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.4) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot),
            aes(x = light_area, y = q50, color = fg)) +
  # Comment out the following line to remove error bars, or change ci_min and ci_max to q25 and q75 to use quantiles instead of the CI of mean.
  geom_errorbar(aes(x = bin_midpoint, ymin = q25, ymax = q75, 
                    group = fg, color = fg, width = 0), #width = error_bar_width * width), 
                position = position_dodge(width = dodge_width)) + 
  geom_point(aes(x = bin_midpoint, y = mean, group = fg, fill = fg),
             size = 4, shape = 21, position = position_dodge(width = dodge_width)) +
  
  scale_x_log10(name = title_x, limits = c(1.5, 412)) + 
  scale_y_log10(name = title_y, position = "right", breaks = c(0.01, 0.03, 0.1, 0.3), 
                labels = c( 0.01, 0.03, 0.1, 0.3)) +
  scale_color_manual(name = 'Functional group', values = guild_fills_nb0, labels = fg_labels) +
  scale_fill_manual(values = guild_fills_nb, labels = fg_labels, guide = FALSE) +
  theme_plant() + theme_no_x()
p2b <- set_panel_size(p_mean_1panel, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p2b)

pdf(file.path(gdrive_path, "Figures/Fig_2/Fig_2b.pdf"))
grid.draw(p2b)
dev.off()

### ################## Fig 2c, Mortality


fg_labels <- c('Fast','LL Pioneer', 'Slow', 'SL Breeder', 'Medium')
guild_fills_nb <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray90")
guild_colors_nb <- c("#3B4403", "#02330A", "#031a49", "#02394F", "gray")

bin_mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_mortalitybins.csv')) # Load binned data
mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_mortalityindividuals.csv')) # Load raw data
fitted_mort <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/mortality_ci_by_fg.csv')) # Load fitted values

library(broom)

mort_slopes <- fitted_mort %>%
  filter(light_per_area > 10) %>%
  filter(light_per_area < 100)%>%
  nest(-fg) %>% #group variable
  mutate(
    fit = map(data, ~ lm(log(q50)/5 ~ log(light_per_area), data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)%>%
  filter(term != '(Intercept)')
mort_slopes


growth_slopes <- obs_light_binned_plotdata %>%
  filter( bin_midpoint > 10,  bin_midpoint < 100) %>%
  nest(-fg) %>% #group variable
  mutate(
    fit = map(data, ~ lm(log(mean) ~ log(bin_midpoint), data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)%>%
  filter(term != '(Intercept)')
growth_slopes


# Truncate fitted mortality lines to not exceed the range of the observed data, using 20 individuals as the cutoff.
obs_range_mort <- bin_mort %>% 
  filter(variable %in% 'light_per_area', lived + died >= 20) %>%
  group_by(fg) %>%
  summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))

fitted_mort_trunc <- fitted_mort %>%
  left_join(obs_range_mort) %>%
  filter(light_per_area >= min_obs & light_per_area <= max_obs)

geom_size = 4
#Mortality
plightarea <- ggplot(data = fitted_mort_trunc %>% mutate(fg = factor(fg, labels = fg_labels))) +
  geom_ribbon(aes(x = light_per_area, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
  geom_line(aes(x = light_per_area, y = q50, group = fg, color = fg)) +
  geom_point(data = bin_mort %>% filter(variable == 'light_per_area', !fg %in% c('all','unclassified'), (lived+died) >= 20)  %>% 
               mutate(fg = factor(fg, labels = fg_labels)),
             aes(x = bin_midpoint, y = mortality, fill = fg),
             shape = 21, size = geom_size) +
  scale_x_log10(name = parse(text = 'Light~per~Crown~Area~(W~m^-2)'), breaks = c(3, 30, 300), limits = c(1.5, 412)) +
  scale_y_continuous(trans = "logit", position = "right", breaks = c(0.03, 0.1, 0.3, 0.6), 
                     labels = c(0.03, 0.1, 0.3, 0.6), limits = c(0.02, 0.65),
                name = expression(paste("Mortality (5 yr"^-1,")"))) +
  scale_color_manual(values = guild_fills_nb) +
  scale_fill_manual(values = guild_fills_nb) +
  theme_plant()


p2c <- set_panel_size(plightarea, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p2c)  #notice that LLP and Slow show steep slow of mortality with light!

pdf(file.path(gdrive_path, "Figures/Fig_2/Fig_2c.pdf"))
grid.draw(p2c)
dev.off()


################################################################################################
# ------------------------Fig 3: Plotting for Piecewise Scaling --------------------------------
################################################################################################

# Plot of slopes in different segments by different functional groups.

fp_plot <- file.path(gdrive_path, 'data/data_forplotting')

# Plot the fitted values on top of the observed histograms.

# Read observed data


for (i in dir(fp_plot, pattern = 'obs_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp_plot, i), stringsAsFactors = FALSE))
}

# Read modeled data (CIs)

for (i in dir(fp_plot, pattern = 'pred_|fitted_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp_plot, i), stringsAsFactors = FALSE))
}

# Create plots.
#Model fit 1 = pareto, 1 segment
#Model Fit 2  = 2 segments, etc

grob_text_a <- grobTree(textGrob("a", x = 0.06, y = 0.92, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
grob_text_b <- grobTree(textGrob("b", x = 0.06, y = 0.88, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
grob_text_c <- grobTree(textGrob("c", x = 0.05, y = 0.95, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
geom_size = 2
obs_dens <- obs_dens %>%
  filter(bin_count >= 20)
plot_dens2 <- function (year_to_plot = 1995, fg_names = c("fg1", "fg2", "fg3", 
                                                          "fg4", "fg5", "all"), model_fit = 1, x_limits, x_breaks = c(1, 
                                                                                                                      3, 10, 30, 100, 300), y_limits, y_breaks, y_labels, fill_names = c("black", 
                                                                                                                                                                                         "#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), color_names = c("black", 
                                                                                                                                                                                                                                                                "#BFE046", "#267038", "#27408b", "#87Cefa", "gray"), x_name = "Diameter (cm)", 
                        y_name = expression(paste("Density (n ha"^-1, "cm"^-1, ")")), 
                        geom_size = 3.5, obsdat = obs_dens, preddat = pred_dens, plot_abline = TRUE, 
                        abline_slope = -2, abline_intercept = 4) 
{
  obsdat <- obsdat %>% dplyr::filter(fg %in% fg_names, year == 
                                       year_to_plot, bin_count > 10) %>% dplyr::filter(bin_value > 
                                                                                         0)
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% dplyr::summarize(min_obs = min(bin_midpoint), 
                                                                    max_obs = max(bin_midpoint))
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% dplyr::filter(dens_model %in% 
                                                                          model_fit, fg %in% fg_names, year == year_to_plot) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% dplyr::filter(dbh >= 
                                                                             min_obs & dbh <= max_obs)
  p <- ggplot2::ggplot() + ggplot2::geom_ribbon(data = preddat, 
                                                ggplot2::aes(x = dbh, ymin = q025, ymax = q975, group = fg, 
                                                             fill = fg), alpha = 0.4)
  if (plot_abline) {
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", linetype = "dashed", 
                                  size = 0.75)
  }
  p + ggplot2::geom_line(data = preddat, ggplot2::aes(x = dbh, 
                                                      y = q50, group = fg, color = fg)) + ggplot2::geom_line(data = preddat[preddat$fg == 
                                                                                                                              "fg5", ], ggplot2::aes(x = dbh, y = q50), color = "gray") + 
    ggplot2::geom_point(data = obsdat, ggplot2::aes(x = bin_midpoint, 
                                                    y = bin_value, group = fg, fill = fg), size = geom_size, 
                        shape = 21, color = "black") + ggplot2::scale_x_log10(name = x_name, 
                                                                              limits = x_limits, breaks = x_breaks) + ggplot2::scale_y_log10(name = y_name, 
                                                                                                                                             limits = y_limits, breaks = y_breaks, labels = y_labels) + 
    ggplot2::scale_color_manual(values = color_names) + ggplot2::scale_fill_manual(values = fill_names) + 
    theme_plant()
}
p <- plot_dens(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = DENS,
          x_limits = c(.8, 230),
          y_limits = c(0.01, 20000),
          x_breaks = c(1, 10, 100),
          y_labels = c(0.001, 0.1, 10,1000),
          y_breaks = c(0.001, 0.1,  10, 1000))
#p <- p +annotation_custom(grob_text_b)

p1 <- set_panel_size(p, width=unit(8,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path,'Figures/Fig_3/Main/Fig_3b_Density.pdf'))
grid.draw(p1)
dev.off()

# Specify dodging with a certain width of error bar
# Model fit 1 = power law
# Model fit 2 = power law exp

plot_prod2 <- function (year_to_plot = 1995, fg_names = c("fg1", "fg2", "fg3", 
                                                          "fg4", "fg5", "all"), model_fit = 1, x_limits, x_breaks = c(1, 
                                                                                                                      3, 10, 30, 100, 300), y_limits, y_labels, y_breaks, fill_names = c("#BFE046", 
                                                                                                                                                                                         "#267038", "#27408b", "#87Cefa", "gray87"), color_names = c("#BFE046", 
                                                                                                                                                                                                                                                     "#267038", "#27408b", "#87Cefa", "gray"), x_name = "Diameter (cm)", 
                        y_name = expression(paste("Growth (kg y"^-1, ")")), average = "mean", 
                        plot_errorbar = FALSE, error_min = "ci_min", error_max = "ci_max", 
                        error_bar_width = 0.1, error_bar_thickness = 0.5, dodge_width = 0.03, 
                        dodge_errorbar = TRUE, geom_size = 4, obsdat = obs_indivprod, 
                        preddat = fitted_indivprod, plot_abline = TRUE, abline_slope = 2, 
                        abline_intercept = -1.3) 
{
  pos <- if (dodge_errorbar) 
    ggplot2::position_dodge(width = dodge_width)
  else "identity"
  obsdat <- obsdat %>% dplyr::filter(fg %in% fg_names, year == 
                                       year_to_plot, !is.na(mean), mean_n_individuals >= 20) %>% 
    dplyr::group_by(bin_midpoint) %>% dplyr::mutate(width = error_bar_width * 
                                                      dplyr::n()) %>% dplyr::ungroup()
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% dplyr::summarize(min_obs = min(bin_midpoint), 
                                                                    max_obs = max(bin_midpoint))
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% dplyr::filter(prod_model %in% 
                                                                          model_fit, fg %in% fg_names, year == year_to_plot) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% dplyr::filter(dbh >= 
                                                                             min_obs & dbh <= max_obs)
  p <- ggplot2::ggplot() + ggplot2::geom_ribbon(data = preddat, 
                                                ggplot2::aes(x = dbh, ymin = q025, ymax = q975, group = fg, 
                                                             fill = fg), alpha = 0.4) + ggplot2::geom_line(data = preddat, 
                                                                                                           ggplot2::aes(x = dbh, y = q50, group = fg, color = fg))
  if (plot_errorbar) {
    p <- p + ggplot2::geom_errorbar(data = obsdat, ggplot2::aes_string(x = "bin_midpoint", 
                                                                       ymin = error_min, ymax = error_max, group = "fg", 
                                                                       color = "fg", width = "width"), position = pos, size = error_bar_thickness)
  }
  p <- p + ggplot2::geom_line(data = preddat[preddat$fg == 
                                               "fg5", ], ggplot2::aes(x = dbh, y = q50), color = "gray") + 
    ggplot2::geom_point(data = obsdat, ggplot2::aes_string(x = "bin_midpoint", 
                                                           y = average, group = "fg", fill = "fg"), size = geom_size, 
                        color = "black", shape = 21, position = pos) + ggplot2::scale_x_log10(name = x_name, 
                                                                                              limits = x_limits, breaks = x_breaks) + ggplot2::scale_y_log10(name = y_name, 
                                                                                                                                                             limits = y_limits, breaks = y_breaks, labels = y_labels) + 
    theme_no_x() + ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")) + 
    ggplot2::scale_color_manual(values = color_names) + ggplot2::scale_fill_manual(values = fill_names) + 
    theme_plant()
  if (plot_abline) {
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", linetype = "dashed", 
                                  size = 0.75)
  }
  p
}

obs_indivprod <- obs_indivprod %>%
  filter(mean_n_individuals >= 20)
p <- plot_prod2(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = PROD,
          x_limits = c(1, 230),
          y_limits = c(0.001, 2000),
          y_breaks = c(0.001,0.1, 10, 1000),
          plot_errorbar = T,
          error_min = 'q25',
          error_max = 'q75',
          error_bar_width = 0,
          y_labels = c(0.001,0.1,10,1000),
          dodge_width = 0.05)
p0 <- p #+ annotation_custom(grob_text_a) 

p1 <- set_panel_size(p0, width=unit(8,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path,'Figures/Fig_3/Main/Fig_3a_Growth.pdf'))
grid.draw(p1)
dev.off()

geom
obs_totalprod <- obs_totalprod %>%
  filter(bin_count >= 20)
grob_text <- grobTree(textGrob("Energy Equivalence: Slope = 0", x = 0.17, y = 0.88, hjust = 0,
                               gp = gpar(col = "gray52", fontsize = 20))) 
grob_text <- grobTree(textGrob("Energy Equivalence: Slope = 0", x = 0.11, y = 0.88, hjust = 0,
                               gp = gpar(col = "gray52", fontsize = 20))) 


plot_totalprod2 <-function (year_to_plot = 1995, fg_names = c("fg1", "fg2", "fg3", 
                                                              "fg4", "fg5", "all"), model_fit_density = 1, model_fit_production = 1, 
                            x_limits, x_breaks = c(1, 3, 10, 30, 100, 300), y_limits = c(0.03, 
                                                                                         100), y_breaks = c(0.01, 0.1, 1, 10, 100, 1000), y_labels, 
                            fill_names = c("black", "#BFE046", "#267038", "#27408b", 
                                           "#87Cefa", "gray87"), color_names = c("black", "#BFE046", 
                                                                                 "#267038", "#27408b", "#87Cefa", "gray"), x_name = "Diameter (cm)", 
                            y_name = expression(paste("Production (kg ha"^-1, " cm"^-1, 
                                                      " y"^-1, ")")), geom_size = 4.5, obsdat = obs_totalprod, 
                            preddat = fitted_totalprod, plot_abline = TRUE, abline_slope = 0, 
                            abline_intercept = 2) 
{
  obsdat <- obsdat %>% dplyr::filter(fg %in% fg_names, year == 
                                       year_to_plot, bin_count > 10) %>% dplyr::filter(bin_value > 
                                                                                         0)
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% dplyr::summarize(min_obs = min(bin_midpoint), 
                                                                    max_obs = max(bin_midpoint))
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% dplyr::filter(dens_model %in% 
                                                                          model_fit_density, prod_model %in% model_fit_production, 
                                                                        fg %in% fg_names, year == year_to_plot) %>% dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                                                                                                                                     dplyr::all_vars(. > min(y_limits))) %>% dplyr::filter(dbh >= 
                                                                                                                                                                                             min_obs & dbh <= max_obs)
  p <- ggplot2::ggplot() + ggplot2::geom_ribbon(data = preddat, 
                                                aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), 
                                                alpha = 0.4) + ggplot2::geom_line(data = preddat, aes(x = dbh, 
                                                                                                      y = q50, group = fg, color = fg)) + ggplot2::geom_point(data = obsdat, 
                                                                                                                                                              aes(x = bin_midpoint, y = bin_value, group = fg, fill = fg), 
                                                                                                                                                              size = geom_size, color = "black", shape = 21) + ggplot2::scale_x_log10(name = x_name, 
                                                                                                                                                                                                                                      limits = x_limits, breaks = x_breaks) + ggplot2::scale_y_log10(name = y_name, 
                                                                                                                                                                                                                                                                                                     limits = y_limits, breaks = y_breaks, labels = y_labels, 
                                                                                                                                                                                                                                                                                                     position = "right") + ggplot2::scale_color_manual(values = color_names) + 
    ggplot2::scale_fill_manual(values = fill_names) + theme_plant() + 
    theme(aspect.ratio = 1)
  if (plot_abline) 
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", linetype = "dashed", 
                                  size = 0.75)
  p
} 
p <- plot_totalprod2(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
               model_fit_density = DENS, 
               model_fit_production = PROD,
               x_limits = c(0.9,230),
               y_limits = c(0.5, 200),
               y_breaks = c(0.1, 1, 10, 100),
               y_labels = c(0.1, 1, 10, 100),
               preddat = fitted_totalprod)
p
p0 <-  p + annotation_custom(grob_text_c) + annotation_custom(grob_text)
#p0
# p0 + geom_abline(intercept = 2, slope = 0, color ="gray72",linetype="dashed", size=.75)

p1 <- set_panel_size(p0, width=unit(12,"cm"), height=unit(14.4,"cm"))

grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/Fig_3/Main/Fig_3c_Total_Production.pdf'))
grid.draw(p1)
dev.off()


#---------------------------------------------------------------------------------------------
# ------------------------ Supplements for growth, density scaling-  ------------------------- 
#---------------------------------------------------------------------------------------------


### Add height  
plot_totalprod3 <- function (year_to_plot = 1995, fg_names = c("fg1", "fg2", "fg3", 
                                                               "fg4", "fg5", "all"), model_fit_density = 1, model_fit_production = 1, 
                             x_limits, x_breaks = c(1, 3, 10, 30, 100, 300), y_limits = c(0.03, 
                                                                                          100), y_breaks = c(0.01, 0.1, 1, 10, 100, 1000), y_labels, 
                             fill_names = c("black", "#BFE046", "#267038", "#27408b", 
                                            "#87Cefa", "gray87"), color_names = c("black", "#BFE046", 
                                                                                  "#267038", "#27408b", "#87Cefa", "gray"), x_name = "Diameter (cm)", 
                             y_name = expression(paste("Production (kg ha"^-1, " cm"^-1, 
                                                       " y"^-1, ")")), geom_size = 3, obsdat = obs_totalprod, 
                             preddat = fitted_totalprod, plot_abline = TRUE, abline_slope = 0, 
                             abline_intercept = 2) 
{
  obsdat <- obsdat %>% dplyr::filter(fg %in% fg_names, year == 
                                       year_to_plot, bin_count >= 20) %>% dplyr::filter(bin_value > 
                                                                                         0)
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% dplyr::summarize(min_obs = min(bin_midpoint), 
                                                                    max_obs = max(bin_midpoint))
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% dplyr::filter(dens_model %in% 
                                                                          model_fit_density, prod_model %in% model_fit_production, 
                                                                        fg %in% fg_names, year == year_to_plot) %>% dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                                                                                                                                     dplyr::all_vars(. > min(y_limits))) %>% dplyr::filter(dbh >= 
                                                                                                                                                                                             min_obs & dbh <= max_obs)
  p <- ggplot2::ggplot() + ggplot2::geom_ribbon(data = preddat, 
                                                aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), 
                                                alpha = 0.4) + ggplot2::geom_line(data = preddat, aes(x = dbh, 
                                                                                                      y = q50, group = fg, color = fg)) + ggplot2::geom_point(data = obsdat, 
                                                                                                                                                              aes(x = bin_midpoint, y = bin_value, group = fg, fill = fg), 
                                                                                                                                                              size = geom_size, color = "black", shape = 21) + ggplot2::scale_x_log10(name = x_name, 
                                                                                                                                                                                                                                      limits = x_limits, breaks = x_breaks) + ggplot2::scale_y_log10(name = y_name, 
                                                                                                                                                                                                                                                                                                     limits = y_limits, breaks = y_breaks, labels = y_labels, 
                                                                                                                                                                                                                                                                                                     position = "right") + ggplot2::scale_color_manual(values = color_names) + 
    ggplot2::scale_fill_manual(values = fill_names) + theme_plant() + 
    theme(aspect.ratio = 1)
  if (plot_abline) 
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", linetype = "dashed", 
                                  size = 0.75)
  p
}
p <- plot_totalprod3(year_to_plot = 1995,
                    fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                    model_fit_density = DENS, 
                    model_fit_production = PROD,
                    x_limits = c(0.5,230),
                    y_limits = c(0.5, 200),
                    y_breaks = c(0.1, 1, 10, 100),
                    y_labels = c(0.1, 1, 10, 100),
                    preddat = fitted_totalprod)
p
# Edit 15 Dec. 2019: secondary height axis now has new all-species allometry on it (though it is now inaccurate for individual species heights)
p0 <- p + scale_x_log10(name = 'Diameter (cm)',
                        breaks = c(1,3,10,30,100,300), limits = c(0.9, 230),
                        sec.axis = sec_axis(~ gMM(., a = 57.17, b = 0.7278, k = 21.57),
                                    name = "Height (m)", breaks = c(2, 3, 5, 10, 20, 40))) +
  scale_y_log10(position = "left", limits = c(0.5, 200), breaks = c(0.1, 1, 10, 100),labels = c(0.1, 1, 10, 100),
                name =expression(atop('Total Production', paste('(kg ha'^-1,' cm'^-1,' yr'^-1,')')))) +
  theme(aspect.ratio = 0.8) + theme_plant2() +geom_point(size = 1)
p0

                               
                        
p1 <- set_panel_size(p0, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/Fig_3/Height/Fig_3c_Total_Production_Height.pdf'))
grid.draw(p1)
dev.off()

# ----------------------- Range of total growth for 5 censuses 1990-2010 ------------------

minmax_prod_bycensus <- obs_totalprod %>%
  filter(bin_value > 0, !fg %in% 'unclassified') %>%
  group_by(fg, bin_midpoint) %>%
  summarize(range_min = min(bin_value), range_max = max(bin_value))

p <- ggplot(minmax_prod_bycensus, aes(x = bin_midpoint, ymin = range_min, ymax = range_max, color = fg)) +
  geom_errorbar(size = 1) +
  scale_x_log10(name = 'Diameter (cm)', breaks = c(1,3,10,30,100,300)) + 
  scale_y_log10(expression(paste('Total Production (kg ha'^-1,' cm'^-1,' yr'^-1,')')),
                breaks = 10^(-2:3), labels = as.character(10^(-2:3)), limits = c(0.1, 200)) +
  scale_color_manual(values = guild_fills2) +
  theme_plant()
p
pdf(file.path(gdrive_path,'Figures/Fig_3/Range/Total_Production_Range.pdf'))
grid.draw(p)
dev.off()

# ------------------------ Individual growth plot using diameter -------------------------

p <- plot_prod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3','fg4','fg5'),
               model_fit = PROD,
               x_limits = c(1, 230),
               y_limits = c(0.02, 1.3),
               y_breaks = c(0.03, 0.1, 0.3, 1),
               y_labels = c(0.03, 0.1, 0.3, 1),
               error_bar_width = 0,
               dodge_width = 0.05,
               obsdat = obs_indivdiamgrowth,
               preddat = fitted_indivdiamgrowth,
               plot_abline = FALSE,
               x_name = 'Diameter (cm)',
               y_name = expression(paste('Diameter growth (cm yr'^-1,')')))

p1 <- p + theme(axis.text.x = element_text(), axis.ticks.x = element_line()) + 
  labs(x = 'Diameter (cm)') + theme(plot.margin=grid::unit(c(1,1,1,1), "mm"))

grid.newpage()
grid.draw(p1)

ggsave(file.path(gdrive_path,'Figures/Fig_3/Diameter/Diam_growth.pdf'), p1, width = 7, height = 5.3)

# ------------------------   WAIC of Piecewise Models  -----------------------------------

# This section was edited by QDR, 20 Jun 2019, for the updated model fits.

fp <- file.path(gdrive_path, 'data/data_piecewisefits')
ics <- read.csv(file.path(fp, 'piecewise_ics_by_fg.csv'), stringsAsFactors = FALSE)
ics$fg <- factor(ics$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))

# Density model
base_size <- 11
p <- ggplot(ics %>% filter(criterion == 'WAIC', 
                           variable == 'density', !fg %in% 'Unclassified'), 
            aes(x = factor(dens_model), y = IC_value, ymin = IC_value - IC_stderr, ymax = IC_value + IC_stderr)) +
  facet_wrap(~ fg, labeller = label_value, scales = 'free_y') +
  geom_pointrange() +
  theme_bw(base_size = base_size, base_family = "",
           base_line_size = base_size/22, base_rect_size = base_size/11)+
  theme(strip.background = element_blank(),panel.grid = element_blank(), 
        text = element_text(size = 14),strip.text = element_text(size=12)) +
  labs(x = 'Segments in Density Function', y = "Widely Applicable Information Criterion (WAIC)")
p

p1 <- set_panel_size(p, width=unit(3.5,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_3/WAIC/WAIC_density.pdf'))
grid.draw(p1)
dev.off()

# Growth model

p <- ggplot(ics %>% filter(criterion == 'WAIC', 
                           variable == 'production', !fg %in% 'Unclassified'), 
            aes(x = factor(prod_model), y = IC_value, ymin = IC_value - IC_stderr, ymax = IC_value + IC_stderr)) +
  facet_wrap(~ fg, labeller = label_value,scales = 'free_y') +
  geom_pointrange() + #theme_facet +
  theme_bw(base_size = base_size, base_family = "",
          base_line_size = base_size/22, base_rect_size = base_size/11)+
  theme(strip.background = element_blank(),panel.grid = element_blank(), 
        text = element_text(size = 14),strip.text = element_text(size=12)) +
  labs(x = 'Segments in Growth Function',  y = "Widely Applicable Information Criterion (WAIC)")


p1 <- set_panel_size(p, width=unit(4,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_3/WAIC/WAIC_growth.pdf'))
grid.draw(p1)
dev.off()

#---------------------------------------------------------------------------------------------

#-------------------   Plot Slopes vs Size for each life history group  -------------------

slopes <- read.csv(file.path(fp, 'piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
slopes$fg <- factor(slopes$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))
slopes$variable <- factor(slopes$variable, labels = c("Density", "Individual Growth", "Production"))
colors <- c("sienna4", "yellowgreen", "springgreen4")

# Using 3 segment density and 1 segment production
p <- ggplot(slopes %>% filter((dens_model == 3 & is.na(prod_model)) | (is.na(dens_model) & prod_model == 1) | (dens_model == 3 & prod_model == 1), !fg %in% 'Unclassified'), 
            aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg, scale = 'free_y', labeller = label_value) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = "springgreen4", size = 0.3) +
  geom_hline(yintercept = -2, linetype = 'dashed', col = "sienna4", size = 0.3) +
  geom_hline(yintercept = 2, linetype = 'dashed', col = "yellowgreen", size = 0.3) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_ribbon(alpha = 0.3, size = 0.2) +
  geom_line(size = 1.25) +
  scale_x_log10(name = 'Diameter (cm)', expand = c(0,0)) +
  theme_bw() + 
  theme(strip.background = element_blank(), panel.grid = element_blank(),strip.text = element_text(size = 12),
        legend.position = 'bottom',legend.spacing.x=unit(.2, "cm"),legend.title=element_blank(),
        legend.text=element_text(size = 12),axis.title = element_text(size = 15),
        axis.text = element_text(color = "black", size = 11)) +
  labs(y = 'Slope') +
  ggtitle('Fitted Slopes: \n 3 Segment Density & 1 Segment Growth Models')+
  theme(plot.title = element_text(hjust=0.5)) #+ theme(legend.title=element_blank())
p 
pdf(file.path(gdrive_path, "Figures/Fig_3/Slopes/fitted_slopes_1_seg_production.pdf"))
p
dev.off()

# Using 3 segment density and 2 segment production
p <- ggplot(slopes %>% filter((dens_model == 3 & is.na(prod_model)) | (is.na(dens_model) & prod_model == 2) | (dens_model == 3 & prod_model == 2), !fg %in% 'Unclassified'), 
            aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg,scale = "free_y", labeller = label_value) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = "springgreen4", size = 0.3) +
  geom_hline(yintercept = -2, linetype = 'dashed', col = "sienna4",size = 0.3) +
  geom_hline(yintercept = 2, linetype = 'dashed', col = "yellowgreen",size = 0.3) +
  geom_ribbon(alpha = 0.3, size = 0.2) +
  geom_line(size = 1.25) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_x_log10(name = 'Diameter (cm)', expand = c(0,0)) + 
  theme_bw() + 
  theme(strip.background = element_blank(), panel.grid = element_blank(),strip.text = element_text(size = 12),
        legend.position = 'bottom',legend.spacing.x=unit(.2, "cm"),legend.title=element_blank(),
        legend.text=element_text(size = 12),axis.title = element_text(size = 15),
        axis.text = element_text(color = "black", size = 11)) +
  labs(y = 'Slope') + theme(axis.title = element_text(size = 15)) +
  ggtitle('Fitted Slopes: \n 3 Segment Density & 2 Segment Growth Models')+
  theme(plot.title = element_text(hjust=0.5)) +
  theme(legend.spacing.x=unit(.2, "cm"))+  theme(legend.title=element_blank())+ theme(legend.text=element_text(size = 12))
p 
pdf(file.path(gdrive_path, "Figures/Fig_3/Slopes/fitted_slopes_2_seg_production.pdf"))
p
dev.off()



#----------------------   Hex Plot of Growth Scaling  ---------------------------


fp <- file.path(gdrive_path,'data/data_forplotting') ## CHANGE PATH AS NEEDED

for (i in dir(fp, pattern = '.csv')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}


# Load raw data
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.R'))

# Process the raw data to get one single data frame with a lot of rows.

# Get only year, func group, dbh, and production (no more is needed to plot right now)
raw_prod <- do.call(rbind, map2(alltreedat, seq(1985,2010,5), function(x, y) cbind(year = y, x %>% filter(!recruit) %>% select(fg, dbh_corr, production))))

raw_prod <- raw_prod %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg', fg), 'unclassified'))


# Original grayscale
hex_scale_1 <- scale_fill_gradient(low = 'gray90', high = 'gray10', guide = FALSE)

hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000,10000), 
                                             labels = c(1,10,100,1000,10000), limits=c(1,10000))
theme_plant2 <- function () {
  ggplot2::theme(panel.grid = ggplot2::element_blank(), aspect.ratio = 0.75, 
                 axis.text = ggplot2::element_text(size = 15, color = "black"), 
                 axis.ticks.length = ggplot2::unit(0.2, "cm"), axis.title = ggplot2::element_text(size = 15), 
                 axis.title.y = ggplot2::element_text(margin = margin(r = 10)), 
                 axis.title.x = ggplot2::element_text(margin = margin(t = 10)), 
                 axis.title.x.top = ggplot2::element_text(margin = margin(b = 5)), 
                 plot.title = ggplot2::element_text(size = 15, face = "plain", 
                                                    hjust = 10), panel.border = ggplot2::element_rect(color = "black", 
                                                                                                      fill = NA, size = 1), panel.background = ggplot2::element_blank(), 
                 legend.position = "none", rect = ggplot2::element_rect(fill = "transparent"), 
                 text = ggplot2::element_text(family = "Helvetica"))
}
plot_prod_withrawdata2 <- function (year_to_plot = 1995, fg_names = c("fg1", "fg2", "fg3", 
                                                                      "fg4", "fg5", "unclassified"), 
                                    full_names = c("Fast", "Slow", "Pioneer", "Breeder", "Medium", "Unclassified"), func_names = c("power law", 
                                                                                                                                                                                    "2-segment power law"), x_limits = c(1, 300), x_breaks = c(1, 
                                                                                                                                                                                                                                               10, 100), y_limits, y_breaks, x_name = "Diameter (cm)", y_name = expression(paste("Growth (kg yr"^-1,                                                                                                                                                                                                                                                                                              ")")), line_types = c("dashed", "solid"), aspect_ratio = 0.75, 
                                    hex_scale = ggplot2::scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 
                                                                                                                                      "RdYlBu")), bias = 1)(50), trans = "log", name = "Individuals", 
                                                                              breaks = c(1, 10, 100, 1000, 10000), labels = c(1, 10, 
                                                                                                                              100, 1000, 10000), limits = c(1, 10000)), obsdat = raw_prod, 
                                    preddat = fitted_indivprod, plot_abline = TRUE, abline_slope = 2, 
                                    abline_intercept = -2.1, plot_fits = FALSE) 
{
  obsdat <- obsdat %>% dplyr::filter(fg %in% fg_names, year == 
                                       year_to_plot)
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% dplyr::summarize(min_obs = min(production), 
                                                                    max_obs = max(production))
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% dplyr::filter(fg %in% 
                                                                          fg_names, year == year_to_plot) %>% dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                                                                                                                               dplyr::all_vars(. > min(y_limits))) %>% dplyr::filter(dbh >= 
                                                                                                                                                                                       min_obs & dbh <= max_obs) %>% dplyr::mutate(prod_model = factor(prod_model, 
                                                                                                                                                                                                                                                       labels = func_names))
  labels <- setNames(full_names, fg_names)
  p <- ggplot2::ggplot() + ggplot2::geom_hex(data = obsdat, 
                                             ggplot2::aes(x = dbh_corr, y = production)) + ggplot2::facet_wrap(~fg, 
                                                                                                               ncol = 2, labeller = ggplot2::labeller(fg = labels)) + 
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, 
                           breaks = x_breaks) + ggplot2::scale_y_log10(name = y_name, 
                                                                       limits = y_limits, breaks = y_breaks) + ggplot2::scale_linetype_manual(name = "Growth fit", 
                                                                                                                                              values = line_types) + hex_scale + theme_plant2() + ggplot2::coord_fixed(ratio = aspect_ratio) + 
    ggplot2::theme(legend.position = c(0.7, 0.15), strip.background = ggplot2::element_blank(), 
                   strip.text = ggplot2::element_text(size = 15), legend.key = ggplot2::element_rect(fill = NA))
  if (plot_abline) 
    p <- p + ggplot2::geom_abline(slope = abline_slope, intercept = abline_intercept, 
                                  linetype = "dashed") + ggplot2::guides(linetype = "none")
  if (plot_fits) 
    p <- p + ggplot2::geom_line(data = preddat, ggplot2::aes(x = dbh, 
                                                             y = q50, group = prod_model, linetype = prod_model), 
                                size = 0.25)
  return(p)
}

p <- plot_prod_withrawdata2(year_to_plot = 1995,
                           fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
                           full_names = c('Fast', 'Pioneer', 'Slow', 'Breeder', 'Medium', 'Unclassified'),
                           x_limits = c(1, 316),
                           x_breaks = c(1,10, 100),
                           y_limits = c(0.001, 1000),
                           y_breaks = c(.001, .1, 10, 1000),
                           line_types = c('dashed', 'solid'),
                           hex_scale = hex_scale_log_colors,
                           plot_abline = FALSE,
                           plot_fits = TRUE)

p <- p + theme(legend.position = 'right') + 
  scale_y_log10(labels = c(0.01, 1, 100), 
                breaks = c(0.01, 1, 100), name = expression(paste('Growth (kg yr'^-1,')'))) 

p 
pdf(file.path(gdrive_path, 'Figures/Fig_3/Growth_hex/growth_hex.pdf'))
p
dev.off()



########################################################################################
# ------------------------------- Fig 4 Light Plots ------------------------------------
########################################################################################


# Section added by QDR 20 June 2019: new plots of total light scalings and total volume scalings, including fits and CIs

# Fitted values for individual light, total light, and total volume
fp_plot <- file.path(gdrive_path, 'data/data_forplotting')
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.R'))
fitted_indivlight <- read.csv(file.path(fp_plot, 'fitted_indivlight.csv'), stringsAsFactors = FALSE)
fitted_totallight <- read.csv(file.path(fp_plot, 'fitted_totallight.csv'), stringsAsFactors = FALSE)
fitted_totalvol <- read.csv(file.path(fp_plot, 'fitted_totalvol.csv'), stringsAsFactors = FALSE)
indivlightbins_fg <- read.csv(file.path(fp_plot, 'obs_indivlight.csv'), stringsAsFactors = FALSE)
totallightbins_fg <- read.csv(file.path(fp_plot, 'obs_totallight.csv'), stringsAsFactors = FALSE)
totalvolbins_fg <- read.csv(file.path(fp_plot, 'obs_totalvol.csv'), stringsAsFactors = FALSE)

#------Fig 4a------
# Plot total volume using the "totalprod" function
totalvolbins_fg <- totalvolbins_fg %>%
  filter(bin_count >= 20)
p <- plot_totalprod(year_to_plot = 1995,
                    fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                    model_fit_density = DENS, 
                    model_fit_production = PROD,
                    x_limits = c(0.9, 150),
                    y_limits = c(5, 5500),
                    y_breaks = c(1, 10, 100, 1000),
                    y_labels = c(1, 10, 100, 1000),
                    y_name = expression(paste('Total Crown Volume (m'^3, ' cm'^-1, ' ha'^-1,')')), 
                    preddat = fitted_totalvol,
                    obsdat = totalvolbins_fg, 
                    plot_abline = FALSE,
                    geom_size = 4)
p
p0 <- p + scale_y_continuous(position = "left", trans = "log10", limits = c(9, 5500),
                            name = expression(atop('Total Crown Volume',paste('(m'^3, ' cm'^-1,' ha'^-1,')')))) +
  theme(aspect.ratio = 0.75)
plot(p0)
grob_text0 <- grobTree(textGrob("a", x = 0.05, y = 0.91, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
p00 <- p0 +  annotation_custom(grob_text0) #+theme(axis.title.x = "Diameter (cm)") #element_blank())#,
                 # axis.text.x=element_blank(),
                 # axis.ticks.x=element_blank())
p00
p_tot_vol <- p00
#p1 <- set_panel_size(p, width=unit(14.3,"cm"), height=unit(14.3,"cm"))
p_tot_vol_x<- set_panel_size(p00, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p_tot_vol_x)

pdf(file.path(gdrive_path,'Figures/Fig_4/Fig_4a_Total_Crown_Vol.pdf'))
grid.draw(p_tot_vol)
dev.off()


#drop values under n = 10; are lines fitting too low

#------Fig 4b------
# Plot total light using the "totalprod" function
totallightbins_fg <- totallightbins_fg %>%
  filter(bin_count >= 20)
tot_light <- plot_totalprod(year_to_plot = 1995,
                            fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                            model_fit_density = DENS, 
                            model_fit_production = PROD,
                            x_limits = c(0.9,150),
                            y_limits = c(100, 200000),
                            geom_size = 3,
                            y_breaks = c(100, 1000, 10000, 100000),
                            y_labels = c("0.1", "1", "10", "100"),
                            y_name = expression(paste('Total Light Intercepted (kW cm'^-1,' ha'^-1,')')),
                            preddat = fitted_totallight,
                            obsdat = totallightbins_fg,
                            plot_abline = FALSE)
tot_light
tot_light1 <- tot_light + scale_y_continuous(position = "left", trans = "log10", breaks = c(100, 1000, 10000, 100000),
                            labels = c("0.1", "1", "10", "100"), limits = c(100, 400000),
                            name = expression(paste('Total Light Intercepted (W cm'^-1,' ha'^-1,')'))) +
  theme(aspect.ratio = 0.75) 
plot(tot_light1 )
grob_text <- grobTree(textGrob("Solar Equivalence", x = 0.27, y = 0.80, hjust = 0,
                               gp = gpar(col = "gray52", fontsize = 18))) 

grob_text2 <- grobTree(textGrob("a", x = 0.06, y = 0.91, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
tot_light2 <- tot_light  + scale_y_continuous(position = "left", trans = "log10", breaks = c(100, 1000, 10000, 100000),
                             labels = c("0.1", "1", "10", "100"), limits = c(100, 450000),
                             name = expression(atop('Total Light Intercepted',paste('(W m'^3, ' cm'^-1,' ha'^-1,')'))))  +
  theme(aspect.ratio = 0.75) + 
  geom_abline(intercept = log10(70000), slope = 0, color ="gray72",linetype="dashed", size=.75) +
  annotation_custom(grob_text) + annotation_custom(grob_text2)
plot(tot_light2)
p_tot_light <- tot_light2


g_tot_light <- ggplotGrob(p_tot_light)
g_tot_vol <- ggplotGrob(p_tot_vol)


g3 <- rbind(g_tot_vol, g_tot_light, size = "first")
g3$widths <- unit.pmax(g_tot_vol$widths, g_tot_light$widths)
grid.newpage()
grid.draw(g3)
ggsave(g3, height = 7.6, width = 6, filename = file.path(gdrive_path,'Figures/Fig_5/fig5_old.pdf'))


# Symmetry plot -----------------------------------------------------------


# Load parameter df to find the point at which to evaluate the fitted slopes, then load the fitted slopes

params <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/piecewise_paramci_by_fg.csv'), stringsAsFactors = FALSE)

xvalues <- params %>% 
  filter(variable == 'density', model == DENS) %>%
  group_by(fg) %>%
  summarize(mid_cutoff = (mean[parameter == 'tau_high'] + mean[parameter == 'tau_low'])/2)

# Load fitted slopes of total growth and total light.

growth_slopes <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
light_slopes <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/light_piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
growth_slopes_atmiddle <- growth_slopes %>% 
  filter(variable == 'total_production', dens_model == DENS, prod_model == PROD) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

light_slopes_atmiddle <- light_slopes %>% 
  filter(variable == 'total_incoming_light', dens_model == DENS, prod_model == PROD) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

fg_full_names <- c('Fast', 'Pioneer', 'Slow', 'Breeder', 'Medium', 'All Trees', 'Unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

allslopes <- rbind(growth_slopes_atmiddle, light_slopes_atmiddle) %>%
  ungroup %>%
  mutate(fg = factor(fg, levels = fgs, labels = fg_full_names))
grob0 <- grobTree(textGrob("b", x = 0.04, y = 0.9,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob1 <- grobTree(textGrob("Light Capture", x = 0.65, y = 0.94, hjust = 0,
                           gp = gpar(col = "gold3", fontsize = 18))) 
grob2 <- grobTree(textGrob("Production", x = 0.65, y = 0.86, hjust = 0,
                           gp = gpar(col = "darkgreen", fontsize = 18)))# fontface="italic"
grob3 <- grobTree(textGrob("Energy Equivalence", x = 0.25, y = 0.52, hjust = 0,
                           gp = gpar(col = "black", fontsize = 18))) #, fontface = "bold")))
# Plot

slopes <- ggplot(allslopes %>% filter(!fg %in% 'Unclassified'), aes(x = fg, y = q50, ymin = q025, ymax = q975, color = variable)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = .75) +
  geom_errorbar(position = position_dodge(width = 0.6), size = 1, width = 0) +
  geom_point(position = position_dodge(width = 0.6), shape = 21, size = 2.5, stroke = 1) +
  labs( x = NULL, y = 'Scaling Slope') +#x = 'Life History Guild', +
  scale_y_continuous(limits = c(-1.05, 1.3)) +#, labels = c("-1", "-0.5", "0", "0.5", "1")) +
  scale_color_manual(values = c('gold2', 'darkgreen'), labels = c('Total Light', 'Total Growth')) +
  theme_plant() + theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "italic", size = 18)) +
  annotation_custom(grob1) + annotation_custom(grob2) + annotation_custom(grob3) +annotation_custom(grob0)

#grid.newpage()
#grid.draw(slopes)
#pdf(file.path(gdrive_path, 'Figures/Fig_5/Growth Light symmetry.pdf'))
#grid.draw(slopes)
#dev.off()


g_tot_light <- ggplotGrob(p_tot_light)
g_tot_vol <- ggplotGrob(p_tot_vol)
g_slopes <- ggplotGrob(slopes)

g3 <- rbind(g_tot_light, g_slopes, size = "first")
g3$widths <- unit.pmax(g_tot_light$widths,g_slopes$widths)
grid.newpage()
grid.draw(g3)
ggsave(g3, height = 8.6, width = 6, filename = file.path(gdrive_path,'Figures/Fig_5/fig5.pdf'))

#-------------------------------------------------------------------------------------------

#----------------------Supplementals Max Growth by Light ----------------------------------------------------

### -----

# 1. Plot of maximum slope by functional group

# Remove all tree and unclassified groups
param_ci$fg <- factor(param_ci$fg ,levels = c("fg1", "fg2", "fg5", "fg3", "fg4"))
fg_labels2 <- c("Fast", "Pioneer", "Intermediate", "Slow", "Breeder")
p <- ggplot(param_ci %>% filter(fg != 'NA', year == year_to_plot, parameter %in% 'log_slope', !fg %in% c('alltree','unclassified')),
            aes(x = fg, y = q50, ymin = q025, ymax = q975)) + 
  #geom_hline(yintercept = 1, linetype = 'dotted', color = 'black', size = 1) + 
  geom_errorbar(width = 0.4) + geom_point(size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_x_discrete(name = 'Life History Strategy', labels = fg_labels2) +
  scale_y_continuous(expression(paste('Max. Growth Rate (kg yr'^-1, ' m'^-2,')')), 
                                limits = c(.6, 1.1),
                                breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2)) +
 theme_plant() + theme(aspect.ratio = 0.75)
p
pdf(file.path(gdrive_path, "Figures/Fig_4/Max_growth_by_fg.pdf"))
p
dev.off()

#---------------------------- Median binned growth by light + fg --------------------------------
p_median_panels <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels2)) +
  #facet_wrap(~ fg, labeller = label_value,scales = 'free_y') +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
              aes(x = light_area, ymin = q025, ymax = q975,group=fg, color=NA, fill = fg),
              alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), 
            aes(x = light_area, y = q50,group = fg, color=fg), size=0.5) +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q25, yend = q75), size = 0.3) +
  #geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q025, yend = q975)) +
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), 
               aes(x = x_max_q50 * 0.5, xend = x_max_q50 * 2, y = y_max_q50 * 0.5, yend = y_max_q50 * 2), color = 'brown1', size = .5) +
  geom_point(shape=21, aes(x = bin_midpoint, y = median)) +
  scale_color_manual(values = guild_fills_nb ) +
  scale_fill_manual(values = fg_colors)+
  scale_x_log10(name = title_x, breaks=c(1,10,100)) + 
  scale_y_log10(name = title_y, labels=signif, breaks=c(0.01,0.1,1,10)) +
  theme_plant() + theme_facet()

p_median_panels

pdf(file.path(gdrive_path, "Figures/Fig_4/Growth by light_Medians.pdf"))
p_median_panels
dev.off()

#---------------------------- Raw growth by light + fg --------------------------------
p_raw_panels <- ggplot(obs_light_raw %>% filter(year == year_to_plot, !fg %in% c('alltree','unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels)) +
  geom_point(shape = 21,aes(x = light_area, y = production_area, alpha = fg)) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
              aes(x = light_area, ymin = q025, ymax = q975, group=fg,color=NA,fill=fg), alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), 
            aes(x = light_area, y = q50, group=fg, color=fg), size=0.25) +
  #scale_alpha_manual(values = c(0.15, 0.15, 0.05, 0.008, 0.008)) +
  scale_alpha_manual(values = c(0.15, 0.15, 0.15, 0.15, 0.15)) +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, breaks=c(0.001,0.01,1),labels=signif) +
  scale_color_manual(values = fg_colors) +
  scale_fill_manual(values = fg_colors) +
  theme_plant() + theme_facet()

p_raw_panels
pdf(file.path(gdrive_path, "Figures/Fig_4/Growth by light/Raw.pdf"))
p_raw_panels
dev.off()


# colored hexagon plot ----------------------------------------------------


obs_light_raw$fg <- factor(obs_light_raw$fg, levels = c(paste0('fg', 1:5), 'unclassified'),
                           labels = c("Fast", "Pioneer", "Slow", "Breeder", "Medium", "Unclassified"))
pred_light_5groups$fg <- factor(pred_light_5groups$fg, levels = c(paste0('fg', 1:5)),
                                labels = c("Fast", "Pioneer", "Slow", "Breeder", "Medium"))

hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000), 
                                             labels = c(1,10,100,1000), limits=c(1,5000))

unique(obs_light_raw$fg)
p_hex_panels <- ggplot(obs_light_raw %>% filter(year == year_to_plot, !fg %in% 'Unclassified')) +
  facet_wrap(~ fg, nrow = 3) +
  geom_hex(aes(x = light_area, y = production_area)) +
  #geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
  #          aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.3) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), 
            aes(x = light_area, y = q50), color = 'black') +
  #geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q025, color = fg),  linetype = 'dotted') +
  #geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q975, color = fg),  linetype = 'dotted') +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, labels=signif) +
  hex_scale_log_colors + 
  #scale_color_manual(name = 'Functional group',values = fg_colors[1:5]) +
  #scale_fill_manual(values = fg_colors, labels = fg_labels, guide = FALSE) +
  theme_plant2() +theme(strip.text = element_text(size=14))+
  guides(color = FALSE) +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA),
        legend.position = c(0.7, 0.15),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=15))
p_hex_panels
#p <-set_panel_size(p_hex_panels, width=unit(10,"cm"), height=unit(6,"cm"))
#plot(p)
pdf(file.path(gdrive_path, "Figures/Fig_4/growth_light_heat_map.pdf"))
p_hex_panels
dev.off()




################################################################################################
# ------------------------------- Fig 5 Light Interception  ------------------------------------
################################################################################################

lightperareafakebin_fg <- read.csv(file.path(fp_plot, 'lightperareafakebin_fg.csv'), stringsAsFactors = FALSE)
lightpervolfakebin_fg <- read.csv(file.path(fp_plot, 'lightpervolfakebin_fg.csv'), stringsAsFactors = FALSE)
unscaledlightbydbhfakebin_fg <- read.csv(file.path(fp_plot, 'unscaledlightbydbhfakebin_fg.csv'), stringsAsFactors = FALSE)

# Plot: raw data ----------------------------------------------------------

exl2 <- expression(paste('Light per Crown Area (W cm'^-1, 'm'^-2, ')', sep = ''))
exl <- expression(atop('Light per Crown Area', paste('(W cm'^-1, 'm'^-2, ')')))
exv <- expression(atop('Light per Crown Volume', paste('(W cm'^-1, 'm'^-3, ')')))
exv2 <- expression(paste('Light per Crown Volume (W cm'^-1, ')', sep = ''))
exd <- 'Diameter (cm)'

fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified')
full_names = c('Fast', 'Pioneer', 'Slow', 'Breeder', 'Medium', 'Unclassified')
labels <- setNames(full_names, fg_names)

#----------------------   Fig 5a: Light per crown area by diameter -----------------------------
# All together
grob0 <- grobTree(textGrob("a", x = 0.04, y = 0.93,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
p_crown <- ggplot() + geom_point(alpha = 0.01, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownarea), color = 'chartreuse3') +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_plant() +annotation_custom(grob0) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


p_crown <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(g_crown)
pdf(file.path(gdrive_path,'Figures/Fig_5/Fig_5a.pdf'))
grid.draw(g_crown)
dev.off()

#name = expression(atop('Total Light Intercepted',paste('(W m'^3, ' cm'^-1,' ha'^-1,')'))))  
# Supp: Each group
alpha_value <- 0.05
ggplot() +
  geom_point(alpha = 0.05, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownarea), color = 'chartreuse3') +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), 
                  aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg, ncol = 2, labeller = labeller(fg = labels)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_facet2()



#----------------------   Fig 5b: Light per crown volume by diameter -----------------------------

# All together
grob1 <- grobTree(textGrob("b", x = 0.04, y = 0.93,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
p_vol <- ggplot() +
  #geom_point(alpha = 0.01, data = alltree_light_95, 
   #          aes(x = dbh_corr, y = light_received/crownvolume), color = 'chartreuse3') +
  geom_pointrange(data = lightpervolfakebin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv, breaks = c(10, 100)) +
  theme_plant() + annotation_custom(grob1)
p_vol <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(g_vol)
pdf(file.path(gdrive_path,'Figures/Fig_5/Fig_5b.pdf'))
grid.draw(g_vol)
dev.off()


grob1 <- grobTree(textGrob("b", x = 0.04, y = 0.93,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 

g_crown <- ggplotGrob(p_crown)
g_vol <- ggplotGrob(p_vol)


g5 <- rbind(g_crown, g_vol, size = "first")
g5$widths <- unit.pmax(g_crown$widths,g_vol$widths)
grid.newpage()
grid.draw(g5)
ggsave(g5, height = 8.6, width = 6, filename = file.path(gdrive_path,'Figures/Fig_4/Fig_4.pdf'))

# Plot: total unscaled light energy by dbh --------------------------------

alpha_value <- 1
hexfill2 <- scale_fill_gradient(low = 'forestgreen', high = 'navy', trans = 'log', breaks = c(1,3,10,30,100,300))
exl <- expression(atop('Intercepted Light', paste('per Individual (W)')))
exd <- 'Diameter (cm)'

labels = trans_format("log10", math_format(10^.x))


p <- ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received)) +
  geom_pointrange(data = unscaledlightbydbhfakebin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, breaks = c(1,100,10000, 1000000), limits = c(1,1000000), 
                labels = trans_format("log10", math_format(10^.x))) +
  theme_plant() + theme(legend.position = "right", legend.text = element_text(size = 15), legend.title = element_text(size = 16))+
  hex_scale_log_colors +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))# +
p
ggsave(file.path(gdrive_path,'Figures/Fig_4/Indiv_light.pdf'), plot = p)
  #geom_abline(intercept = 0.903027, slope = 2.343608)
p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p2)
pdf(file.path(gdrive_path,'Figures/Fig_4/Indiv_light.pdf'))
grid.draw(p)
dev.off()



################################################################################################
# ------------------------------- Fig 6 Life History Ratios ------------------------------------
################################################################################################


# ------------------------------- Fig 6 Relative Abundance & Production ------------------------------------

# Load 1995 bin data
load(file.path(gdrive_path, 'data/data_binned/bin_object_singleyear.RData'))

prod_ratio_diam <- breeder_stats_bydiam_byyear %>% 
  mutate(ID = 'Breeder-Pioneer') %>%
  rename(production_ratio = breeder_production_ratio, density_ratio = breeder_density_ratio) %>%
  rbind(fastslow_stats_bydiam_byyear %>% 
          mutate(ID = 'Fast-Slow') %>%
          rename(production_ratio = fastslow_production_ratio, density_ratio = fastslow_density_ratio)) %>%
  filter(year == 1995)

prod_ratio_light <- breeder_stats_bylight_byyear %>% 
  mutate(ID = 'Breeder-Pioneer') %>%
  rename(production_ratio = breeder_production_ratio, density_ratio = breeder_density_ratio) %>%
  rbind(fastslow_stats_bylight_byyear %>% 
          mutate(ID = 'Fast-Slow') %>%
          rename(production_ratio = fastslow_production_ratio, density_ratio = fastslow_density_ratio)) %>%
  filter(year == 1995)

# ---------------------------------   Fig 6a  Production by light ------------------------------------

grob1_1a <- grobTree(textGrob("1:1", x = 0.5, y = 0.52,  hjust = 0,
                             gp = gpar(col = "black", fontsize = 19))) 

prod_ratio <- prod_ratio_light   %>%
  filter(n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio, fill = ID)) +
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  theme_plant() + #annotation_custom(grob0) + #annotation_custom(grob1_1a) + 
  scale_x_log10(name = expression(paste('Light per Crown Area (W m'^-2,')')), limits=c(2,330), breaks=c(3,  30,  300)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.01,200),
                name = NULL) #name = expression("Ratio"))

prod_ratio
p1 <- set_panel_size(prod_ratio, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/Fig_6a_Production_light.pdf'))
grid.draw(p1)
dev.off()

# -------------------- ----- Supp: Density ratio by light -----------------------------

p <- prod_ratio_diam   %>%
  filter(n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio, fill = ID)) +
  geom_point(shape = 21, size = 4,  stroke = .5, color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey")) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  theme_plant() +
  scale_x_log10(name = 'Diameter (cm)', limits=c(1,150), breaks=c(1, 10, 100)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.01,200),
                name = expression("Production Ratio")) +
  theme_no_y()

p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/Sup_Production_diam.pdf'))
grid.draw(p1)
dev.off()

#----------------------------- Fig 6B Abundance by Diameter ----------------------------
grob1_1 <- grobTree(textGrob("1:1", x = 0.5, y = 0.54,  hjust = 0,
                             gp = gpar(col = "black", fontsize = 19))) 

grob_r <- grobTree(textGrob("Replacement", x = 0.55, y = 0.35,  hjust = 0,
                             gp = gpar(col = "black", fontsize = 20))) 

dens_ratio <- prod_ratio_diam %>% 
  filter(density_ratio > 0) %>%
  filter(n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  
  scale_x_log10(limits=c(1,100),breaks=c(1,10, 100), name = expression(paste('Diameter (cm)'))) + 
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.01,200),
                name = expression("Ratio")) + 
  theme_plant() + #annotation_custom(grob1) + annotation_custom(grob1_1) + 
  #annotation_custom(grob_r) + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
dens_ratio 

p1 <- set_panel_size(dens_ratio , width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/Fig_6b_Density.pdf'))
grid.draw(p1)
dev.off()




# -------------------- ----- Supp: Production ratio by diameter  -----------------------------

p <- prod_ratio_diam  %>%
  filter(production_ratio > 0) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  theme_plant() +
  scale_x_log10(name = expression(paste('Diameter (cm)')), 
                limits=c(1,150), breaks=c(1,3,10,30,100,300)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100), limits=c(0.01,200),
                name = expression("Production Ratio")) 
p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/Sup Ratios/Production_diam.pdf'))
grid.draw(p1)
dev.off()



grob1 <- grobTree(textGrob("b", x = 0.04, y = 0.93,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob1_1 <- grobTree(textGrob("1:1", x = 0.5, y = 0.5,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25))) 

g_dens  <- ggplotGrob(dens_ratio )
g_prod <- ggplotGrob(prod_ratio)
g6 <- cbind(g_prod, g_dens , size = "first")
g6$heights <- unit.pmax(g_crown$heights, g_vol$heights)
grid.newpage()
grid.draw(g6)
ggsave(g6, height = 8, width = 11, filename = file.path(gdrive_path,'Figures/Fig_6/Fig_6ab.pdf'))

# -------------------- ------   Supp: PCA score   -----------------------------
## by Diameter

error_bar_width <- 0.15

# Combine
fastslowscore_bin_bydiam<- fastslowscore_bin_bydiam_byyear %>%
  left_join(fastslow_stats_bydiam_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Fast-Slow")
breederscore_bin_bydiam <- breederscore_bin_bydiam_byyear %>%
  left_join(breeder_stats_bydiam_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = 'Breeder-Pioneer')
score_bin_bydiam <- rbind(fastslowscore_bin_bydiam, breederscore_bin_bydiam)
grob_a <- grobTree(textGrob("a", x = 0.04, y = 0.93,  hjust = 0,
                            gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 

PCA_diam <- score_bin_bydiam %>%
  filter(year == 1995, n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  geom_errorbar(width = error_bar_width) +theme_plant() +
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  scale_x_log10(name = 'Diameter (cm)', limits = c(.8,100)) + 
  annotation_custom(grob_a) +
  scale_y_continuous(limits=c(-1.5,1.25),breaks=c(-1,0,1),name = "PCA Score") #+
  #theme(axis.text.y=element_blank()) + theme(axis.ticks.y = element_blank())
PCA_diam
p1 <- set_panel_size(PCA_diam, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/PCA_Scores/PCA_diam.pdf'))
grid.draw(p1)
dev.off()


#-------------------------------  PCA score by light ----------------------------------

# Combine
fastslowscore_bin_bylight <- fastslowscore_bin_bylight_byyear %>%
  left_join(fastslow_stats_bylight_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Fast-Slow")
breederscore_bin_bylight<- breederscore_bin_bylight_byyear  %>%
  left_join(breeder_stats_bylight_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Breeder-Pioneer")
PCA_score_by_light <- rbind(fastslowscore_bin_bylight,breederscore_bin_bylight)
grob_b <- grobTree(textGrob("b", x = 0.04, y = 0.93,  hjust = 0,
                            gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob_p <- grobTree(textGrob("Pioneer:Breeder", x = 0.13, y = 0.93,  hjust = 0,
                            gp = gpar(col = "black", fontface = "italic", fontface = "bold", fontsize = 20))) 

grob_f <- grobTree(textGrob("Fast:Slow", x = 0.13, y = 0.83,  hjust = 0,
                            gp = gpar(col = "gray43",  fontface = "italic",fontsize = 20))) 
error_bar_width <- .15
PCA_light <- PCA_score_by_light %>%
  filter(year == 1995, n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, fill = ID)) +
  geom_errorbar(width = error_bar_width) +theme_plant() +
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+ annotation_custom(grob_a) +
  annotation_custom(grob_p) + annotation_custom(grob_f) +
  scale_x_log10(limits=c(1.8,450), breaks = c(3, 30, 300),
                name = expression(paste('Light per Crown Area (W m'^-2,')'))) + 
  scale_y_continuous(limits=c(-1.5,1.25),breaks=c(-1,0,1),name = 'PCA Score') 
PCA_light




g_light  <- ggplotGrob(PCA_light)
g_diam <- ggplotGrob(PCA_diam)
g6 <- rbind(g_light, g_diam , size = "first")
g6$widths <- unit.pmax(g_light$widths,g_diam $widths)
grid.newpage()
grid.draw(g6)
ggsave(g6, height = 8, width = 11, filename = file.path(gdrive_path,'Figures/Fig_6/PCA.pdf'))

