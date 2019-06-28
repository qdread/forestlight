

### Plotting of main and supplemental figures

# Set paths to google drive forest light folder, and github forest light folder
# All paths below will be relative to these paths.

#Q Dog
gdrive_path <- '~/google_drive/ForestLight'
github_path <- '~/Documents/GitHub/forestlight'

#Grady
gdrive_path <- '/Users/jgradym/Google Drive/ForestLight'
github_path <- '/Users/jgradym/Documents/GitHub/forestlight'

library(tidyverse)
library(egg)
library(scales)

# Some Plotting Code


theme_plant <- theme(panel.grid = element_blank(), #for Total Production
                      aspect.ratio = .75,
                      axis.text = element_text(size = 19, color = "black"), 
                      axis.ticks.length=unit(0.2,"cm"),
                      axis.title = element_text(size = 19),
                      axis.title.y = element_text(margin = margin(r = 10)),
                      axis.title.x = element_text(margin = margin(t = 10)),
                      axis.title.x.top = element_text(margin = margin(b = 5)),
                      plot.title = element_text(size = 19, face = "plain", hjust = 10),
                      panel.border = element_rect(color = "black", fill=NA,  size=1),
                      panel.background = element_blank(),
                      legend.position = "none",
                      rect = element_rect(fill = "transparent"),
                      text = element_text(family = 'Helvetica')) 

theme_facet <- theme(strip.background = element_rect(fill=NA),
      panel.border = element_rect(color = "black", fill=NA,  size=.75),legend.position = 'none',
      panel.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text = element_text(size = 15, color = "black"), 
      axis.ticks.length=unit(0.2,"cm"),
      axis.title = element_text(size = 15)) 

theme_facet2 <- theme(strip.background = element_rect(fill=NA),
                     panel.border = element_rect(color = "black", fill=NA,  size=.75),legend.position = 'none',
                     panel.background = element_blank(),
                     strip.text.x = element_blank(),
                     axis.text = element_text(size = 15, color = "black"), 
                     axis.ticks.length=unit(0.2,"cm"),
                     axis.title = element_text(size = 15)) 

theme_no_y <- theme(axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank())

theme_no_x <- theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
################################################################################################
# ------------------------------ Fig 1: hand drawn schematics ---------------------------------
################################################################################################



################################################################################################
# ------------------------------ Fig 2: PCA Axes ---------------------------------
################################################################################################

# Classification description: 

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

guild_fills <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "ivory")#RColorBrewer::brewer.pal(5, 'Set1')
guild_colors <- c("black", "#99AF3C", "#1D5128", "#1e3160", "#6BA3BF", "#D6D6D6")
guild_fills_nb <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "ivory")
guild_fills_nb0 <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")
guild_colors_nb0 <- c("#99AF3C", "#1D5128", "#1e3160", "#6BA3BF", "#D6D6D6")
guild_colors_nb <- c("#3B4403", "#02330A", "#031a49", "#02394F", "#595A5B")
fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('Fast','LL Pioneer', 'Slow', 'SL Breeder', 'Medium')

p <- ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, fill = factor(fg5))) +
  geom_point(shape = 21, size = 4.5, color = "black") + theme_plant +
  labs(x = 'Slow to Fast', y = 'Breeders to Pioneers') +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_color_manual(values = guild_colors_nb, labels = fg_labels, name = 'functional group')+
  scale_fill_manual(values = guild_fills_nb)
p
pdf(file.path(gdrive_path, 'Figures/Fig_1/Life_Histories.pdf'))
p
dev.off()


################################################################################################
# ------------------------Fig 3: Plotting for Piecewise Scaling --------------------------------
################################################################################################
# Define plotting functions -----

geom_size = 4
# Plot single model fit with multiple functional groups for density
plot_dens <- function(year_to_plot = 1995,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                      model_fit = 'pareto',
                      x_limits,
                      x_breaks = c(1, 3, 10, 30, 100,300),
                      y_limits,
                      y_breaks,
                      y_labels,
                      fill_names = c("black","#BFE046", "#267038", "#27408b", "#87Cefa", "ivory"),
                      fill_names0 = c("black","#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),
                      x_name = 'Diameter (cm)',
                      y_name = expression(paste('Density (n ha'^-1,'cm'^-1,')')),
                      obsdat = obs_dens,
                      preddat = pred_dens
) {
  
  require(dplyr)
  require(ggplot2)
  
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot, bin_count > 10) %>%
    filter(bin_value > 0)
  
  # Get minimum and maximum observed bin value for each group to be plotted
  # Delete points on the predicted line that are outside of this range
  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(dens_model %in% model_fit, fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs) 
  
  ggplot() +
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
    geom_abline(intercept=4, slope = -2, color ="gray72",linetype="dashed", size=.75)+   
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    geom_line(data = preddat[preddat$fg == "fg5",], aes(x = dbh, y = q50), color = "gray")+ # white circles get gray line
    geom_ribbon(data = preddat[preddat$fg == "fg5",], aes(x = dbh, ymin = q025, ymax = q975), fill = "gray", alpha = 0.4) +
    geom_point(data = obsdat, aes(x = bin_midpoint, y = bin_value, group = fg, fill=fg), size = geom_size, shape=21,color="black") +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks,labels = y_labels) +
    scale_color_manual(values = fill_names0) +theme_plant + 
    scale_fill_manual(values = fill_names) 
  
  
}

# Plot single model fit with multiple functional groups for individual production
plot_prod <- function(year_to_plot = 1995,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                      model_fit = 'powerlaw',
                      x_limits,
                      x_breaks = c(1, 3, 10, 30, 100,300),
                      y_limits,
                      y_labels,
                      y_breaks,
                      fill_names = c("#BFE046", "#267038", "#27408b", "#87Cefa", "ivory"),
                      fill_names0 = c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),
                      x_name = 'Diameter (cm)',
                      y_name = expression(paste('Growth (kg y'^-1,')')),
                      average = 'mean',
                      error_quantiles = c('ci_min', 'ci_max'),
                      error_bar_width = 0.03,
                      dodge_width = 0.03,
                      dodge_errorbar = TRUE,
                      obsdat = obs_indivprod,
                      preddat = fitted_indivprod
) {
  
  require(dplyr)
  require(ggplot2)
  
  pos <- if (dodge_errorbar) position_dodge(width = dodge_width) else 'identity'
  
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot, !is.na(mean), mean_n_individuals > 10) %>%
    group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup
  
  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(prod_model %in% model_fit, fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs)
  
  ggplot() +
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
    geom_abline(intercept= -1.6, slope = 2, color ="gray72",linetype="dashed", size=.75)+   
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    #geom_errorbar(data = obsdat, aes_string(x = 'bin_midpoint', ymin = error_quantiles[1], ymax = error_quantiles[2], group = 'fg', color = 'fg', width = 'width'), position = pos) +
    geom_line(data = preddat[preddat$fg == "fg5",], aes(x = dbh, y = q50), color = "gray")+ # white circles get gray line
    geom_ribbon(data = preddat[preddat$fg == "fg5",], 
                aes(x = dbh, ymin = q025, ymax = q975), 
                fill = "gray", alpha = 0.4) +
    geom_point(data = obsdat, aes_string(x = 'bin_midpoint', y = average, group = 'fg', fill = 'fg'),
               size = geom_size,color="black",shape=21, position = pos) +
    scale_x_log10()+#name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) +
    theme_no_x +theme(rect = element_rect(fill = "transparent"))+ # all rectangles
    scale_color_manual(values = fill_names0) +
    scale_fill_manual(values = fill_names) + theme_plant
  
  
}

# Plot single model fit with multiple functional groups for total production
# Specify both density and production type

plot_totalprod <- function(year_to_plot = 1995,
                           fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                           model_fit_density = 'pareto',
                           model_fit_production = 'powerlaw',
                           x_limits,
                           x_breaks = c(1, 3, 10, 30, 100,300),
                           y_limits = c(0.03,100),
                           y_breaks = c(0.01,0.1, 1, 10,100, 1000),
                           y_labels,
                           fill_names = c("black","#BFE046", "#267038", "#27408b", "#87Cefa", "ivory"),
                           fill_names0 = c("black","#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),
                           x_name = 'Diameter (cm)',
                           y_name = expression(paste('Total Production (kg ha'^-1,' cm'^-1,' yr'^-1,')')),
                           obsdat = obs_totalprod,
                           preddat = fitted_totalprod,
                           plot_abline = TRUE
) {
  
  require(dplyr)
  require(ggplot2)
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot, bin_count > 10) %>%
    filter(bin_value > 0)
  
  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(dens_model %in% model_fit_density, prod_model %in% model_fit_production, fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs)
  
  p <- ggplot() +
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
       
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    #geom_line(data = preddat[preddat$fg == "fg5",], aes(x = dbh, y = q50), color = "gray")+ # white circles get gray line
    geom_ribbon(data = preddat[preddat$fg == "fg5",], aes(x = dbh, ymin = q025, ymax = q975), fill = "gray", alpha = 0.4) +
    geom_point(data = obsdat, aes(x = bin_midpoint, y = bin_value,group = fg, fill=fg), size = geom_size, color = "black",shape=21) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels, position = "right") +
    scale_color_manual(values = fill_names0) + theme_plant + theme(aspect.ratio = 1) +
    scale_fill_manual(values = fill_names) 
  
  if (plot_abline) p <- p + geom_abline(intercept= 2, slope = 0, color ="gray72",linetype="dashed", size=.75)
  p
  
}  

# Plot of slopes in different segments by different functional groups.


fp <- file.path(gdrive_path, 'data/data_piecewisefits')
fpfig <- file.path(gdrive_path, 'figs/piecewiseplots_27sep2018')


# Plot the fitted values on top of the observed histograms.

# Read observed data
fp_obs <- file.path(gdrive_path, 'data/data_forplotting_aug2018')

for (i in dir(fp_obs, pattern = 'obs_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp_obs, i), stringsAsFactors = FALSE))
}

# Read modeled data (CIs)

for (i in dir(fp, pattern = 'pred_|fitted_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}

#source(file.path(github_path, 'stan/piecewise_workflow/plottingfunctionspiecewise.r'))
# Create plots.
#Model fit 1 = pareto, 1 segment
#Model Fit 2  = 2 segments, etc


p <- plot_dens(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = 3,
          x_limits = c(1, 260),
          y_limits = c(0.003, 3000),
          y_labels = c(0.001, 0.1, 10,1000),
          y_breaks = c(0.001, 0.1,  10, 1000))
p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path,'Figures/Fig_3/Main/Fig_3b_Density.pdf'))
plot(p1)
dev.off()

# Specify dodging with a certain width of error bar
# Model fit 1 = power law
# Model fit 2 = power law exp
p <- plot_prod(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = 2,
          x_limits = c(1, 280),
          y_limits = c(0.001, 2000),
          y_breaks = c(0.001,0.1, 10, 1000),
          y_labels = c(0.001,0.1,10,1000),
          error_bar_width = 0.01,
          dodge_width = 0.05)
p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path,'Figures/Fig_3/Main/Fig_3a_Growth.pdf'))
plot(p1)
dev.off()

p <- plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
               model_fit_density = 3, 
               model_fit_production = 2,
               x_limits = c(0.9,250),
               y_limits = c(0.03, 200),
               y_breaks = c(0.1, 1, 10, 100),
               y_labels = c(0.1, 1, 10, 100),
               preddat = fitted_totalprod)
p
p <- p + geom_abline(intercept= 2, slope = 0, color ="gray72",linetype="dashed", size=.75)
p1 <- set_panel_size(p, width=unit(14.3,"cm"), height=unit(14.3,"cm"))
#p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)

pdf(file.path(gdrive_path,'Figures/Fig_3/Main/Fig_3c_Total_Production.pdf'))
plot(p1)
dev.off()


#---------------------------------------------------------------------------------------------
# ------------------------ Supplements for growth, density scaling-  ------------------------- 
#---------------------------------------------------------------------------------------------



# ------------------------   WOOIC of Piecewise Models  -----------------------------------

# This section was edited by QDR, 20 Jun 2019, for the updated model fits.
# We are now using WAIC as the information criterion because the LOOIC unfortunately returned errors for some groups. Results are the same though.

ics <- read.csv(file.path(fp, 'newpiecewise_ics_by_fg.csv'), stringsAsFactors = FALSE)
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
        text = element_text(size = 14), strip.text = element_text(size=12)) +
  labs(x = 'Segments in Density Function', y = "Widely Applicable Information Criterion (WAIC)")
p
ggsave(file.path(gdrive_path, 'Figures/Fig_3/WOOIC/WAIC_density.pdf'), p)
# pdf(file.path(gdrive_path, "Figures/LOOIC/LOOIC_density.pdf"))
# p
# dev.off()
#ggsave(file.path(fpfig, 'density_model_information_criteria.pdf'), height = 6, width = 9)

# Production model

p <- ggplot(ics %>% filter(criterion == 'WAIC', 
                           variable == 'production', !fg %in% 'Unclassified'), 
            aes(x = factor(prod_model), y = IC_value, ymin = IC_value - IC_stderr, ymax = IC_value + IC_stderr)) +
  facet_wrap(~ fg, labeller = label_value,scales = 'free_y') +
  geom_pointrange() +
  theme_bw(base_size = base_size, base_family = "",
           base_line_size = base_size/22, base_rect_size = base_size/11)+
  theme(strip.background = element_blank(),panel.grid = element_blank(), 
        text = element_text(size = 14),strip.text = element_text(size=12)) +
  labs(x = 'Segments in Production Function',  y = "Widely Applicable Information Criterion (WAIC)")
p
ggsave(file.path(gdrive_path, 'Figures/Fig_3/WOOIC/WAIC_production.pdf'), p)
# pdf(file.path(gdrive_path, "Figures/LOOIC/LOOIC_production.pdf"))
# p
# dev.off()
#ggsave(file.path(fpfig, 'production_model_information_criteria.pdf'), height = 6, width = 9)
#---------------------------------------------------------------------------------------------

#-------------------   Plot Slopes vs Size for each life history group  -------------------

slopes <- read.csv(file.path(fp, 'newpiecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
slopes$fg <- factor(slopes$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))
slopes$variable <- factor(slopes$variable, labels = c("Density", "Individual Growth", "Production"))
colors <- c("sienna4", "yellowgreen", "springgreen4")

# Using 3 segment density and 1 segment production
p <- ggplot(slopes %>% filter((dens_model == 3 & is.na(prod_model)) | (is.na(dens_model) & prod_model == 1) | (dens_model == 3 & prod_model == 1), !fg %in% 'Unclassified'), 
            aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg,labeller = label_value) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = "springgreen4", size = 0.3) +
  geom_hline(yintercept = -2, linetype = 'dashed', col = "sienna4", size = 0.3) +
  geom_hline(yintercept = 2, linetype = 'dashed', col = "yellowgreen", size = 0.3) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_ribbon(alpha = 0.4, size = 0.3) +
  geom_line() +
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
  geom_ribbon(alpha = 0.4, size = 0.3) +
  geom_line() +
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
#ggsave(file.path(fpfig, 'fitted_slopes_3partdensity_2partproduction.pdf'), height = 6, width = 9)



#----------------------   Hex Plot of Growth Scaling  ---------------------------


fp <- '~/google_drive/ForestLight/data/data_forplotting_aug2018' ## CHANGE PATH AS NEEDED
fp <- '~/Google Drive/ForestLight/data/data_forplotting_aug2018' ## CHANGE PATH AS NEEDED

for (i in dir(fp, pattern = '.csv')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}


# Load raw data
load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')
load('~/Google Drive/ForestLight/data/rawdataobj_alternativecluster.r')


# Process the raw data to get one single data frame with a lot of rows.
library(dplyr)
library(purrr)

# Get only year, func group, dbh, and production (no more is needed to plot right now)
raw_prod <- do.call(rbind, map2(alltreedat, seq(1985,2010,5), function(x, y) cbind(year = y, x %>% select(fg, dbh_corr, production))))

raw_prod <- raw_prod %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg', fg), 'unclassified'))

# Function to plot production with raw data -------

plot_prod <- function(year_to_plot = 1995,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
                      full_names = c('Fast', 'Slow', 'Pioneer', 'Breeder', 'Medium', 'Unclassified'),
                      func_names = c('power law', 'Power Law\ntimes Exponential'),
                      x_limits = c(1, 300),
                      x_breaks = c(1, 10, 100),
                      y_limits,
                      y_breaks,
                      x_name = 'Diameter (cm)',
                      y_name = expression(paste('Growth (kg yr'^-1,')')), line_types = c('dashed', 'solid'),
                      aspect_ratio = 0.75,
                      hex_scale = hex_scale_log_colors,
                      obsdat = raw_prod,
                      preddat = fitted_indivprod
) {
  
  require(dplyr)
  require(ggplot2)
  #require(pracma)
  
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot)
  
  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(production), max_obs = max(production))
  
  # Add more identical rows of the f. groups with less data so that the fill scales appear the same on all plots.
  #multiples <- Reduce(Lcm, table(obsdat$fg)) / table(obsdat$fg)
  #obsdat <- obsdat[rep(1:nrow(obsdat), ceiling(multiples/(min(multiples)/3))[obsdat$fg]), ]
  
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(dens_model == 'pareto', fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs) %>%
    mutate(prod_model = factor(prod_model, labels = func_names))
  
  labels <- setNames(full_names, fg_names)
  
  ggplot() +
    geom_hex(data = obsdat, aes(x = dbh_corr, y = production)) +
    #geom_line(data = preddat, aes(x = dbh, y = q50, group = prod_model, linetype = prod_model), size=0.25) +
    geom_abline(slope = 2, intercept = -2.1, linetype = "dashed")+
    facet_wrap(~ fg, ncol = 2, labeller = labeller(fg = labels)) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
  scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels=signif) +
    #scale_linetype_manual(values = line_types, name = 'Functional form') +
    scale_linetype_manual(values = line_types) +
    hex_scale + theme_plant+
    coord_fixed(ratio = aspect_ratio) + guides(linetype = 'none')+
    theme(legend.position = c(0.7, 0.15), strip.background = element_blank(),
          strip.text = element_text(size=12))#, legend.text = element_blank())
  
  
}
# Original grayscale
hex_scale_1 <- scale_fill_gradient(low = 'gray90', high = 'gray10', guide = FALSE)

# Biased grayscale to emphasize the variation in the hexagons with less data.
# Bias >1 does this.
#hex_scale_2 <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=10)(50))

# Biased color scale of red yellow and blue to do the same, using a brewer palette
# the bias is <1 this time because I had to reverse the scale.
#hex_scale_3 <- scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(9,'RdYlBu'), bias=0.3)(50)), guide = FALSE)

# Biased and customized color scale
#hex_scale_4 <- scale_fill_gradientn(colours = colorRampPalette(c('forestgreen', 'goldenrod', 'indianred'), bias = 3)(50), guide = FALSE)

# Edit the hex scale argument to draw this with other color scales.
#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=1)(50), trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))
hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000,10000), 
                                             labels = c(1,10,100,1000,10000), limits=c(1,10000))


p <- plot_prod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
               full_names = c('Fast', 'Pioneer', 'Slow', 'Breeder', 'Medium', 'Unclassified'),
               x_limits = c(1, 316),
               x_breaks = c(1,10, 100),
               y_limits = c(0.001, 1000),
               y_breaks = c(.001, .1, 10, 1000),
               line_types = c('dashed', 'solid'),
               hex_scale = hex_scale_log_colors)
p 
pdf(file.path(gdrive_path, 'Figures/Growth_hex/growth_hex.pdf'))
p
dev.off()
#pdf(file.path(gdrive_path, "Plots_J/New/Supplementals/raw_growth_heat.pdf"))
#p
#dev.off()       



########################################################################################
# ------------------------------- Fig 4 Light Plots ------------------------------------
########################################################################################


# Section added by QDR 20 June 2019: new plots of total light scalings and total volume scalings, including fits and CIs

# Fitted values for individual light, total light, and total volume
fp <- file.path(gdrive_path, 'data/data_piecewisefits')
fitted_indivlight <- read.csv(file.path(fp, 'fitted_indivlight.csv'), stringsAsFactors = FALSE)
fitted_totallight <- read.csv(file.path(fp, 'fitted_totallight.csv'), stringsAsFactors = FALSE)
fitted_totalvol <- read.csv(file.path(fp, 'fitted_totalvol.csv'), stringsAsFactors = FALSE)

# Observed (binned) values for 1995 for individual light, total light, and total volume
#load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))
source(file.path(github_path, 'code/allfunctions27july.r'))
fp_obs <- file.path(gdrive_path, 'data/data_forplotting_aug2018')
binedgedata <- read.csv(file.path(fp_obs,'obs_dens.csv'),stringsAsFactors = FALSE) %>% filter(fg == 'all', year == 1995) 
area_core <- 42.84

totallightbins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = light_received, edges = binedgedata))
totalvolbins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = crownvolume, edges = binedgedata))

totallightbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$light_received, edges = binedgedata)) %>%
  ungroup %>%
  rbind(data.frame(fg = 'all', totallightbins_all)) %>%
  mutate(bin_value = bin_value / area_core, year = 1995)
totalvolbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$crownvolume, edges = binedgedata)) %>%
  ungroup %>%
  rbind(data.frame(fg = 'all', totalvolbins_all)) %>%
  mutate(bin_value = bin_value / area_core, year = 1995)

indivlightbins_all <- alltree_light_95 %>%
  mutate(indivlight_bin = cut(dbh_corr, breaks = c(binedgedata$bin_min[1], binedgedata$bin_max), labels = binedgedata$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(indivlight_bin) %>%
  do(c(n = nrow(.), quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))
indivlightbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  mutate(indivlight_bin =cut(dbh_corr, breaks = c(binedgedata$bin_min[1], binedgedata$bin_max), labels = binedgedata$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, indivlight_bin) %>%
  do(c(n = nrow(.), quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))
indivlightbins_fg <- data.frame(fg = 'all', indivlightbins_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(indivlightbins_fg)) %>%
  mutate(indivlight_bin = as.numeric(as.character(indivlight_bin))) %>%
  rename(bin_midpoint = indivlight_bin)

#------Fig 4a------
# Plot total volume using the "totalprod" function
p <- plot_totalprod(year_to_plot = 1995,
                    fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                    model_fit_density = 3, 
                    model_fit_production = 2,
                    x_limits = c(0.9, 150),
                    y_limits = c(5, 2500),
                    y_breaks = c(1, 10, 100, 1000),
                    y_labels = c(1, 10, 100, 1000),
                    y_name = expression(paste('Total Crown Volume (m'^3, ' cm'^-1, ' ha'^-1,')')), 
                    preddat = fitted_totalvol %>% mutate(prod_model = 2),
                    obsdat = totalvolbins_fg, 
                    plot_abline = FALSE)
p
p0 <- p + scale_y_continuous(position = "left", trans = "log10", 
                            name = expression(atop('Total Crown Volume',paste('(m'^3, ' cm'^-1,' ha'^-1,')')))) +
  theme(aspect.ratio = 0.75)
plot(p0)
p00 <- p0 + theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())
p00
p1 <- set_panel_size(p, width=unit(14.3,"cm"), height=unit(14.3,"cm"))
p1 <- set_panel_size(p0, width=unit(10.25,"cm"), height=unit(7,"cm"))
p1 <- set_panel_size(p00, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path,'Figures/Fig_4/Fig_4a_Total_Crown_Vol.pdf'))
plot(p1)
dev.off()


#drop values under n = 10; are lines fitting too low

#------Fig 4b------
# Plot total light using the "totalprod" function
p <- plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
               model_fit_density = 3, 
               model_fit_production = 2,
               x_limits = c(0.9,150),
               y_limits = c(100, 100000),
               y_breaks = c(100, 1000, 10000, 100000),
               y_labels = c("0.1", "1", "10", "100"),
               y_name = expression(paste('Total Light Intercepted (kW cm'^-1,' ha'^-1,')')),
               preddat = fitted_totallight,
               obsdat = totallightbins_fg,
               plot_abline = FALSE)
p
p0 <- p + scale_y_continuous(position = "left", trans = "log10", breaks = c(100, 1000, 10000, 100000),
                            labels = c("0.1", "1", "10", "100"), #limits = c(100, 100000),
                            name = expression(paste('Total Light Intercepted (W cm'^-1,' ha'^-1,')'))) +
  theme(aspect.ratio = 0.75)
plot(p0)

p1 <- p + scale_y_continuous(position = "left", trans = "log10", breaks = c(100, 1000, 10000, 100000),
                             labels = c("0.1", "1", "10", "100"), #limits = c(100, 100000),
                             name = expression(atop('Total Light Intercepted',paste('(W m'^3, ' cm'^-1,' ha'^-1,')'))))  +
  theme(aspect.ratio = 0.75)
plot(p1)


p3 <- set_panel_size(p, width=unit(14.3,"cm"), height=unit(14.3,"cm"))
plot(p3)
p4 <- set_panel_size(p1, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p4)
pdf(file.path(gdrive_path,'Figures/Fig_4/Fig_4b_Total_light.pdf'))
plot(p4)
dev.off()

# Plot individual light using the modified "prod" function
source(file.path(github_path, 'code/plotting/plot_prod_fixed.R'))

p <- plot_prod_fixed(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = 2,
          x_limits = c(1, 150),
          y_limits = c(10, 1e6),
          y_breaks = c(10, 1000,1e5),
          y_labels =  c("0.01", "1", "100"),
          x_name = 'Diameter (cm)',
          y_name = 'Individual Light Intercepted (kW)',
          error_bar_width = 0.01,
          dodge_width = 0.05,
          preddat = fitted_indivlight,
          obsdat = indivlightbins_fg %>% mutate(year = 1995, mean = q50) %>% rename(mean_n_individuals = bin_count))
p
p1 <- p + scale_y_continuous(position = "left", trans = "log10", breaks = c(10, 1000,1e5),
                             labels = c("0.01", "1", "100"), #limits = c(100, 100000),
                             name = expression(atop('Individual Light',paste('Intercepted (kW)'))))  +
  theme(aspect.ratio = 0.75, axis.ticks.x = element_line())
plot(p1)  
p2 <- set_panel_size(p1, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p2)
pdf(file.path(gdrive_path,'Figures/Fig_4/Indiv_light.pdf'))
plot(p1)
dev.off()


#----------------------Supplementals Max Growth by Light ----------------------------------------------------
year_to_plot <- 1995 ### CHANGE THIS IF YOU WANT TO PLOT 1990

# Load data ----

fp <- file.path(gdrive_path, 'data/data_forplotting_light_june2018')

# New File path needed
obs_light_binned <- read.csv(file.path(fp, 'obs_light_binned.csv'), stringsAsFactors = FALSE)
obs_light_raw <- read.csv(file.path(fp, 'obs_light_raw.csv'), stringsAsFactors = FALSE)
pred_light <- read.csv(file.path(fp, 'pred_light.csv'), stringsAsFactors = FALSE)
param_ci <- read.csv(file.path(fp, 'lightbyarea_paramci_by_fg.csv'), stringsAsFactors = FALSE)

# Get rid of the predicted points that are outside the limits of the observed data for each FG
obs_limits <- obs_light_binned %>%
  group_by(fg, year) %>%
  summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))

pred_light <- pred_light %>%
  left_join(obs_limits) %>%
  filter(light_area >= min_obs & light_area <= max_obs)

pred_light_5groups <- pred_light %>% filter(!fg %in% c('alltree','unclassified'))

# Do some additional computation to correct the error bar width for the number of groups in each bin
obs_light_binned <- obs_light_binned %>%
  group_by(bin_midpoint, year) %>% mutate(width = sum(c('fg1','fg2','fg3','fg4','fg5') %in% fg)) %>% ungroup


# Axis titles
title_x <- expression(paste('Light per Crown Area (W m'^-2,')',sep=''))
title_y <- expression(paste('Growth (kg yr'^-1, ' m'^-2,')', sep=''))

# Colors
colors <- c('black',"#BFE046","#267038" , "#27408b", "#87Cefa", "ivory" ) # #BFE046= light green,#267038=dark green,27408b=dark blue,"#87Cefa" = light blue,    
year_to_plot = 1995
fg_colors <- c(fg1 = "#BFE046", fg2 = "#27408b" , fg3 = "#267038", fg4 =  "#87Cefa", fg5 = "gray70",  alltree = 'black') #unclassified = 'brown',
fg_colors2 <- c(fg1 = "#BFE046", fg2 = "#27408b" , fg3 = "#267038", fg4 =  "#87Cefa", fg5 = "ivory",  alltree = 'black') #unclassified = 'brown',

### -----

# 1. Plot of maximum slope by functional group

# Remove all tree and unclassified groups
param_ci$fg <- factor(param_ci$fg ,levels = c("fg1", "fg2", "fg5", "fg3", "fg4"))
fg_labels2 <- c("Fast", "LL Pioneer", "Intermediate", "Slow", "SL Breeder")
p <- ggplot(param_ci %>% filter(fg != 'NA', year == year_to_plot, parameter %in% 'log_slope', !fg %in% c('alltree','unclassified')),
            aes(x = fg, y = q50, ymin = q25, ymax = q75)) +theme_facet +
  geom_hline(yintercept = 1, linetype = 'dotted', color = 'black', size = 1) + 
  geom_errorbar(width = 0.4) + geom_point(size = 4) +theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_x_discrete(name = 'Life History Strategy', labels = fg_labels2) +
  scale_y_continuous(expression(paste('Max. Growth Rate (kg yr'^-1, ' m'^-2,')'), limits=c(.6,1.2),
                                breaks = seq(0, 1.5, 0.2), labels = seq(0, 1.5, 0.2))) +
  theme_plant + theme(aspect.ratio = 0.75)
p
pdf(file.path(gdrive_path, "Figures/Fig_3/Slopes/Max_growth_by_fg.pdf"))
p
dev.off()

library(reshape2)
melt_pars <- melt(param_ci, id.vars=1:3)
#melt_pars$fg <- factor(melt_pars$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))
cast_pars <- dcast(melt_pars, fg+year~parameter+variable)

# Version from 27 Apr. Manually find x and y locations for the slope segment to be plotted.

segment_location <- pred_light_5groups %>%
  group_by(fg, year) %>%
  summarize(xmax = sum(light_area[which.max(diff(q50)):(1+which.max(diff(q50)))])/2,
            ymax = sum(q50[which.max(diff(q50)):(1+which.max(diff(q50)))])/2)

cast_pars <- left_join(cast_pars, segment_location)

dodge_width <- 0.03
error_bar_width <- 0.04
#---------------------------- Combined mean binned growth by light + fg --------------------------------
p_mean_1panel <- ggplot(obs_light_binned %>% filter(year == year_to_plot, mean_n_individuals > 10, !fg %in% c('alltree', 'unclassified'))) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50, color = fg)) +
  geom_errorbar(aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max, group = fg, color = fg, width = error_bar_width * width), position = position_dodge(width = dodge_width)) +
  geom_point(aes(x = bin_midpoint, y = mean, group = fg, fill = fg), size = 4, shape = 21, position = position_dodge(width = dodge_width)) +
  
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y) +
  scale_color_manual(name = 'Functional group', values = guild_fills_nb0, labels = fg_labels) +
  scale_fill_manual(values = guild_fills_nb, labels = fg_labels, guide = FALSE) +
  theme_plant #+

p_mean_1panel 
pdf(file.path(gdrive_path, "Figures/Fig_4/Growth by light/Combined mean growth.pdf"))
p_mean_1panel 
dev.off()

#---------------------------- Median binned growth by light + fg --------------------------------
p_median_panels <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels)) +theme_plant+
  #facet_wrap(~ fg, labeller = label_value,scales = 'free_y') +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
              aes(x = light_area, ymin = q025, ymax = q975,group=fg,color=NA,fill=fg), alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), 
            aes(x = light_area, y = q50,group=fg, color=fg), size=0.5) +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q25, yend = q75), size = 0.3) +
  #geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q025, yend = q975)) +
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), 
               aes(x = xmax * 0.5, xend = xmax * 2, y = ymax * 0.5, yend = ymax * 2), color = 'brown1', size = .5) +
  geom_point(shape=21, aes(x = bin_midpoint, y = median)) +
  scale_color_manual(values = guild_fills_nb ) +
  scale_fill_manual(values = fg_colors)+
  scale_x_log10(name = title_x, breaks=c(1,10,100)) + 
  scale_y_log10(name = title_y, labels=signif, breaks=c(0.01,0.1,1,10)) +
  theme_facet 

p_median_panels

pdf(file.path(gdrive_path, "Figures/Fig_4/Growth by light/Medians.pdf"))
p_median_panels
dev.off()

#---------------------------- Raw growth by light + fg --------------------------------
p_raw_panels <- ggplot(obs_light_raw %>% filter(year == year_to_plot, !fg %in% 'unclassified')) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels)) +theme_plant+
  geom_point(shape=21,aes(x = light_area, y = production_area, alpha = fg)) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975, group=fg,color=NA,fill=fg), alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50, group=fg, color=fg), size=0.25) +
  #scale_alpha_manual(values = c(0.15, 0.15, 0.05, 0.008, 0.008)) +
  scale_alpha_manual(values = c(0.15, 0.15, 0.15, 0.15, 0.15)) +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, breaks=c(0.001,0.01,1),labels=signif) +
  scale_color_manual(values = fg_colors) +
  scale_fill_manual(values = fg_colors)+
  theme_facet
#theme(strip.background = element_rect(fill=NA),
#    panel.border = element_rect(color = "black", fill=NA,  size=.75),legend.position = 'none',
#   panel.background = element_blank(),
#  strip.text.x = element_blank(),
# axis.text = element_text(size = 12, color = "black"), 
#axis.ticks.length=unit(0.2,"cm"),
#axis.title = element_text(size = 12)) 

p_raw_panels
pdf(file.path(gdrive_path, "Figures/Fig_4/Growth by light/Raw.pdf"))
p_raw_panels
dev.off()

#---------------------------- Extra growth by light + fg --------------------------------

p_mean_segments <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  #facet_wrap(~ fg, labeller = fg_labels) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels)) + 
  
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50), color = 'red') +
  geom_errorbar(aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max)) +
  geom_point(aes(x = bin_midpoint, y = mean), shape = 21, size = 2) +
  
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = xmax * 0.5, xend = xmax * 2, y = ymax * 0.5, yend = ymax * 2), color = 'green', size = 1) +
  scale_x_log10(name = title_x) + 
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = xmax * 0.5, xend = xmax * 2, y = ymax * 0.5^log_slope_q50 , yend = ymax * 2^log_slope_q50) , color = 'blue', size = 1) +
  scale_y_log10(name = title_y) +
  theme_plant +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA))
p_mean_segments

# 2. Plot with different panels for each functional group, and raw data

# I attempted to set an alpha scale so that the amount of transparency is roughly the same but the numbers may need to be tweaked
#pdf(file.path(gdrive_path, "Plots_J/New/Supplementals/p_raw_panels.pdf"))


#dev.off()
# 3. Plot with different panels for each functional group, and quantiles

# Problem!!
# Corrected 28 June 2019 by QDR because format of new "slopes" DF has changed.
ggplot(slopes %>% filter((dens_model == 3 & prod_model == 2) | (dens_model == 3 & is.na(prod_model)) | (is.na(dens_model) & prod_model == 2), !fg %in% 'Unclassified'), 
       aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg,scale = "free_y", labeller = label_value) +
  geom_ribbon() + geom_line()

unique(obs_light_binned$fg)



# 4. Plot with different panels for each functional group, and means

p_mean_panels <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels)) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975), alpha = 0.5, fill = 'red') +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50), color = 'red') +
  geom_errorbar(aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max)) +
  geom_point(aes(x = bin_midpoint, y = mean)) +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA))
p_mean_panels 



# 6. Plot with all functional groups on the same panel, and means



#7. Plot line segments of the maximum slope at correct location, and segments with slope=1 for isometry



# 8. Raw data plot converted to hexbin instead of plain scatter


#col1 <-colorRampPalette(c('forestgreen', 'goldenrod', 'indianred'))

#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=3)(50),
#                                  trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))

#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(c('khaki1', 'gold', 'red3'), bias=3)(50),
#                                  trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))

#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(c('aliceblue', 'forestgreen','khaki1','red3'), bias=3)(50),
#                        trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))
#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(c('forestgreen', 'red3'), bias=3)(50),
#                        trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))
#hex_scale_3b <- scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(9,'RdYlBu'), bias=0.1)(50)), guide = F)
hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000), 
                                             labels = c(1,10,100,1000), limits=c(1,5000))


p_hex_panels <- ggplot(obs_light_raw %>% filter(year == year_to_plot, !fg %in% 'unclassified', )) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels)) +
  geom_hex(aes(x = light_area, y = production_area)) +
  #geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
  #          aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.3) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), 
            aes(x = light_area, y = q50, color = fg)) +
  #geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q025, color = fg),  linetype = 'dotted') +
  #geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q975, color = fg),  linetype = 'dotted') +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, labels=signif) +
  hex_scale_log_colors + 
  scale_color_manual(name = 'Functional group',values = fg_colors, labels = fg_labels) +
  #scale_fill_manual(values = fg_colors, labels = fg_labels, guide = FALSE) +
  theme_plant +
  guides(color = FALSE) +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA),
        legend.position = c(0.8, 0.19))
p_hex_panels
#p <-set_panel_size(p_hex_panels, width=unit(10,"cm"), height=unit(6,"cm"))
#plot(p)
pdf(file.path(gdrive_path, "Figures/Fig_4/Growth by light/growth_light_heat_map.pdf"))
p_hex_panels
dev.off()

# 5. Plot with all functional groups on the same panel, and quantiles

dodge_width <- 0.00

p_median_1panel <- ggplot(obs_light_binned %>% 
                            filter(year == year_to_plot,mean_n_individuals >= 10, 
                                   !fg %in% c('alltree', 'unclassified'))) +
  geom_ribbon(data = pred_light_5groups %>% 
                filter(year == year_to_plot),
              aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.3) +
  geom_line(data = pred_light_5groups %>% 
              filter(year == year_to_plot), 
            aes(x = light_area, y = q50, color = fg)) +
  geom_line(data = pred_light_5groups[pred_light_5groups$fg == "fg5",]  %>% 
              filter(year == year_to_plot), aes(x = light_area, y = q50), color = "gray")+ 
  #geom_errorbar(aes(x = bin_midpoint, ymin = q25, ymax = q75, group = fg, color = fg), size = 0.75, width = 0, position = position_dodge(width = dodge_width)) +
  geom_errorbar(aes(x = bin_midpoint, ymin =ci_min, ymax = ci_max, group = fg,
                    color = fg),size = 0.5, width = 0.1, position = position_dodge(width = dodge_width)) +
  
  #geom_errorbar(aes(x = bin_midpoint, ymin = q025, ymax = q975, group = fg, color = fg), width = 0, position = position_dodge(width = dodge_width)) +
  geom_point(size=3.5, shape=21,color='black',aes(x = bin_midpoint,
                                                  y = median, group = fg, fill = fg),
             position = position_dodge(width = dodge_width)) +
  #geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = xmax * 0.5, xend = xmax * 2, y = ymax * 0.5, yend = ymax * 2), color = 'brown1', size = 1) +
  scale_x_log10(name = title_x, breaks=c(1,3,10,30,100,300))+ 
  scale_y_log10(name =  expression(atop('Growth per Crown Area',
                                        paste('(kg y'^-1, ' m'^-2,')')))) +
  scale_color_manual(name = 'Functional group', values = fg_colors, labels = fg_labels) +
  scale_fill_manual(values = fg_colors2, labels = fg_labels, guide = FALSE) +
  theme_plant +theme(aspect.ratio = 0.7) #+
# theme(panel.border = element_rect(fill=NA),
#      legend.position = c(0.2, 0.8))


p_median_1panel
p<-set_panel_size(p_median_1panel, width=unit(11.8,"cm"), height=unit(8.26,"cm"))
plot(p)



################################################################################################
# ------------------------------- Fig 5 Light Interception  ------------------------------------
################################################################################################


alltree_light_95 <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), as.character(NA)))

dbhbin1995 <- with(alltree_light_95, logbin(x = dbh_corr, n = 20))

lightperareafakebin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  do(quantile(.$light_received/.$crownarea, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightperareafakebin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  do(quantile(.$light_received/.$crownarea, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightperareafakebin_fg <- data.frame(fg = 'all', lightperareafakebin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightperareafakebin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

lightpervolfakebin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  do(quantile(.$light_received/.$crownvolume, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightpervolfakebin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  do(quantile(.$light_received/.$crownvolume, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightpervolfakebin_fg <- data.frame(fg = 'all', lightpervolfakebin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightpervolfakebin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

# Plot: raw data ----------------------------------------------------------

exl2 <- expression(paste('Light per Crown Area (W cm'^-1, 'm'^-2, ')', sep = ''))
exl <- expression(atop('Light per Crown Area', paste('(W cm'^-1, 'm'^-2, ')')))
exv <- expression(atop('Light per Crown Volume', paste('(W cm'^-1, ')')))
exv2 <- expression(paste('Light per Crown Volume (W cm'^-1, ')', sep = ''))
exd <- 'Diameter (cm)'

fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified')
full_names = c('Fast', 'Pioneer', 'Slow', 'Breeder', 'Medium', 'Unclassified')
labels <- setNames(full_names, fg_names)

#----------------------   Fig 5a: Light per crown area by diameter -----------------------------
# All together
p <- ggplot() + geom_point(alpha = 0.01, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownarea), color = 'chartreuse3') +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_plant +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p

p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p2)
pdf(file.path(gdrive_path,'Figures/Fig_5/Fig_5a.pdf'))
plot(p2)
dev.off()

#name = expression(atop('Total Light Intercepted',paste('(W m'^3, ' cm'^-1,' ha'^-1,')'))))  
# Supp: Each group
ggplot() +
  geom_point(alpha = 0.05, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownarea), color = 'chartreuse3') +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), 
                  aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg, ncol = 2, labeller = labeller(fg = labels)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  #theme_bw()
  theme_facet2



#----------------------   Fig 5b: Light per crown area by diameter -----------------------------

# All together
p <- ggplot() +
  geom_point(alpha = 0.01, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownvolume), color = 'chartreuse3') +
  geom_pointrange(data = lightpervolfakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv, breaks = c(10, 100)) +
  theme_plant

p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p2)
pdf(file.path(gdrive_path,'Figures/Fig_5/Fig_5b.pdf'))
plot(p2)
dev.off()



# ------------------- Each group
ggplot() +
  geom_point(alpha = 0.05, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownvolume), color = 'chartreuse3') +
  geom_pointrange(data = lightpervolfakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg, ncol = 2, labeller = labeller(fg = labels)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv) +
  theme_facet2
# Plot: hexagon plot ------------------------------------------------------

alpha_value <- 0.6
hexfill <- scale_fill_gradient(low = 'blue', high = 'red', trans = 'log', breaks = c(1,10,100,1000))
hex_scale_log_colors 
####### by area #######
# Each group
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownarea)) +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg, ncol = 2, labeller = labeller(fg = labels)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, breaks = c(1,10,100,1000), limits = c(1,1000)) +
  theme_facet2 +
  hex_scale_log_colors  +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# All together
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownarea)) +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, breaks = c(1,10,100,1000), limits = c(1,1000)) +
  theme_plant +
  hex_scale_log_colors  +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

####### by vol #######
# Each group
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownvolume)) +
  geom_pointrange(data = lightpervolfakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg, ncol = 2, labeller = labeller(fg = labels)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv, breaks = c(1,10,100,1000), limits = c(1,1000)) +
  theme_facet2 +
  hex_scale_log_colors  +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# All together
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownvolume)) +
  geom_pointrange(data = lightpervolfakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv, breaks = c(1,10,100,1000), limits = c(1,1000)) +
  theme_plant +
  hex_scale_log_colors  +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# Plot: total unscaled light energy by dbh --------------------------------


unscaledlightbydbhfakebin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  do(quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

unscaledlightbydbhfakebin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  do(quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

unscaledlightbydbhfakebin_fg <- data.frame(fg = 'all', unscaledlightbydbhfakebin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(unscaledlightbydbhfakebin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

alpha_value <- 0.6
hexfill2 <- scale_fill_gradient(low = 'forestgreen', high = 'navy', trans = 'log', breaks = c(1,3,10,30,100,300))
exl <- expression(atop('Intercepted Light', paste('per Individual (W)')))
exd <- 'Diameter (cm)'

labels = trans_format("log10", math_format(10^.x))



p <- ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received)) +
  geom_pointrange(data = unscaledlightbydbhfakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, breaks = c(1,100,10000, 1000000), limits = c(1,1000000), labels = trans_format("log10", math_format(10^.x))) +
  theme_plant +
  hex_scale_log_colors +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value))) +
  geom_abline(intercept = 0.903027, slope = 2.343608)
p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p2)
pdf(file.path(gdrive_path,'Figures/Fig_5/Indiv_light.pdf'))
plot(p2)
dev.off()


# Each group
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received)) +
  geom_pointrange(data = unscaledlightbydbhfakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg, ncol = 2, labeller = labeller(fg = fg_labels)) +theme_plant+
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, breaks = c(1,100,10000, 1000000), limits = c(1,1000000), labels = trans_format("log10", math_format(10^.x))) +
  theme_facet2 +
  hex_scale_log_colors +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# All together

# Very simple model
loglogregressions <- alltree_light_95 %>%
  group_by(fg) %>%
  do(model = lm(log10(light_received) ~ log10(dbh_corr), data = .))
library(broom)
loglogregressions$model
lapply(loglogregressions$model, summary)
#(Intercept)     0.903027   0.005262   171.6   <2e-16 ***
 #log10(dbh_corr) 2.343608   0.007530   311.2   <2e-16 ***


# Symmetry plot
# Compare total growth fitted slope (evaluated halfway between the high and low cutoffs) with total light fitted slope at the same point. 
# QDR / Forestlight / 20 June 2019



# Load parameter df to find the point at which to evaluate the fitted slopes, then load the fitted slopes

params <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/newpiecewise_paramci_by_fg.csv'), stringsAsFactors = FALSE)

xvalues <- params %>% 
  filter(variable == 'density', model == 3) %>%
  group_by(fg) %>%
  summarize(mid_cutoff = (mean[parameter == 'tau_high'] + mean[parameter == 'tau_low'])/2)

# Load fitted slopes of total growth and total light.

growth_slopes <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/newpiecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
light_slopes <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/totallightscaling/light_piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)

growth_slopes_atmiddle <- growth_slopes %>% 
  filter(variable == 'total_production', dens_model == 3, prod_model == 2) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

light_slopes_atmiddle <- light_slopes %>% 
  filter(variable == 'total_incoming_light', dens_model == 3, prod_model == 2) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

fg_full_names <- c('Fast', 'LL Pioneer', 'Slow', 'SL Breeder', 'Medium', 'All Trees', 'Unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

allslopes <- rbind(growth_slopes_atmiddle, light_slopes_atmiddle) %>%
  ungroup %>%
  mutate(fg = factor(fg, levels = fgs, labels = fg_full_names))
library(grid)
grob1 <- grobTree(textGrob("Light Capture", x = 0.75, y = 0.95, hjust = 0,
                           gp = gpar(col = "gold3", fontsize = 18))) 
grob2 <- grobTree(textGrob("Production", x = 0.75, y = 0.89, hjust = 0,
                           gp = gpar(col = "darkgreen", fontsize = 18)))# fontface="italic"
grob3 <- grobTree(textGrob("Energy Equivalence", x = 0.39, y = 0.47, hjust = 0,
                           gp = gpar(col = "black", fontsize = 18))) #, fontface = "bold")))
# Plot

p <- ggplot(allslopes %>% filter(!fg %in% 'Unclassified'), aes(x = fg, y = q50, ymin = q025, ymax = q975, color = variable)) +
  geom_hline(yintercept = 0, linetype = 'dotted', size = 1) +
  geom_errorbar(position = position_dodge(width = 0.5), size = 1, width = 0.7) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  labs(x = 'Life History Guild', y = 'Scaling Slopes') +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  scale_color_manual(values = c('gold2', 'darkgreen'), labels = c('Total Light', 'Total Growth')) +
  theme_plant + theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "italic", size = 15)) +
  annotation_custom(grob1) + annotation_custom(grob2) + annotation_custom(grob3)

p
pdf(file.path(gdrive_path, 'Figures/Fig_5/Growth Light symmetry.pdf'))
p
dev.off()
#-------------------------------------------------------------------------------------------


################################################################################################
# ------------------------------- Fig 6 Life History Ratios ------------------------------------
################################################################################################

# Loop through all the csv files and load them into R
fpdata <- file.path(gdrive_path, 'data/data_june2018_alternativecluster')
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'breeder_stats_bydiam_2census', 'breederscore_bin_bydiam_2census', 'breeder_stats_bylight_2census', 'breederscore_bin_bylight_2census', 'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census', 'fastslow_stats_bylight_2census', 'fastslowscore_bin_bylight_2census')

for (i in file_names) {
  assign(i, read.csv(file.path(fpdata, paste0(i,'.csv')), stringsAsFactors = FALSE))
}


error_bar_width <- 0.15
dodge_width <- 0

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.
global_diam_xlimits <- c(1, 316) # Maximum and minimum for x axis if it's diameter.

indivproductionbin_5census$fg <- factor(indivproductionbin_5census$fg, levels=c('fg1','fg2', 'fg3', 'fg4', 'fg5', 'unclassified', 'all'))
densitybin_5census$fg <- factor(densitybin_5census$fg, levels=c('fg1','fg2', 'fg3', 'fg4', 'fg5', 'unclassified', 'all'))

unique(densitybin_5census $fg)
fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('Fast','LL Pioneer', 'Slow', 'SL Breeder', 'Intermediate')
fg_labels2 <- c('Fast','LL Pioneer', 'Slow', 'SL Breeder', 'Intermediate','All')

p_dodge <- position_dodge(width = dodge_width)

geom_size = 3.5


###### Read in Data ######

# Load the 5 diameter census ones, can be used if needed
fastslow_stats_bydiam_5census <- read_csv(file.path(gdrive_path, 'data/data_june2018_alternativecluster/fastslow_stats_bydiam_5census.csv'))
breeder_stats_bydiam_5census <- read_csv(file.path(gdrive_path, 'data/data_june2018_alternativecluster/breeder_stats_bydiam_5census.csv'))

# Load the 2 yr diameter
fastslow_stats_bydiam_2census <- read_csv(file.path(gdrive_path, 'data/data_june2018_alternativecluster/fastslow_stats_bydiam_2census.csv'))
breeder_stats_bydiam_2census <- read_csv(file.path(gdrive_path, 'data/data_june2018_alternativecluster/breeder_stats_bydiam_2census.csv'))

# Load the 2 yr light
fastslow_stats_bylight_2census <-  read_csv(file.path(gdrive_path, 'data/data_june2018_alternativecluster/fastslow_stats_bylight_2census.csv'))
breeder_stats_bylight_2census <-  read_csv(file.path(gdrive_path, 'data/data_june2018_alternativecluster/breeder_stats_bylight_2census.csv'))


# ------------------------------- Fig 6 Relative Abundance & Pr=oduction ------------------------------------

# 2 samples by diameter

    # Rearrange Order so going up 
breeder_stats_bydiam_2censusb <- breeder_stats_bydiam_2census
breeder_stats_bydiam_2censusb[,2:7] <- 1/breeder_stats_bydiam_2censusb[,2:7]
breeder_stats_bydiam_2censusb <- as.matrix(breeder_stats_bydiam_2censusb )
    # Get rid of infinite values
breeder_stats_bydiam_2censusb[!is.finite(breeder_stats_bydiam_2censusb)] <- NA
breeder_stats_bydiam_2censusb <- as.data.frame(breeder_stats_bydiam_2censusb )
    # Add ID for combining
breeder_stats_bydiam_2census_id <- breeder_stats_bydiam_2censusb %>%
  mutate(ID = "Breeder-Pioneer")
fastslow_bydiam_2census_id <- fastslow_stats_bydiam_2census %>%
  mutate(ID = "Fast-Slow")
    # Combine
fastslow_stats_breeder_bydiam_2census <-rbind(breeder_stats_bydiam_2census_id,fastslow_bydiam_2census_id)
     # Plot by diameter
fill_ <- c("grey", "grey25")
fastslow_stats_breeder_bydiam_2census$ID <- as.factor(fastslow_stats_breeder_bydiam_2census$ID)
str(fastslow_stats_breeder_bydiam_2census)


# 5 samples by diameter
breeder_stats_bydiam_5censusb <- breeder_stats_bydiam_5census
breeder_stats_bydiam_5censusb[,2:7] <- 1/breeder_stats_bydiam_5censusb[,2:7]
breeder_stats_bydiam_5censusb <- as.matrix(breeder_stats_bydiam_5censusb )
breeder_stats_bydiam_5censusb[!is.finite(breeder_stats_bydiam_5censusb)] <- NA
breeder_stats_bydiam_5censusb <- as.data.frame(breeder_stats_bydiam_5censusb )
breeder_stats_bydiam_5census_id <- breeder_stats_bydiam_5censusb %>%
  mutate(ID = "Breeder")
fastslow_bydiam_5census_id <- fastslow_stats_bydiam_5census %>%
  mutate(ID = "Fast-Slow")
fastslow_stats_breeder_bydiam_5census <-rbind(breeder_stats_bydiam_5census_id,fastslow_bydiam_5census_id)
fill_ <- c("grey", "grey25")
fastslow_stats_breeder_bydiam_5census$ID <- as.factor(fastslow_stats_breeder_bydiam_5census$ID)
str(fastslow_stats_breeder_bydiam_5census)


# 2 samples by light 

fastslow_stats_bylight_2census_b <- fastslow_stats_bylight_2census
fastslow_stats_bylight_2census_b <- fastslow_stats_bylight_2census
fastslow_stats_bylight_2census_b[,2:7] <- 1/fastslow_stats_bylight_2census_b[,2:7] 
fastslow_stats_bylight_2census_b <- as.matrix(fastslow_stats_bylight_2census_b)
fastslow_stats_bylight_2census_b[!is.finite(fastslow_stats_bylight_2census_b)] <- NA
fastslow_stats_bylight_2census_b <- as.data.frame(fastslow_stats_bylight_2census_b)
fastslow_stats_bylight_2census_id <- fastslow_stats_bylight_2census %>%
  mutate(ID = "Fast-Slow")

breeder_stats_bylight_2census_b <- breeder_stats_bylight_2census
breeder_stats_bylight_2census_b <- breeder_stats_bylight_2census
breeder_stats_bylight_2census_b[,2:7] <- 1/breeder_stats_bylight_2census_b[,2:7] 
breeder_stats_bylight_2census_b <- as.matrix(breeder_stats_bylight_2census_b)
breeder_stats_bylight_2census_b[!is.finite(breeder_stats_bylight_2census_b)] <- NA
breeder_stats_bylight_2census_b <- as.data.frame(breeder_stats_bylight_2census_b)
breeder_stats_bylight_2census_id <- breeder_stats_bylight_2census_b %>%
  mutate(ID = "Breeder-Pioneer")

prod_ratio_light <- as.data.frame(rbind(fastslow_stats_bylight_2census_id,breeder_stats_bydiam_2census_id))

fastslow_stats_bydiam_2census_id <- fastslow_stats_bydiam_2census %>%
  mutate(ID = "Fast-Slow")

breeder_stats_bydiam_2census_b <- breeder_stats_bydiam_2census
breeder_stats_bydiam_2census_b[, 2:7] <- 1/breeder_stats_bydiam_2census_b[, 2:7]
breeder_stats_bydiam_2census_b <- as.matrix(breeder_stats_bydiam_2census_b)
breeder_stats_bydiam_2census_b[!is.finite(breeder_stats_bydiam_2census_b )] <- NA
breeder_stats_bydiam_2census_b <- as.data.frame(breeder_stats_bydiam_2census_b)
breeder_stats_bydiam_2census_id <- breeder_stats_bydiam_2census_b  %>%
  mutate(ID = "Breeder-Pioneer")

prod_ratio_diam <- as.data.frame(rbind(fastslow_stats_bydiam_2census_id, breeder_stats_bydiam_2census_id))



# ---------------------------------   Fig 6a  Production by light ------------------------------------
error_bar_width <- 0.15

p <- prod_ratio_light   %>%
  filter(mean_n_individuals > 10) %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, 
             ymin =production_ratio_min, ymax = production_ratio_max, fill = ID)) +
  geom_errorbar(width = error_bar_width) +
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  theme_plant+
  scale_x_log10(name = expression(paste('Light per Crown Area (W m'^-2,')')), limits=c(1,330), breaks=c(1, 10, 100)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.06,100),
                name = expression("Production Ratio"))+
  scale_y_log10(limits=c(0.006,100)) 
p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/Fig_6b_Production_light.pdf'))
plot(p1)
dev.off()

# -------------------- ----- Supp: Density ratio by light -----------------------------

p <- prod_ratio_diam   %>%
  filter(mean_n_individuals > 10) %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, 
             ymin =production_ratio_min, ymax = production_ratio_max, fill = ID)) +
  geom_errorbar(width = error_bar_width) +
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  theme_plant+
  scale_x_log10(name = expression(paste('Light per Crown Area (W m'^-2,')')), limits=c(1,330), breaks=c(1, 10, 100)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.06,100),
                name = expression("Production Ratio"))+
  scale_y_log10(limits=c(0.006,100)) + theme(axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             axis.ticks.y = element_blank())

p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/Sup_Production_diam.pdf'))
plot(p1)
dev.off()

# alt code
p <- bylight_2census    %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, 
             ymin = density_ratio_min, ymax = density_ratio_max, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  geom_errorbar(width = error_bar_width) + theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  
  scale_x_log10(limits=c(1,400),breaks=c(1,10,100), name = expression(paste('Light per Crown Area (W m'^-2,')'))) + 
  scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels=signif, limits=c(0.003,500),
                name = expression("Abundance Ratio"))
#scale_y_log10() + theme(axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#    axis.ticks.y = element_blank())
p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)

#----------------------------- Fig 6B Abundance by Diameter ----------------------------
# 5 yr samples

   # 5 sampling periods 5 yrs each
p <- fastslow_stats_breeder_bydiam_5census %>% 
  filter(density_ratio_mean > 0) %>%
  filter(mean_n_individuals > 10) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, 
             ymin = density_ratio_min, ymax = density_ratio_max, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  geom_errorbar(width = error_bar_width) + theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Breeder" = "black", "Fast-Slow" = "grey"))+
  
  scale_x_log10(limits=c(.7,330),breaks=c(1,10, 100), name = expression(paste('Diameter (cm)'))) + 
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.006,100),
                name = expression("Ratio")) +
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/Fig_6b_Density_5yr.pdf'))
plot(p1)
dev.off()





#----------------------------- Fig 6B alt: Abundance by Diameter for 2 sampling periods ----------------------------
#not plotting black points! (?)
p <- fastslow_stats_breeder_bydiam_2census %>% #fastslow_stats_breeder_bydiam_2census  %>%
  filter(density_ratio_mean > 0) %>%
  filter(mean_n_individuals > 10) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, 
             ymin = density_ratio_min, ymax = density_ratio_max, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  
  scale_x_log10(limits=c(.7,330),breaks=c(1,10, 100), name = expression(paste('Diameter (cm)'))) + 
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.006,100),
                name = expression("Ratio")) #+
  #theme(axis.title.y = element_blank(),axis.text.y = element_blank(),
  #      axis.ticks.y = element_blank())


p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/Fig_6b_Density_2yr.pdf'))
plot(p1)
dev.off()

# -------------------- ----- Supp: Production ratio by diameter  -----------------------------


p <- prod_ratio_diam  %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, 
             ymin =production_ratio_min, ymax = production_ratio_max, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_errorbar(width = error_bar_width) +
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  
  theme_plant+
  scale_x_log10(name = expression(paste('Diameter (cm)')), limits=c(1,150), breaks=c(1,3,10,30,100,300)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100), limits=c(0.01,200),
                name = expression("Production Ratio")) 
p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/Sup Ratios/Production_diam.pdf'))
plot(p1)
dev.off()


# -------------------- ------   Supp: PCA Rank   -----------------------------
## by Diameter

# Combine
fastslowscore_bin_bydiam_2census_id <- fastslowscore_bin_bydiam_2census %>%
  mutate(ID = "Fast-Slow")
breederscore_bin_bydiam_2census_id <- breederscore_bin_bydiam_2census %>%
  mutate(ID = "Breeder-Pioneer")
score_bin_bydiam <- as.data.frame(rbind(fastslowscore_bin_bydiam_2census_id, breederscore_bin_bydiam_2census_id))

p <- score_bin_bydiam %>%
  filter(mean_n_individuals > 10) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  #geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  
  scale_x_log10(name = 'Diameter (cm)', limits = c(1,300)) + 
  scale_y_continuous(limits=c(-2,2.5),breaks=c(-2,0,2),name = 'Slow-Fast PCA') 
p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/PCA_Scores/PCA_diam.pdf'))
plot(p1)
dev.off()


#-------------------------------  PCA rank by light ----------------------------------

# Combine
fastslowscore_bin_bylight_2census_id <- fastslowscore_bin_bylight_2census %>%
  mutate(ID = "Fast-Slow")
breederscore_bin_bylight_2census_id <- breederscore_bin_bylight_2census  %>%
  mutate(ID = "Breeder-Pioneer")
PCA_score_by_light <- as.data.frame(rbind(fastslowscore_bin_bylight_2census_id,breederscore_bin_bylight_2census_id))

p <- PCA_score_by_light %>%
  filter(mean_n_individuals > 10) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, fill = ID)) +
  #geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  scale_x_log10(limits=c(1,450),name = expression(paste('Light per Crown Area (W m'^-2,')'))) + 
  scale_y_continuous(limits=c(-2.2,1.5),name = 'Slow-Fast') 
p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/PCA_Scores/PCA_light.pdf'))
plot(p1)
dev.off()





#-----------------------------------------------------------------------------------
###################### Additional Light Plots  ########################



# Plots of "raw light scaling" aka incoming energy scaling
  # Fitted and predicted values, as well as comparing slopes to the production scaling
  # Use 1995 data
  


### 


# Fitted slopes
# light scaling
# QDR edit 27 June to add density slopes back in

rawlightslopes <- read.csv(file.path(gdrive_path,'data/data_piecewisefits/totallightscaling/light_piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
slopes <- read.csv(file.path(gdrive_path,'data/data_piecewisefits/newpiecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
rawlightslopes <- rbind(slopes %>% filter(variable=='density'), rawlightslopes)
#rawlightslopes <- read.csv(file.path(gdrive_path,'data/data_piecewisefits/rawlightpiecewise/rawlightpiecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
rawlightslopes$fg <- factor(rawlightslopes$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))
rawlightslopes$variable <- factor(rawlightslopes$variable, labels = c("Density", "Individual Incoming Energy", "Total Incoming Energy"))
colors <- c("sienna4", "yellowgreen", "springgreen4")

# Using 3 segment density and 2 segment production
ggplot(rawlightslopes %>% filter((dens_model == 3 & is.na(prod_model)) | (prod_model == 2 & is.na(dens_model)) | (dens_model == 3 & prod_model == 2), !fg %in% 'Unclassified'), 
       aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg,labeller = label_value) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = "springgreen4", size = 0.3) +
  geom_hline(yintercept = -2, linetype = 'dashed', col = "sienna4", size = 0.3) +
  geom_hline(yintercept = 2, linetype = 'dashed', col = "yellowgreen", size = 0.3) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_ribbon(alpha = 0.4, size = 0.3) +
  geom_line() +
  scale_x_log10(name = 'Diameter (cm)', expand = c(0,0)) +
  theme_bw() + theme(axis.text = element_text(color = "black", size = 12))+
  theme(strip.background = element_blank(), panel.grid = element_blank(), legend.position = 'bottom') +
  labs(y = 'Fitted Slope') +  #coord_fixed(ratio = .1)+
  ggtitle('Fitted slopes \n 3 segment density model and 2 segment incoming energy model')+
  theme(plot.title = element_text(hjust=0.5)) #+ theme(legend.title=element_blank())

slopes <- read.csv(file.path(gdrive_path,'data/data_piecewisefits/newpiecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
slopes$fg <- factor(slopes$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))
slopes$variable <- factor(slopes$variable, labels = c("Density", "Individual Growth", "Total Growth"))
colors <- c("sienna4", "yellowgreen", "springgreen4")

# Using 3 segment density and 2 segment production
ggplot(slopes %>% filter((dens_model == 3 & is.na(prod_model)) | (prod_model == 1 & is.na(dens_model)) | (dens_model == 3 & prod_model == 1), !fg %in% 'Unclassified'), 
       aes(x = dbh, y = q50, ymin = q025, ymax = q975, color = variable, fill = variable)) +
  facet_wrap(~ fg,labeller = label_value) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = "springgreen4", size = 0.3) +
  geom_hline(yintercept = -2, linetype = 'dashed', col = "sienna4", size = 0.3) +
  geom_hline(yintercept = 2, linetype = 'dashed', col = "yellowgreen", size = 0.3) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_ribbon(alpha = 0.4, size = 0.3) +
  geom_line() +
  scale_x_log10(name = 'Diameter (cm)', expand = c(0,0)) +
  theme_bw() + theme(axis.text = element_text(color = "black"))+
  theme(strip.background = element_blank(), panel.grid = element_blank(), legend.position = 'bottom') +
  labs(y = 'Fitted Slope') +  #coord_fixed(ratio = .1)+
  ggtitle('Fitted slopes \n 3 segment density model and 2 segment production model')+
  theme(plot.title = element_text(hjust=0.5)) #+ theme(legend.title=element_blank())

# Get the slope + ci for the middle portion
prod_slopes_10cm <- slopes %>% filter(variable == 'Total Growth', dbh == unique(slopes$dbh)[38], dens_model == 3, prod_model == 2)
light_slopes_10cm <- rawlightslopes %>% filter(variable == 'Total Incoming Energy', dbh == unique(rawlightslopes$dbh)[38], dens_model == 3, prod_model == 2)

ggplot(rbind(prod_slopes_10cm, light_slopes_10cm) %>% filter(!fg %in% c('All','Unclassified')), aes(x = fg, y = q50, ymin = q025, ymax = q975, group = variable, color = variable)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_errorbar(position = position_dodge(width = 0.2), width = 0.2) +
  geom_point(aes(color = c("red", "black"), position = position_dodge(width = 0.1), size = 6, shape = 21)) +
  theme_plant +
  ggtitle('Symmetry of total growth and total incoming energy patterns', 'Slope of scaling relationship at 10 cm dbh')

###
# Observed values with fitted lines superimposed (update to Figure 4)

# Load fitted values
rawlight_ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/rawlightpiecewise/rawlightpiecewise_ci_by_fg.csv', stringsAsFactors = FALSE) 

area_core <- 42.84

rawlight_ci_df$fg[rawlight_ci_df$fg == 'alltree'] <- 'all'

rawlight_pred_dens <- rawlight_ci_df %>%
  filter(variable == 'density') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

rawlight_fitted_indivprod <- rawlight_ci_df %>%
  filter(variable == 'production_fitted') %>%
  select(-variable)

rawlight_fitted_totalprod <- rawlight_ci_df %>%
  filter(variable == 'total_production_fitted') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

rawlight_pred_indivprod <- rawlight_ci_df %>%
  filter(variable == 'production') %>%
  select(-variable)

rawlight_pred_totalprod <- rawlight_ci_df %>%
  filter(variable == 'total_production') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

# Load observed values
rawlight_obs_totalprod <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/rawlightpiecewise/observed_light_bin.csv', stringsAsFactors = FALSE) %>% mutate(bin_value = bin_value/area_core)


# The below df is from the "total workflow" (Do not run)
# rawlight_observed_bin <- cbind(fg = 'all', lightreceivedbin_alltree_byyear[[2]]) %>%
#   rbind(map2_dfr(lightreceivedbin_fg_byyear, c('fg1','fg2','fg3','fg4','fg5','unclassified'), ~ cbind(fg = .y, .x[[2]])))
# write.csv(rawlight_observed_bin, file = '~/google_drive/ForestLight/data/data_piecewisefits/rawlightpiecewise/observed_light_bin.csv', row.names = FALSE)

# Modify fitted values so that they do not go above the data.
bin_maxes <- rawlight_obs_totalprod %>% group_by(fg) %>% summarize(max = max(bin_midpoint[bin_count > 10]))
rawlight_fitted_totalprod <- rawlight_fitted_totalprod %>%
  left_join(bin_maxes) %>%
  filter(dbh <= max)

prawlightfits <- ggplot(rawlight_obs_totalprod %>% filter(!is.na(fg), !fg %in% 'unclassified', bin_count > 10)) +
  geom_ribbon(data = rawlight_fitted_totalprod %>% filter(prod_model == 2, !fg %in% 'unclassified'), aes(x = dbh, ymin = q025, ymax = q975), fill = 'gray80') +
  geom_line(data = rawlight_fitted_totalprod %>% filter(prod_model == 2, !fg %in% 'unclassified'), aes(x = dbh, y = q50)) +
  geom_point(aes(x = bin_midpoint, y = bin_value)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = expression(paste('Total incoming light (W ha'^-1,')', sep=''))) +
  theme_bw() 

mycols <- c('black', RColorBrewer::brewer.pal(5,'Set1'))

prawlightfits_onefig <- ggplot(rawlight_obs_totalprod %>% filter(!is.na(fg), !fg %in% 'unclassified', bin_count > 10)) +
  geom_ribbon(data = rawlight_fitted_totalprod %>% filter(prod_model == 2, !fg %in% 'unclassified'), aes(x = dbh, ymin = q025, ymax = q975, fill = fg), alpha = 0.4) +
  geom_line(data = rawlight_fitted_totalprod %>% filter(prod_model == 2, !fg %in% 'unclassified'), aes(x = dbh, y = q50, color = fg)) +
  geom_point(aes(x = bin_midpoint, y = bin_value, color = fg)) +
  scale_x_log10(name = 'Diameter (cm)') +
  scale_y_log10(name = expression(paste('Total incoming light (W ha'^-1,')', sep='')), limits = c(1e2,1e5)) +
  theme_bw() +
  scale_color_manual(values=mycols) + scale_fill_manual(values=mycols)


fpfig <- '~/google_drive/ForestLight/figs/lightpowerlaws_feb2019'
ggsave(file.path(fpfig, 'totallightscaling_separateplots.png'), prawlightfits, height = 5, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'totallightscaling_oneplot.png'), prawlightfits_onefig, height = 5, width = 6, dpi = 300)

# light per volume

  # Light Vs Size Plot
  # Edited 21 March: include volume in addition to area.

#load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))
source(file.path(github_path, 'code/allfunctions27july.r'))
  
alltree_light_95 <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), as.character(NA)))

dbhbin1995 <- with(alltree_light_95, logbin(x = dbh_corr, n = 20))

lightperareafakebin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  do(quantile(.$light_received/.$crownarea, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightperareafakebin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  do(quantile(.$light_received/.$crownarea, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightperareafakebin_fg <- data.frame(fg = 'all', lightperareafakebin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightperareafakebin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

lightpervolfakebin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  do(quantile(.$light_received/.$crownvolume, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightpervolfakebin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  do(quantile(.$light_received/.$crownvolume, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightpervolfakebin_fg <- data.frame(fg = 'all', lightpervolfakebin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightpervolfakebin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

# Plot: raw data ----------------------------------------------------------

exl <- expression(paste('Light received per crown area (W m'^-2, ')', sep = ''))
exl <- expression(atop('Light per Crown Area',paste('(W m'^2, ' ha'^-1,')')))  

exv <- expression(atop('Light per Crown Volume',paste('(W m'^3, ' ha'^-1,')')))  
exd <- 'Diameter (cm)'

####### by area ######
# Each group
label_fg <- labeller(fg = c(all = "All", fg1 = 'Fast', fg2 = 'LL Pioneer',fg3 = "Slow", fg4 = "SL Breeder", fg5 = "Medium"))

p <- ggplot() +
  geom_point(alpha = 0.05, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownarea, color = fg)) +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg, labeller = label_fg)+ #labeller(fg = c(all = "all", fg1 = 'Fast', fg2 = 'LL Pioneer',fg3 = "Slow", fg4 = "SL Breeder", fg5 = "Medium"))) +#labeller = label_value) +
  #labeller(lightperareafakebin_fg$fg = c(fg1 = 'Fast', fg2 = 'LL Pioneer', fg3 = "Slow", fg4 = "SL Breeder", fg5 = "Medium")) +
  scale_color_manual(values = guild_fills_nb2)+
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, limits = c(1, 1000), breaks = c(1, 10, 100, 1000)) +
  theme(axis.title = element_text(size = 12)) +
  theme_bw() + 
  theme(axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 13)) + 
  theme(legend.position="none") +
  theme(strip.background = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = .75),
        strip.text.x = element_text(size = 12, face = "italic"),
        panel.grid = element_blank(), legend.position = 'none')
p
lapdf(file.path(gdrive_path, 'Figures/Growth_light/light_crown_area_fg.pdf'))
p
dev.off()
lightperareafakebin_fg$fg
# All together
p <- ggplot() +
  geom_point(alpha = 0.01, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownarea), color = 'chartreuse3') +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, limits = c(1, 1000), breaks = c(1, 10, 100, 1000)) +
  theme_plant 

p
pdf(file.path(gdrive_path, 'Figures/Growth_light/light_crown_area_all.pdf'))
p
dev.off()

####### by vol #######
# Each group
alltree_light_95$fg <- factor(alltree_light_95$fg , labels = c("Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium"))
lightpervolfakebin_fg$fg <- factor(lightpervolfakebin_fg$fg , labels = c( "all", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium"))

axis.title = element_text(size = 12)
p <- ggplot() +
  geom_point( alpha = 0.05, data = alltree_light_95 %>% filter(!is.na(fg)), 
             aes(x = dbh_corr, y = light_received/crownvolume,color = fg)) +
  geom_pointrange(data = lightpervolfakebin_fg %>% 
                    filter(!fg %in% 'all', !is.na(fg)), size = 0.3, 
                  aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg,labeller = label_value) +
  scale_x_log10(name = exd, breaks = c(1, 3, 10, 30, 100)) +
  scale_y_log10(name = exv) +
  theme(axis.title = element_text(size = 12)) +
  theme_bw() + 
  theme(axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 13)) + 
          theme(legend.position="none") +
  theme(strip.background = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = .75),
        strip.text.x = element_text(size = 12, face = "italic"),
        panel.grid = element_blank(), legend.position = 'none') +
scale_color_manual(values = guild_fills_nb2)
p
  
pdf(file.path(gdrive_path, 'Figures/Growth_light/light_crown_vol_fg.pdf'))
p
dev.off()

colors <- c("sienna4", "yellowgreen", "springgreen4")



# All together
p <- ggplot() +
  geom_point(alpha = 0.01, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownvolume), color = 'chartreuse3') +
  geom_pointrange(data = lightpervolfakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd, breaks = c(1, 3, 10, 30, 100)) +
  scale_y_log10(name = exv) +
  theme_plant 
p
pdf(file.path(gdrive_path, 'Figures/Growth_light/light_crown_vol_all.pdf'))
p
dev.off()



# Plot: hexagon plot ------------------------------------------------------

alpha_value <- 0.6
hexfill <- scale_fill_gradient(low = 'forestgreen', high = 'navy', trans = 'log', breaks = c(1,3,10,30,100,300))

####### by area #######
# Each group
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownarea)) +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_bw() +
  hexfill +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# All together
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownarea)) +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_bw() +
  hexfill +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

####### by vol #######
# Each group
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownvolume)) +
  geom_pointrange(data = lightpervolfakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv) +
  theme_bw() +
  hexfill +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# All together
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownvolume)) +
  geom_pointrange(data = lightpervolfakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv) +
  theme_bw() +
  hexfill +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# Plot: total unscaled light energy by dbh --------------------------------


unscaledlightbydbhfakebin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  do(quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

unscaledlightbydbhfakebin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  do(quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

unscaledlightbydbhfakebin_fg <- data.frame(fg = 'all', unscaledlightbydbhfakebin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(unscaledlightbydbhfakebin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))

alpha_value <- 0.6
hexfill2 <- scale_fill_gradient(low = 'forestgreen', high = 'navy', trans = 'log', breaks = c(1,3,10,30,100,300))
exl <- 'Light received (W)'
exd <- 'Diameter (cm)'

# Each group
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received)) +
  geom_pointrange(data = unscaledlightbydbhfakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_bw() +
  hexfill2 +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# All together
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received)) +
  geom_pointrange(data = unscaledlightbydbhfakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_bw() +
  hexfill2 +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# Very simple model
loglogregressions <- alltree_light_95 %>%
  group_by(fg) %>%
  do(model = lm(log10(light_received) ~ log10(dbh_corr), data = .))

lapply(loglogregressions$model, summary)


############### ############### Extra or Deprecated ##################################

####### Fig 4a Total Crown Volume ############
light_growth <- source(file.path(github_path,'code/plotting/lightdistributionplots.r'))
load(file.path(github_path, 'code/plotting/lightdistributionplots.r'))
###  Crown Area
error_bar_width <- 0.13

p <- crownvolumebins1995 %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_value), bin_value > 0) %>%
  filter(bin_count > 10) %>%
  mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
  mutate(bin_value = bin_value/area_core, bin_min = bin_min/area_core, bin_max = bin_max/area_core) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = bin_value, ymin = bin_min, ymax = bin_max, group = fg, fill=fg,color = fg)) +
  theme_plant+#geom_errorbar(aes(width = width), position=p_dodge) +
  #geom_abline(intercept=4, slope = -2/3, color ="darkgray",linetype="dashed", size=1.5)+
  geom_point(shape = 21, size = geom_size,  stroke = .5, color = "black")+
  scale_x_log10(name = 'Diameter (cm)', limits = c(1, 200), breaks=c(1,3,10,30,100,300)) +  
  scale_y_log10(limits = c(3, 3000),breaks=c(1, 10, 100, 1000, 10000), labels = signif,
                #scale_y_log10(labels = trans_format("log10", math_format(10^.x)), #limits = c(.2, 0.6),breaks=c(0.1, 0.2, 0.3, 0.4, 0.6), labels = signif,
                name = expression(atop('Total Crown Volume',paste('(m'^3, ' ha'^-1,')'))))  +  
  scale_color_manual(values = c('black', guild_fills_nb), labels = c('All', fg_labels), name = 'Functional group') +
  scale_fill_manual(values = c('black', guild_fills_nb), labels = c('All', fg_labels), name = 'Functional group') 

p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_4/total_crown_vol.pdf'))
plot(p1)
dev.off()

##### Fig 4b
# Load precalculated bin data.



###  Individual Growth vs Light [deprecated]

p <- indivprodperareabin_2census %>%
  filter(fg %in% fg_names, !is.na(mean), mean > 0) %>%
  filter(mean_n_individuals >= 10)%>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, group = fg, fill=fg,color = fg)) +
  theme_plant + geom_errorbar(aes(width = width), position = p_dodge) +
  geom_point(shape = 21, size = geom_size, color = "black", stroke = .5)+
  scale_x_log10(limits = c(1, 450),breaks=c(1,10,100,1000),  
                name = expression(paste('Light per Crown Area (W m'^-2,')')))+ 
  scale_y_log10(limits = c(.003, .15), breaks = c(0.003, 0.01, .03, 0.1, 0.3), labels = c(0.003, 0.01, .03, 0.1, 0.3),
                name = expression(atop('Growth per Crown',paste('Area (kg y'^-1, ' m'^-2,')')))) +
  scale_color_manual(values = guild_fills_nb0, labels = fg_labels, name = 'Functional Froup') +
  geom_abline(intercept = -3.85, slope = 1, color ="darkgray",linetype="dashed", size=1)+
  scale_fill_manual(values = guild_fills_nb, labels = fg_labels, name = 'Functional Froup') 
p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_4/Growth_per_Light.pdf'))
plot(p1)
dev.off()

# growth by light, normalized for volume

light_growth <- source(file.path(github_path,'code/plotting/lightdistributionplots.r'))

load(file.path(github_path, 'code/plotting/lightdistributionplots.r'))
###  Crown Area
error_bar_width <- 0.13
p <- crownareabin_2census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  filter(mean_n_individuals > 10) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, fill=fg,color = fg)) +
  theme_plant+geom_errorbar(aes(width = width), position=p_dodge) +
  #geom_abline(intercept=4, slope = -2/3, color ="darkgray",linetype="dashed", size=1.5)+
  geom_point(shape = 21, size = geom_size,  stroke = .5, color = "black")+
  scale_x_log10(name = 'Diameter (cm)', limits = c(1, 350), breaks=c(1,3,10,30,100,300)) +  
  scale_y_log10(limits = c(1, 10000),breaks=c(1, 10, 100, 1000, 10000), labels = signif,
                name = expression(atop('Total Crown Area',paste('(m'^2, ' ha'^-1,')'))))  +  
  scale_color_manual(values = c('black', guild_fills_nb), labels = c('All', fg_labels), name = 'Functional group') +
  scale_fill_manual(values = c('black', guild_fills_nb), labels = c('All', fg_labels), name = 'Functional group') 

p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_4/extra_crown_area_size.pdf'))
plot(p1)
dev.off()


crownareabin_2census_mid <- crownareabin_2census %>%
  filter(bin_midpoint > 4) %>%
  filter(bin_midpoint < 40) 
fitted_models = crownareabin_2census_mid %>% group_by(fg) %>% do(model = lm(log10(bin_yvalue) ~ log10(bin_midpoint), data = .))

fitted_models %>% tidy(model)
fitted_models %>% glance(model)
summary(lm1)


load(file.path(gdrive_path, 'data/area_and_volume_bins_1995.RData'))


library(broom)
crownvolumebins1995_mid <- crownvolumebins1995 %>%
  filter(bin_midpoint > 3) %>%
  filter(bin_midpoint < 30) %>%
  filter(fg == "all")

lm1 <- lm(log10(crownvolumebins1995_mid$bin_value)~ log10(crownvolumebins1995_mid$bin_midpoint))
summary(lm1)
plot(log10(crownvolumebins1995_mid$bin_value)~ log10(crownvolumebins1995_mid$bin_midpoint))

fitted_models = crownvolumebins1995_mid %>% group_by(fg) %>% 
  do(model = lm(log10(bin_value) ~ log10(bin_midpoint), data = .))

fitted_models %>% tidy(model)
fitted_models %>% glance(model)
summary(lm1)



