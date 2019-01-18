
# Some Plotting Code

theme_plant_0.6 <- theme(panel.grid = element_blank(), #for Density and Growth
                         aspect.ratio = .60,
                         axis.text = element_text(size = 19, color = "black"), 
                         axis.ticks.length=unit(0.2,"cm"),
                         axis.title = element_text(size = 19),
                         axis.title.y = element_text(margin = margin(r = 10)),
                         axis.title.x = element_text(margin = margin(t = 10)),
                         axis.title.x.top = element_text(margin = margin(b = 5)),
                         plot.title = element_text(size = 19, face = "plain", hjust = 10),
                         panel.border = element_rect(color = "black", fill=NA,  size=1),
                         panel.background = element_blank(),
                         #plot.margin = unit(c(1, 1, 1,1), "cm"),
                         legend.position = "none",
                         legend.key = element_rect(fill="transparent"),
                         
                         text = element_text(family = 'Helvetica')) 

theme_plant2 <- theme(panel.grid = element_blank(), #for Total Production
                      aspect.ratio = 1,
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
                      text = element_text(family = 'Helvetica')) 

 ### Plotting of main and supplemental figures

# Set paths to google drive forest light folder, and github forest light folder
# All paths below will be relative to these paths.
gdrive_path <- '~/google_drive/ForestLight'
github_path <- '~/Documents/GitHub/forestlight'

gdrive_path <- '/Users/jgradym/Google Drive/ForestLight'
github_path <- '/Users/jgradym/Documents/GitHub/forestlight'

library(tidyverse)
#---------------------------------------------------------------------------------------------
# Fig 1: hand drawn schematics 
#---------------------------------------------------------------------------------------------
# Fig 2: PCA Axes

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
guild_colors_nb <- c("#99AF3C", "#1D5128", "#1e3160", "#6BA3BF", "#D6D6D6")
fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('Fast','LL Pioneer', 'Slow', 'SL Breeder', 'Medium')

ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, fill = factor(fg5))) +
  geom_point(shape = 21, size = 3.5, color = "black") + theme_plant_0.6 +theme(aspect.ratio = 0.75)+
  labs(x = 'Slow to Fast', y = 'Breeders to Pioneers') +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group')+
  scale_fill_manual(values = guild_fills)


#---------------------------------------------------------------------------------------------
# Fig 3: Plotting functions for piecewise fits

# Define plotting functions -----


# Plot single model fit with multiple functional groups for density
plot_dens <- function(year_to_plot = 1995,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                      model_fit = 'pareto',
                      x_limits,
                      x_breaks = c(1, 3, 10, 30, 100,300),
                      y_limits,
                      y_breaks,
                      y_labels,
                      color_names = c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "ivory"),
                      x_name = 'Diameter (cm)',
                      y_name = expression(paste('Density (trees ha'^-1,'cm'^-1,')')),
                      obsdat = obs_dens,
                      preddat = pred_dens
) {
  
  require(dplyr)
  require(ggplot2)
  
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot, bin_count > 1) %>%
    filter(bin_value > 0)
  
  # Get minimum and maximum observed bin value for each group to be plotted
  # Delete points on the predicted line that are outside of this range
  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(dens_model %in% model_fit, prod_model == 2, fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs) 
  
  ggplot() +
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    geom_line(data = preddat[preddat$fg == "fg5",], aes(x = dbh, y = q50), color = "gray")+ # white circles get gray line
    geom_ribbon(data = preddat[preddat$fg == "fg5",], aes(x = dbh, ymin = q025, ymax = q975), fill = "gray", alpha = 0.4) +
    geom_point(data = obsdat, aes(x = bin_midpoint, y = bin_value, group = fg, fill=fg,size=2), shape=21,color="black") +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks,labels = y_labels) +
    scale_color_manual(values = color_names) +theme_plant+
    scale_fill_manual(values = color_names) 
  
  
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
                      color_names = c("#BFE046", "#267038", "#27408b", "#87Cefa", "ivory"),
                      x_name = 'Diameter (cm)',
                      y_name = expression(paste('Production (kg y'^-1,')')),
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
    filter(fg %in% fg_names, year == year_to_plot, !is.na(mean), mean_n_individuals > 1) %>%
    group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup
  
  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(prod_model %in% model_fit, dens_model == 1, fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs)
  
  ggplot() +
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    #geom_errorbar(data = obsdat, aes_string(x = 'bin_midpoint', ymin = error_quantiles[1], ymax = error_quantiles[2], group = 'fg', color = 'fg', width = 'width'), position = pos) +
    geom_line(data = preddat[preddat$fg == "fg5",], aes(x = dbh, y = q50), color = "gray")+ # white circles get gray line
    geom_ribbon(data = preddat[preddat$fg == "fg5",], aes(x = dbh, ymin = q025, ymax = q975), fill = "gray", alpha = 0.4) +
    geom_point(data = obsdat, aes_string(x = 'bin_midpoint', y = average, group = 'fg', fill = 'fg'),size=4,color="black",shape=21,position = pos) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels=y_labels) +
    scale_color_manual(values = color_names) +theme_plant+
    scale_fill_manual(values = color_names) 
  
  
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
                           color_names = c("black","#BFE046", "#267038", "#27408b", "#87Cefa", "ivory"),
                           x_name = 'Diameter (cm)',
                           y_name = expression(paste('Total production (kg ha'^-1,' y'^-1,')')),
                           obsdat = obs_totalprod,
                           preddat = fitted_totalprod
) {
  
  require(dplyr)
  require(ggplot2)
  obsdat <- obsdat %>%
    filter(fg %in% fg_names, year == year_to_plot) %>%
    filter(bin_value > 0)
  
  obs_limits <- obsdat %>%
    group_by(fg) %>%
    summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))
  
  preddat <- preddat %>%
    left_join(obs_limits) %>%
    filter(dens_model %in% model_fit_density, prod_model %in% model_fit_production, fg %in% fg_names, year == year_to_plot) %>%
    filter_at(vars(starts_with('q')), all_vars(. > min(y_limits))) %>%
    filter(dbh >= min_obs & dbh <= max_obs)
  
  ggplot() +
    geom_ribbon(data = preddat, aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
    geom_line(data = preddat, aes(x = dbh, y = q50, group = fg, color = fg)) +
    geom_line(data = preddat[preddat$fg == "fg5",], aes(x = dbh, y = q50), color = "gray")+ # white circles get gray line
    geom_ribbon(data = preddat[preddat$fg == "fg5",], aes(x = dbh, ymin = q025, ymax = q975), fill = "gray", alpha = 0.4) +
    geom_point(data = obsdat, aes(x = bin_midpoint, y = bin_value, size=2,group = fg, fill=fg), color = "black",shape=21) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) +
    scale_color_manual(values = color_names) +theme_plant2+
    scale_fill_manual(values = color_names) 
  
  
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

source(file.path(github_path, 'stan/piecewise_workflow/plottingfunctionspiecewise.r'))
# Create plots.
#Model fit 1 = pareto, 1 segment
#Model Fit 2  = 2 segments, etc

plot_dens(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = 3,
          x_limits = c(1, 260),
          y_limits = c(0.001, 3000),
          y_labels = c(0.001, 0.1, 10,1000),
          y_breaks = c(0.001, 0.1,  10, 1000))
ggsave(file.path(fpfig, 'fits_3partdensity.pdf'), height = 5, width = 6)

# Specify dodging with a certain width of error bar
# Model fit 1 = power law
# Model fit 2 = power law exp
plot_prod(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5'),
          model_fit = 2,
          x_limits = c(1, 280),
          y_limits = c(0.001, 2000),
          y_breaks = c(0.001,0.1, 10, 1000),
          y_labels = c(0.001,0.1,10,1000),
          error_bar_width = 0.01,
          dodge_width = 0.05)
ggsave(file.path(fpfig, 'fits_2partproduction.pdf'), height = 5, width = 6)

plot_totalprod(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
               model_fit_density = 3, 
               model_fit_production = 2,
               x_limits = c(0.9,250),
               y_limits = c(0.03, 200),
               y_breaks = c(0.1, 1, 10, 100),
               y_labels = c(0.1, 1, 10, 100),
               preddat = fitted_totalprod)
ggsave(file.path(fpfig, 'fits_3by2_totalproduction.pdf'), height = 5, width = 6)




# ------------------------------- Fig 4 Light Plots ------------------------------------

# Load precalculated bin data.
# Loop through all the csv files and load them into R
fpdata <- file.path(gdrive_path, 'data/data_june2018_alternativecluster')
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'breeder_stats_bydiam_2census', 'breederscore_bin_bydiam_2census', 'breeder_stats_bylight_2census', 'breederscore_bin_bylight_2census', 'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census', 'fastslow_stats_bylight_2census', 'fastslowscore_bin_bylight_2census')

for (i in file_names) {
  assign(i, read.csv(file.path(fpdata, paste0(i,'.csv')), stringsAsFactors = FALSE))
}

### Light Capture
# Set options for error bar widths and dodge amounts for all the plots
error_bar_width <- 0.13
dodge_width <- 0

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.
global_diam_xlimits <- c(1, 316) # Maximum and minimum for x axis if it's diameter.
guild_colors2 <- c("firebrick2","lightpink" , "royalblue4", "lightskyblue", "ivory" )
guild_colors3 <- c("firebrick2","lightpink" , "lightskyblue","royalblue4", "ivory","black" )

indivproductionbin_5census$fg <- factor(indivproductionbin_5census$fg, levels=c('fg1','fg2', 'fg3', 'fg4', 'fg5', 'unclassified', 'all'))
densitybin_5census$fg <- factor(densitybin_5census$fg, levels=c('fg1','fg2', 'fg3', 'fg4', 'fg5', 'unclassified', 'all'))

unique(densitybin_5census $fg)
guild_colors <- RColorBrewer::brewer.pal(5, 'Set1')
fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('Fast','Long-lived Pioneer', 'Slow', 'Short-lived Breeder', 'Intermediate')

fg_names <- paste('fg', 1:5, sep = '')
fg_labels2 <- c('Fast','Long-lived Pioneer', 'Slow', 'Short-lived Breeder', 'Intermediate','All')

p_dodge <- position_dodge(width = dodge_width)



lightreceivedbin_2census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = (bin_yvalue/area_core)/1000, bin_ymin = (bin_ymin/area_core)/1000, bin_ymax = (bin_ymax/area_core)/1000) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg,fill=fg)) +
  theme_plant2+ geom_errorbar(aes(width = width), position=p_dodge)+#geom_pointrange()  +
  geom_point(position=p_dodge,shape = 21, size = 4.5,  stroke = .5, color = "black")+
  scale_x_log10(name = 'Diameter (cm)', limits = c(1, 350), breaks=c(1,3,10,30,100,300)) +  
  scale_y_log10(position="right", name = expression(atop('Total Light Capture',paste('(kW ha'^-1,')'))),labels = signif,
                limits = c(0.04, 100), breaks=c(0.1, 1, 10, 100)) +
  scale_color_manual(values = c('black', guild_colors_nb), labels = c('All', fg_labels), name = 'Functional group') +
  scale_fill_manual(values = c('black', guild_fills_nb), labels = c('All', fg_labels), name = 'Functional group') 

###  Growth vs Light
indivprodperareabin_2census %>%
  filter(fg %in% fg_names, !is.na(mean), mean > 0) %>%
  filter(mean_n_individuals >= 10)%>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, group = fg, fill=fg,color = fg)) +
  theme_plant+geom_errorbar(aes(width = width), position=p_dodge) +
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black")+
  scale_x_log10(limits = c(1, 450),breaks=c(1,10,100,1000),  
                name = expression(paste('Light per Crown Area (W m'^-2,')')))+ 
  scale_y_log10(limits = c(.003, .5), breaks = c(0.003, 0.01, .03, 0.1, 0.3), labels = c(0.003, 0.01, .03, 0.1, 0.3),
               name = expression(atop('Production per Crown',paste('Area (kg y'^-1, ' m'^-2,')')))) +
  scale_color_manual(values = guild_colors_nb, labels = fg_labels, name = 'Functional Froup') +
  scale_fill_manual(values = guild_fills_nb, labels = fg_labels, name = 'Functional Froup') 
#dev.off()

###  Crown Area
crownareabin_2census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, fill=fg,color = fg)) +
  theme_plant+geom_errorbar(aes(width = width), position=p_dodge) +
  #geom_abline(intercept=4, slope = -2/3, color ="darkgray",linetype="dashed", size=1.5)+
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black")+
  scale_x_log10(name = 'Diameter (cm)', limits = c(1, 350), breaks=c(1,3,10,30,100,300)) +  
  scale_y_log10(limits = c(.1, 10000),breaks=c(0.1,1, 10, 100, 1000, 10000), labels = signif,
                name = expression(atop('Total Crown Area',paste('(m'^2, ' ha'^-1,')'))))  +  
  scale_color_manual(values = c('black', guild_colors_nb), labels = c('All', fg_labels), name = 'Functional group') +
  scale_fill_manual(values = c('black', guild_fills_nb), labels = c('All', fg_labels), name = 'Functional group') 

# ------------------------------- Fig 5 Relative Abundance  ------------------------------------


###### Life history Relative Abundance by Size 

# Fast vs Slow

# Filepath for fast/slow ratios

fastslowscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  #geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black", fill = "grey25")+
  scale_x_log10(name = 'Diameter (cm)', limits = c(1,300)) + 
  scale_y_continuous(limits=c(-2,2.5),breaks=c(-2,0,2),name = 'Slow-Fast') 
pdf(file.path(gdrive_path, "Plots_J/New/Ranks/Diameter/Fast Slow.pdf"))
dev.off()

### Breeder vs Pioneer

# Relative Abundance
breederscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean*-1, ymin = ci_min*-1, ymax = ci_max*-1)) +
  #geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black", fill = "grey25")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_log10(limits=c(1,300), name = 'Diameter (cm)') + 
  scale_y_continuous(limits=c(-4,2),name = 'Breeder-Performer') 
pdf(file.path(gdrive_path, "Plots_J/New/Ranks/Diameter/Breeders.pdf"))
dev.off()



###### Life history Relative Abundance by Light

### Fast vs Slow
fastslow_stats_bylight_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black", fill = "grey25")+
  
  scale_x_log10(limits=c(1,450),breaks=c(1,10,100), name = expression(paste('Light per Crown Area (W m'^-2,')'))) + 
  scale_y_log10(breaks = c(0.01,0.1,1,10), labels=signif, limits=c(0.003,2),
                name = expression(paste(frac("Fast","Slow"))))  
pdf(file.path(gdrive_path, "Plots_J/New/Abun_Light/Slow to Fast Abundance by Light.pdf"))
dev.off()


### Breeder vs Pioneer

breeder_stats_bylight_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = 1/density_ratio_mean, ymin = 1/density_ratio_min, ymax = 1/density_ratio_max)) +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black", fill = "grey25")+
  scale_x_log10(limits = c(1,450),breaks=c(1,10,100),
                name = expression(paste('Light per Crown Area (W m'^-2,')')))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  
  scale_y_log10(labels=signif, limits=c(0.1,400),breaks=c(0.01, 0.1, 1, 10, 100),
                name = expression(paste(frac("Short-Lived Breeder","Long-lived Pioneer")))) 
#pdf(file.path(gdrive_path, "Plots_J/New/Abun_Light/Breeder Abundance Ratio by light.pdf"))
dev.off()

# Ran





############################# SUPPLEMENTAL PLOTS #############################

#### Relative Abundance & Life history Ranks #####


###### Life history PCA Scores

## by Diameter

# Fast vs Slow

#pdf(file.path(gdrive_path, "Plots_J/New/Ranks/Diameter/Fast Slow.pdf"))
fastslowscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  #geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black", fill = "grey25")+
  scale_x_log10(name = 'Diameter (cm)', limits = c(1,300)) + 
  scale_y_continuous(limits=c(-2,2.5),breaks=c(-2,0,2),name = 'Slow-Fast') 
dev.off()

# Breeder vs Pioneer

breederscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean*-1, ymin = ci_min*-1, ymax = ci_max*-1)) +
  #geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black", fill = "grey25")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_log10(limits=c(1,300), name = 'Diameter (cm)') + 
  scale_y_continuous(limits=c(-4,2),name = 'Breeder-Performer') 


## by Light
# Fast vs Slow
fastslowscore_bin_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  #geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black", fill = "grey25")+
  scale_x_log10(limits=c(1,450),name = expression(paste('Light per Crown Area (W m'^-2,')'))) + 
  scale_y_continuous(limits=c(-2.2,0.5),name = 'Slow-Fast') 

# Breeder vs Pioneer
breederscore_bin_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean*-1, ymin = ci_min*-1, ymax = ci_max*-1)) +
  #geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black", fill = "grey25")+
  scale_x_log10(limits=c(1,450),name = expression(paste('Light per Crown Area (W m'^-2,')'))) +
  scale_y_continuous(breaks=c(-1,0,1),name = 'Breeder-Performer') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
#---------------------------------------------------------------------------------------------


# Life history Production Ratio by diameter

    # Fast vs Slow
#pdf(file.path(gdrive_path, "Plots_J/New/Prod_Diam/Fast vs Slow.pdf"))
fastslow_stats_bydiam_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin =production_ratio_min, ymax = production_ratio_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black", fill = "grey25")+
  theme_plant+
  scale_x_log10(name = expression(paste('Diameter (cm)')), limits=c(1,150), breaks=c(1,3,10,30,100,300)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100), limits=c(0.01,15),
                name = expression(paste(frac("Fast","Slow")))) 
#dev.off()

   # Breeder vs Pioneer 
breeder_stats_bydiam_2census
str(breeder_stats_bydiam_2census)
#pdf(file.path(gdrive_path, "Plots_J/New/Prod_Diam/Breeder Production by Diameter.pdf"))
breeder_stats_bydiam_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = 1/production_ratio_mean, ymin = 1/production_ratio_min, ymax = 1/production_ratio_max)) +
  geom_errorbar(width = error_bar_width) +
  theme_plant+ theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())+
  geom_point(shape = 21, size = 5.5,  stroke = .5, color = "black", fill = "grey25")+
  scale_x_log10(name = 'Diameter (cm)', limits=c(1,150), breaks=c(1, 3,10,30, 100,300)) + 
  scale_y_log10(labels = signif, limits=c(10^-2,10^3.35),breaks=c(0.01, 1, 100),
                name = expression(paste(frac("Long-lived Pioneer","Short-Lived Breeder")))) 
#dev.off()


## Life history Production Ratio by Light

      # Fast vs Slow
#pdf(file.path(gdrive_path, "Plots_J/New/Prod_Light/Fast to slow Production by light.pdf"))
fastslow_stats_bylight_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 5.5,  stroke = .5, color = "black", fill = "grey25")+
  scale_x_log10(limits=c(1,450),breaks=c(1,10,100),
                name = expression(paste('Light per Crown Area (W m'^-2,')')))+
  scale_y_log10(expression(paste(frac("Fast","Slow"))),
                limits=c(0.003,2),breaks=c(0.01,.1,1, 10), labels=signif) 
#dev.off()

     # Breeder vs Pioneer

#pdf(file.path(gdrive_path, "Plots_J/New/Prod_Light/Breeder Production Ratio by light2.pdf"))
breeder_stats_bylight_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = 1/production_ratio_mean, ymin = 1/production_ratio_min, ymax = 1/production_ratio_max)) +
  geom_errorbar(width = error_bar_width) +theme_plant+
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black", fill = "grey25")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_log10(limits = c(1,450),breaks=c(1,10,100),
                name = expression(paste('Light per Crown Area (W m'^-2,')')))+
  scale_y_log10(labels = signif, limits=c(10^-2.2,10^6),breaks=c(10^-2, 10^0, 10^2, 10^4, 10^6),
                name = expression(paste(frac("Long-lived Pioneer","Short-Lived Breeder")))) 
#dev.off()
#---------------------------------------------------------------------------------------------

######################## LOOIC of Piecewise Models  ########################

ics <- read.csv(file.path(fp, 'piecewise_ics_by_fg.csv'), stringsAsFactors = FALSE)
ics$fg <- factor(ics$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))

   # Density model

ggplot(ics %>% filter(prod_model == 1, criterion == 'LOOIC', variable == 'density', !fg %in% 'unclassified'), 
       aes(x = factor(dens_model), y = ic, ymin = ic - se_ic, ymax = ic + se_ic)) +
  facet_wrap(~ fg, labeller = label_value, scales = 'free_y') +
  geom_pointrange() +
  theme_bw() +
  theme(strip.background = element_blank(),panel.grid = element_blank()) +
  labs(x = 'Number of Segments in Density Function')
ggsave(file.path(fpfig, 'density_model_information_criteria.pdf'), height = 6, width = 9)

  # Production model

ggplot(ics %>% filter(dens_model == 1, criterion == 'LOOIC', variable == 'production', !fg %in% 'unclassified'), aes(x = factor(prod_model), y = ic, ymin = ic - se_ic, ymax = ic + se_ic)) +
  facet_wrap(~ fg, labeller = label_value,scales = 'free_y') +
  geom_pointrange() +
  theme_bw() +
  theme(strip.background = element_blank(), panel.grid = element_blank()) +
  labs(x = 'Number of Segments in Production Function')
ggsave(file.path(fpfig, 'production_model_information_criteria.pdf'), height = 6, width = 9)
#---------------------------------------------------------------------------------------------

######################## Plot Slopes vs Size for each life history group  ####################### 

slopes <- read.csv(file.path(fp, 'piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
slopes$fg <- factor(slopes$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))
slopes$variable <- factor(slopes$variable, labels = c("Density", "Individual Growth", "Total Growth"))
colors <- c("sienna4", "yellowgreen", "springgreen4")

# Using 3 segment density and 1 segment production
ggplot(slopes %>% filter(dens_model == 3, prod_model == 1, !fg %in% 'Unclassified'), 
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
ggtitle('Fitted slopes \n 3 segment density model and 2 segment production model')+
  theme(plot.title = element_text(hjust=0.5)) #+ theme(legend.title=element_blank())


# Using 3 segment density and 1 segment production
ggplot(slopes %>% filter(dens_model == 3, prod_model == 2, !fg %in% 'Unclassified'), 
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
  ggtitle('Fitted Slopes \n 3 Segment Density Model & 2 Segment Production Model')+
  theme(plot.title = element_text(hjust=0.5)) +
  theme(legend.spacing.x=unit(.2, "cm"))+  theme(legend.title=element_blank())+ theme(legend.text=element_text(size = 12))
  ggsave(file.path(fpfig, 'fitted_slopes_3partdensity_2partproduction.pdf'), height = 6, width = 9)
#---------------------------------------------------------------------------------------------


  
  theme_plant <- theme(panel.grid = element_blank(), 
                       aspect.ratio = .70,
                       axis.text = element_text(size = 17, color = "black"), 
                       axis.ticks.length=unit(0.15,"cm"),
                       axis.title = element_text(size = 17),
                       axis.title.y = element_text(margin = margin(r = 15)),
                       axis.title.x = element_text(margin = margin(t = 15)),
                       axis.title.x.top = element_text(margin = margin(b = 5)),
                       plot.title = element_text(size = 19, face = "plain", hjust = 10),
                       panel.border = element_rect(color = "black", fill=NA,  size=1),
                       panel.background = element_blank(),
                       plot.margin = unit(c(1, 1, 1,1), "cm"),
                       legend.position = "none",
                       legend.key = element_rect(fill="transparent"),
                       text = element_text(family = 'Helvetica')) 
  theme_plant1.1 <- theme(panel.grid = element_blank(), 
                          aspect.ratio = .70,
                          axis.text = element_text(size = 14, color = "black"), 
                          axis.ticks.length=unit(0.15,"cm"),
                          axis.title = element_text(size = 14),
                          axis.title.y = element_text(margin = margin(r = 0)),
                          axis.title.x = element_text(margin = margin(t = 10)),
                          axis.title.x.top = element_text(margin = margin(b = 5)),
                          plot.title = element_text(size = 14, face = "plain", hjust = 10),
                          panel.border = element_rect(color = "black", fill=NA,  size=.7),
                          panel.background = element_blank(),
                          plot.margin = unit(c(1, 1, 1,1), "cm"),
                          legend.position = "none",
                          legend.key = element_rect(fill="transparent"),
                          text = element_text(family = 'Helvetica')) 
  
  
  
  
######################## Hex Plot of Growth Scaling  ######################## 
# note 'object 'alltreedat' not found'
# Appropriate path needed 

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

plot_prod <- function(year_to_plot = 1990,
                      fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
                      full_names = c('Fast', 'Slow', 'Pioneer', 'Breeder', 'Medium', 'Unclassified'),
                      func_names = c('power law', 'power law\ntimes exponential'),
                      x_limits = c(1, 316),
                      x_breaks = c(1, 3, 10, 30, 100),
                      y_limits,
                      y_breaks,
                      x_name = 'Diameter (cm)',
                      y_name = expression(paste('Growth (kg yr'^-1,')')), line_types = c('dashed', 'solid'),
                      aspect_ratio = 1,
                      hex_scale = scale_fill_gradient(low = 'gray90', high = 'gray10', guide = FALSE),
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
    geom_line(data = preddat, aes(x = dbh, y = q50, group = prod_model, linetype = prod_model),size=0.25) +
    facet_wrap(~ fg, labeller = labeller(fg = labels)) +
    scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) +
    scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels=signif) +
    #scale_linetype_manual(values = line_types, name = 'Functional form') +
    scale_linetype_manual(values = line_types) +
    hex_scale +theme_plant1.1+
    coord_fixed(ratio = aspect_ratio) + guides(linetype = 'none')+
    theme(legend.position = c(0.85, 0.15), strip.background = element_blank())#, legend.text = element_blank())
  
  
}
# Original grayscale
#hex_scale_1 <- scale_fill_gradient(low = 'gray90', high = 'gray10', guide = FALSE)

# Biased grayscale to emphasize the variation in the hexagons with less data.
# Bias >1 does this.
#hex_scale_2 <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=10)(50))

# Biased color scale of red yellow and blue to do the same, using a brewer palette
# the bias is <1 this time because I had to reverse the scale.
#hex_scale_3 <- scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(9,'RdYlBu'), bias=0.3)(50)), guide = FALSE)

# Biased and customized color scale
#hex_scale_4 <- scale_fill_gradientn(colours = colorRampPalette(c('skyblue', 'goldenrod', 'indianred'), bias = 3)(50), guide = FALSE)

# Edit the hex scale argument to draw this with other color scales.
#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=1)(50), trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))
hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000), 
                                             labels = c(1,10,100,1000), limits=c(1,5000))


p <-plot_prod(year_to_plot = 1990,
              fg_names = c('fg1','fg2','fg3','fg4','fg5'),
              full_names = c('Fast', 'Slow', 'Pioneer', 'Breeder', 'Medium', 'Unclassified'),
              x_limits = c(1, 316),
              x_breaks = c(1, 10, 100),
              y_limits = c(5e-03, 1e04),
              y_breaks = c(.001, .1, 10, 1000),
              line_types = c('dashed', 'solid'),
              hex_scale = hex_scale_log_colors,
              aspect_ratio = 0.7) 
p 
#pdf(file.path(gdrive_path, "Plots_J/New/Supplementals/raw_growth_heat.pdf"))
#p
#dev.off()       

#-----------------------------------------------------------------------------------
###################### Additional Light Plots  ########################
library(tidyverse)
library(egg)
library(scales)


year_to_plot <- 1990 ### CHANGE THIS IF YOU WANT TO PLOT 1990

# Load data ----
fp <- '/Users/jgradym/Google Drive/ForestLight/data/data_forplotting_light_june2018'
#fp <- 'data/data_forplotting_light_june2018'
# New File path needed
obs_light_binned <- read.csv(file.path(fp, 'obs_light_binned.csv'), stringsAsFactors = FALSE)
obs_light_raw <- read.csv(file.path(fp, 'obs_light_raw.csv'), stringsAsFactors = FALSE)
pred_light <- read.csv(file.path(fp, 'pred_light.csv'), stringsAsFactors = FALSE)
param_ci <- read.csv(file.path(fp, 'lightbyarea_paramci_by_fg.csv'), stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)

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

# Create plots ---------

### -----------
### SET THESE OPTIONS

# Names of functional groups to display


# Axis titles
title_x <- expression(paste('Light per Crown Area (W m'^-2,')',sep=''))
title_y <- expression(paste('Growth (kg yr'^-1, ' m'^-2,')', sep=''))

# Colors
# these are shit colors but I just put them in as a placeholder
colors <- c('black',"#BFE046","#267038" , "#27408b", "#87Cefa", "ivory" ) # #BFE046= light green,#267038=dark green,27408b=dark blue,"#87Cefa" = light blue,    
year_to_plot = 1995
fg_colors <- c(fg1 = "#BFE046", fg2 = "#27408b" , fg3 = "#267038", fg4 =  "#87Cefa", fg5 = "gray70",  alltree = 'black') #unclassified = 'brown',
fg_colors2 <- c(fg1 = "#BFE046", fg2 = "#27408b" , fg3 = "#267038", fg4 =  "#87Cefa", fg5 = "ivory",  alltree = 'black') #unclassified = 'brown',

### -----

# 1. Plot of maximum slope by functional group

# Remove all tree and unclassified groups

p<-ggplot(param_ci %>% filter(year == year_to_plot, parameter %in% 'log_slope', !fg %in% c('alltree','unclassified')),
          aes(x = fg, y = q50, ymin = q025, ymax = q975)) +theme_plant+
  geom_hline(yintercept = 1, linetype = 'dotted', color = 'dodgerblue', size = 1) + 
  geom_errorbar(width = 0.1) + geom_point(size = 3) +
  scale_x_discrete(name = 'Life History Strategy', labels = fg_labels) +
  scale_y_continuous(name = 'Maximum Slope', limits=c(.25,1.1),breaks = seq(0, 1.5, 0.25), labels = seq(0, 1.5, 0.25)) +theme_plant
p
pdf(file.path(gdrive_path, "Plots_J/New/Supplementals/slopes.pdf"))
p
dev.off()
# Plot of intercept by FG
ggplot(param_ci %>% filter(year == year_to_plot, parameter %in% 'G', !fg %in% c('alltree','unclassified')),
       aes(x = fg, y = q50, ymin = q025, ymax = q975)) +theme_plant+
  geom_errorbar(width = 0.1) + geom_point() +
  scale_x_discrete(name = 'functional group', labels = fg_labels) +
  scale_y_continuous(name = 'growth vs light intercept', 
                    breaks = seq(0, 1.25, 0.25))+#, labels = seq(0, 1.25, 0.25)) +
  theme_plant
#panel_border(colour = 'black')

library(reshape2)
melt_pars <- melt(param_ci, id.vars=1:3)
cast_pars <- dcast(melt_pars, fg+year~parameter+variable)

# Version from 27 Apr. Manually find x and y locations for the slope segment to be plotted.

segment_location <- pred_light_5groups %>%
  group_by(fg, year) %>%
  summarize(xmax = sum(light_area[which.max(diff(q50)):(1+which.max(diff(q50)))])/2,
            ymax = sum(q50[which.max(diff(q50)):(1+which.max(diff(q50)))])/2)

cast_pars <- left_join(cast_pars, segment_location)

p_mean_segments <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels)) +
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
p_mean_segments
# 2. Plot with different panels for each functional group, and raw data

# I attempted to set an alpha scale so that the amount of transparency is roughly the same but the numbers may need to be tweaked
#pdf(file.path(gdrive_path, "Plots_J/New/Supplementals/p_raw_panels.pdf"))
p_raw_panels <- ggplot(obs_light_raw %>% filter(year == year_to_plot, !fg %in% 'unclassified')) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels)) +theme_plant+
  geom_point(shape=21,aes(x = light_area, y = production_area, alpha = fg)) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975, group=fg,color=NA,fill=fg), alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50, group=fg, color=fg), size=0.25) +
  scale_alpha_manual(values = c(0.15, 0.15, 0.05, 0.008, 0.008)) +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, breaks=c(0.01,0.1,1,10),labels=signif) +
  scale_color_manual(values = fg_colors) +
  scale_fill_manual(values = fg_colors)+
  theme(strip.background = element_rect(fill=NA),
        panel.border = element_rect(color = "black", fill=NA,  size=.75),legend.position = 'none',
        panel.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 12, color = "black"), 
        axis.ticks.length=unit(0.2,"cm"),
        axis.title = element_text(size = 12))
p_raw_panels
#dev.off()
# 3. Plot with different panels for each functional group, and quantiles
p_median_panels <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels)) +theme_plant+
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975,group=fg,color=NA,fill=fg), alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50,group=fg, color=fg), size=0.25) +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q25, yend = q75), size = 0.3) +
  #geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q025, yend = q975)) +
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), 
               aes(x = xmax * 0.5, xend = xmax * 2, y = ymax * 0.5, yend = ymax * 2), color = 'brown1', size = .5) +
  geom_point(shape=21, aes(x = bin_midpoint, y = median)) +
  scale_color_manual(values = fg_colors) +
  scale_fill_manual(values = fg_colors)+
  scale_x_log10(name = title_x, breaks=c(1,10,100)) + 
  scale_y_log10(name = title_y, labels=signif, breaks=c(0.01,0.1,1,10)) +
  theme(strip.background = element_rect(fill=NA),
        panel.border = element_rect(color = "black", fill=NA,  size=.75),legend.position = 'none',
        panel.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 12, color = "black"), 
        axis.ticks.length=unit(0.2,"cm"),
        axis.title = element_text(size = 12))
p_median_panels
pdf(file.path(gdrive_path, "Plots_J/New/Supplementals/light_growth.pdf"))
p_median_panels
dev.off()
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

dodge_width <- 0.03
error_bar_width <- 0.03

p_mean_1panel <- ggplot(obs_light_binned %>% filter(year == year_to_plot, mean_n_individuals > 10, !fg %in% c('alltree', 'unclassified'))) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50, color = fg)) +
  geom_errorbar(aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max, group = fg, color = fg, width = error_bar_width * width), position = position_dodge(width = dodge_width)) +
  geom_point(aes(x = bin_midpoint, y = mean, group = fg, fill = fg), size = 3, shape = 21, position = position_dodge(width = dodge_width)) +
  
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y) +
  scale_color_manual(name = 'Functional group', values = fg_colors, labels = fg_labels) +
  scale_fill_manual(values = fg_colors, labels = fg_labels, guide = FALSE) +
  theme_plant +
  theme(panel.border = element_rect(fill=NA),
        legend.position = c(0.2, 0.8))
p_mean_1panel 
# 7. Plot line segments of the maximum slope at correct location, and segments with slope=1 for isometry




# 8. Raw data plot converted to hexbin instead of plain scatter
#col1 <-colorRampPalette(c('skyblue', 'goldenrod', 'indianred'))

#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=3)(50),
    #                                  trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))

#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(c('khaki1', 'gold', 'red3'), bias=3)(50),
    #                                  trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))

#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(c('aliceblue', 'skyblue','khaki1','red3'), bias=3)(50),
              #                        trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))
#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(c('skyblue', 'red3'), bias=3)(50),
              #                        trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))
#hex_scale_3b <- scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(9,'RdYlBu'), bias=0.1)(50)), guide = F)
hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000), 
                                             labels = c(1,10,100,1000), limits=c(2,5000))


p_hex_panels <- ggplot(obs_light_raw %>% filter(year == year_to_plot, !fg %in% 'unclassified')) +
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
  theme_plant1.1+
  guides(color = FALSE) +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA),
        legend.position = c(0.85, 0.15))
p_hex_panels
#p<-set_panel_size(p_hex_panels, width=unit(10,"cm"), height=unit(6,"cm"))
#plot(p)
pdf(file.path(gdrive_path, "Plots_J/New/Supplementals/growth_light_heat.pdf"))
p_hex_panels
dev.off()

# 5. Plot with all functional groups on the same panel, and quantiles

dodge_width <- 0.00

p_median_1panel <- ggplot(obs_light_binned %>% filter(year == year_to_plot, mean_n_individuals >= 10, !fg %in% c('alltree', 'unclassified'))) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.3) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50, color = fg)) +
  #geom_errorbar(aes(x = bin_midpoint, ymin = q25, ymax = q75, group = fg, color = fg), size = 0.75, width = 0, position = position_dodge(width = dodge_width)) +
  geom_errorbar(aes(x = bin_midpoint, ymin =ci_min, ymax = ci_max, group = fg, color = fg),size = 0.5, width = 0.75, position = position_dodge(width = dodge_width)) +
  
  #geom_errorbar(aes(x = bin_midpoint, ymin = q025, ymax = q975, group = fg, color = fg), width = 0, position = position_dodge(width = dodge_width)) +
  geom_point(size=5, shape=21,color='black',aes(x = bin_midpoint, y = median, group = fg, fill = fg), position = position_dodge(width = dodge_width)) +
  #geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), aes(x = xmax * 0.5, xend = xmax * 2, y = ymax * 0.5, yend = ymax * 2), color = 'brown1', size = 1) +
  scale_x_log10(name = title_x, breaks=c(1,3,10,30,100,300))+ 
  scale_y_log10(name =  expression(atop('Growth per Crown Area',paste('(kg y'^-1, ' m'^-2,')')))) +
  scale_color_manual(name = 'Functional group', values = fg_colors, labels = fg_labels) +
  scale_fill_manual(values = fg_colors2, labels = fg_labels, guide = FALSE) +
  theme_plant #+
# theme(panel.border = element_rect(fill=NA),
#      legend.position = c(0.2, 0.8))


p_median_1panel
p<-set_panel_size(p_median_1panel, width=unit(11.8,"cm"), height=unit(8.26,"cm"))
plot(p)
pdf(file.path(gdrive_path, "Plots_J/New/Fig4/fig 4a.pdf"))
plot(p)
dev.off()




