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
gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/My Drive/ForestLight'))
github_path <- ifelse(Sys.info()['user'] == 'qread', '~/Documents/GitHub/MSU_repos', file.path('/Users/jgradym/Documents/GitHub/'))

gdrive_path2 <- file.path('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/My\\ Drive/ForestLight')

library(broom)
library(egg)
library(scales)
library(RColorBrewer)
library(gtable)
library(grid)
library(reshape2)
library(hexbin)
library(Hmisc, pos = 100)
library(rstan)
library(lme4)
library(sjstats)
library(rstanarm)
library(forestscaling) # Packaged all the functions and ggplot2 themes here!
library(tidyverse)
#devtools::install_github('qdread/forestscaling')

# Define color schemes and labels
guild_fills_all = c(fg1 = "#BFE046", fg2 =  "#267038", fg3 = "#27408b", fg4 = "#87Cefa", fg5 = "gray93", all = "black")
guild_fills_fg = c(fg1 = "#BFE046", fg2 =  "#267038", fg3 = "#27408b", fg4 = "#87Cefa", fg5 = "gray93")
guild_fills = c( "#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")

guild_colors_fg <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")
guild_colors_all <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray")

fg_labels <- c("fg1", "fg2","fg3", "fg4", "fg5")
fg_labels2 <- c("all", "fg1", "fg2","fg3", "fg4", "fg5")
fg_labels3 <- c("all", "fg1", "fg2","fg3", "fg4", "fg5", "unclassified")

guild_names <- paste('fg', 1:5, sep = '')
guild_labels <- c('Fast','Tall', 'Slow', 'Short', 'Medium')
guild_labels2 <- c('All', 'Fast','Tall', 'Slow', 'Short', 'Medium')
guild_labels3 <- c('All', 'Fast','Tall', 'Slow', 'Short', 'Medium',  "Unclassified")

guild_lookup <- data.frame(fg = c('fg1','fg2','fg3','fg4','fg5','all','unclassified'), 
                           fg_name = c('Fast','Tall','Slow','Short','Medium','All','Unclassified'))
guild_lookup
year_to_plot = 1995
geom_size <- 4

guide <- guides(fill = guide_legend(title = NULL), override.aes = list(fill = NA), color = F)
guide2 <- guides(color = guide_legend(title = NULL))


grob_a <- grobTree(textGrob("a", x = 0.04, y = 0.93,  hjust = 0,
                            gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob_b <- grobTree(textGrob("b", x = 0.04, y = 0.93,  hjust = 0,
                            gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 

#guide <- guides(color = guide_legend(title=NULL), fill = guide_legend(title=NULL)) #, color = F, override.aes=list(fill=NA))
grob_fast <- grobTree(textGrob("Fast", x = 0.04, y = 0.95,  hjust = 0,
                               gp = gpar(col = "#BFE046", fontsize = 15, fontface = "italic"))) 

grob_tall <- grobTree(textGrob("Tall", x = 0.04, y = 0.88,  hjust = 0,
                               gp = gpar(col = "#267038", fontsize = 15, fontface = "italic"))) 

grob_medium <- grobTree(textGrob("Medium", x = 0.04, y = 0.81,  hjust = 0,
                                 gp = gpar(col = "gray70", fontsize = 15, fontface = "italic"))) 

grob_slow <- grobTree(textGrob("Slow", x = 0.04, y = 0.74,  hjust = 0,
                               gp = gpar(col = "#27408b", fontsize = 15, fontface = "italic"))) 

grob_short <- grobTree(textGrob("Short", x = 0.04, y = 0.67,  hjust = 0,
                                gp = gpar(col = "#87Cefa", fontsize = 15, fontface = "italic"))) 

grob_all <- grobTree(textGrob("All", x = 0.04, y = 0.60,  hjust = 0,
                              gp = gpar(col = "black", fontsize = 15, fontface = "italic"))) 

fast_slow_fill <- c("#27408b", "#BFE046" )
fast_slow_fill2 <- c("#27408b", "#939E6C" ) #  #A4B662
fast_slow_fill3 <- c("#27408b", "#AABF5D" ) #  AABF5D
tall_slow_fill <- c("#27408b", "#267038")  #74B8CC; 9FAF65

scale_tall_slow <- scale_fill_gradientn(colours = tall_slow_fill,
                                        trans = 'log')#,# name = 
scale_fast_slow <- scale_fill_gradientn(colours = fast_slow_fill,
                                        trans = 'log')#,# name = 
scale_fast_slow3 <- scale_fill_gradientn(colours = fast_slow_fill3,
                                         trans = 'log')#,# name = 

# To add back the legend
theme_plant2 <- theme_plant() + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
theme_plant <- theme_plant() + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position = "right", legend.text = element_text(size = 14 ), 
        legend.key = element_blank())


#---------------------------- Data ---------------------------- 
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.R'))
fp <- file.path(gdrive_path, 'data/data_forplotting')

for (i in dir(fp, pattern = 'obs_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}


for (i in dir(fp, pattern = 'pred_|fitted_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}
#### Extract model output to get the fitted values, slopes, etc.
load(file.path(gdrive_path, 'data/data_piecewisefits/fits_bylight_forratio.RData'))
load(file.path(gdrive_path, 'data/data_binned/bin_object_singleyear.RData'))
load(file.path(gdrive_path, 'data/data_forplotting/light_scaling_plotting_data.RData'))
binned_data <- read_csv(file.path(gdrive_path, "data/data_binned/additional_bins_fg_year.csv"))

# source the extra extraction functions that aren't in the package
source(file.path(github_path, 'forestlight/stan/clean_workflow/model_output_extraction_functions.r'))
source(file.path(github_path, 'forestlight/stan/get_ratio_slopes_fromfit.R'))


##########################################################################################
############################ Life History PCA Plots ######################################
##########################################################################################

fgbci <- read.table(file.path(gdrive_path, 'data/Ruger/fgroups_dynamics_new.txt'), stringsAsFactors = FALSE)
# Correct functional groups so that: 1 fast, 2 tall, 3 slow, 4 short, 5 intermediate
# Old 1,2,3,4,5 --> New 2,3,1,4,5
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5))

# Currently X1new is high for slow species and low for fast species
# Currently X2new is high for pioneer species and low for breeder species
# Correct these
fgbci$PC_slow_to_fast <- -fgbci$X1new
fgbci$PC_breeder_to_pioneer <- fgbci$X2new


binomial <- paste(fgbci$genus, fgbci$species, sep = " ")
binomial
unique(binomial)
guild_lookup <- data.frame(fg = c('fg1','fg2','fg3','fg4','fg5','all','unclassified'), 
                           fg_name = c('Fast','Tall','Slow','Short','Medium','All','Unclassified'))


####--------- Fig 1A PCA Plot -----------------------
geom_size = 3
#geom_size = 4
Fig_3a <- ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, fill = factor(fg5))) +
  # geom_point(shape = 24, size = geom_size, color = "black") + 
  geom_point(shape = 24, size = geom_size, stroke = 0.3, color = "black") + 
  labs(x = 'Survivorship–Growth Tradeoff', y = 'Stature—Recruitment Tradeoff') +
  #theme_plant() + 
  theme_plant_small() + theme(aspect.ratio = 0.75) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_x_continuous(limits = c(-6,7), breaks = seq(-6,6,3))+
  scale_color_manual(values = guild_colors_fg, labels = fg_labels, name = 'functional group')+
  scale_fill_manual(values = guild_fills) 
Fig_3a 

p2  <- set_panel_size(Fig_3a , width=unit(14.3,"cm"), height=unit(14.3,"cm"))
p2 <- set_panel_size(Fig_3a , width=unit(10.25,"cm"), height=unit(7,"cm"))
p2 <- set_panel_size(Fig_3a , width=unit(10.25,"cm"), height=unit(8,"cm"))

grid.newpage()
grid.draw(p2)
#pdf(file.path(gdrive_path, 'Figures/New_main/life_history/life_history.pdf'))
#grid.draw(p2)
#dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/life_history/life_history.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/life_history/life_history.pdf')) 
)

#peering into things
guild_lookup
fgbci[fgbci$genus == "Luehea",]
fgbci[fgbci$genus == "Cecropia",]
fgbci[fgbci$genus == "Cavanillesia",]

#compare to Kitajima
fgbci[fgbci$genus == "Anacardium",] # Tall; in paper early-late; max height 30 m
fgbci[fgbci$genus == "Luehea",] # tall (almost medium); early-late; max height 30 m
fgbci[fgbci$genus == "Cecropia",]  #fast; high pioneer
fgbci[fgbci$genus == "Antirrhoea",] #NA; early successional;  max height <20m
fgbci[fgbci$genus == "Castilla",] #NA;  early successional; max height <20m

Antirrhoea

fgbci[fgbci$X2new > 2,] # tall
fgbci[fgbci$X2new < -2,] # short
fgbci[fgbci$X1new < -2,] # fast

fgbci$binomial <- paste(fgbci$genus, fgbci$species, sep = " ") 
slow <- fgbci[fgbci$fg5 == 3,] %>% dplyr::select(fg5, binomial, family, X1new) # slow
slow <- slow[order(desc(slow$X1new)),]
slow
#write_csv(slow, "~/Desktop/slow_tree_list.csv")


#anacardium excelsum  = slow
#Luehea seemannii = slow


########################################################################################
# ------------------------------- Fig x, stacked histogram ------------------------------
########################################################################################

#---------------------------
#---------- Percent abundance
#---------------------------
obs_dens$fg2 <- factor(obs_dens$fg, levels = c("fg3", "fg4","fg2", "fg1", "fg5", "all", "unclassified")) # Loo
obs_dens$fg2 <- factor(obs_dens$fg, levels = c("fg3", "fg4", "fg1", "fg2", "fg5", "all", "unclassified")) # Loo
obs_dens$fg2 <- factor(obs_dens$fg, levels = c("fg5", "fg3", "fg4", "fg1", "fg2",  "all", "unclassified")) # Loo

obs_dens$fg2<- factor(obs_dens$fg, levels = c("fg3", "fg1","fg2",  "fg4", "fg5", "all", "unclassified")) # Loo
obs_dens$bin_count2<- obs_dens$bin_count
obs_dens$bin_count2[is.na(obs_dens$bin_count2)] = 0
unique(obs_dens$fg)
str(obs_richnessbydiameter)
str(obs_dens)

#area plot

obs_richnessbydiameter$fg2 <- factor(obs_richnessbydiameter$fg, levels = c("fg3", "fg4","fg2","fg1",  "fg5", "all", "unclassified")) # Loo
obs_richnessbydiameter = obs_richnessbydiameter %>%
  group_by(fg) %>%
  mutate(perc_abun = n_individuals/obs_richnessbydiameter$n_individuals[obs_richnessbydiameter$fg == "all"])
(p = ggplot(data = obs_richnessbydiameter %>% filter(!fg2 %in% c("all", "unclassified")),
             aes(x = bin_midpoint, y = abundance_by_bin_width, fill = fg2)) +
    geom_area(position = position_fill(reverse = TRUE)) +
    scale_fill_manual(values = guild_fills_fg) +
    theme_plant()+ theme(legend.position = "none") +
    #theme(aspect.ratio = 0.5) +
    #theme_no_x() +
    scale_x_log10(name = "Stem Diameter (cm)", breaks = c(1, 3, 10, 30, 100, 300), 
                  #limits = c(0.8, 400),
                  expand = c(0,0)) +
    scale_y_continuous(labels = scales::percent,
                       position = "left", 
                       name = "Percent Abundance", expand = c(0,0)) 
  
)
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path,'Figures/New_main/percent/percent_abundance.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/percent/percent_abundance'), 
                    file.path(gdrive_path2,'Figures/New_main/percent/percent_abundance')) 
)

#bar plot
ggplot(data = obs_dens %>%  filter(!fg2 %in% c("all", "unclassified")),
       aes(x = bin_midpoint, y = bin_count, fill = fg2)) + # could use bin_value instead of bin_count - no difference
  geom_col(position = position_fill(reverse = TRUE), width = 0.125) +
  scale_fill_manual(values = guild_fills_fg) +
  theme_plant +
  scale_x_log10(name = "Stem Diameter (cm)", breaks = c(1, 3, 10, 30, 100, 300)) +
  scale_y_continuous(labels = scales::percent, name = "Relative Abundance") 


#-----------------------------
#---------- Percent Richness
#---------------------------
obs_richnessbydiameter$fg2 <- factor(obs_richnessbydiameter$fg, levels = c("fg3", "fg4","fg2","fg1",  "fg5", "all", "unclassified")) # Loo
#obs_richnessbydiameter$fg2 <- factor(obs_richnessbydiameter$fg, levels = c("fg3", "fg2","fg1", "fg4", "fg5", "all", "unclassified")) # Loo


#area plot
ggplot(data = obs_richnessbydiameter%>% filter(!fg2 %in% c("all", "unclassified")),
       aes(x = bin_midpoint, y = richness, fill = fg2)) +
  geom_area(position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = guild_fills_fg) +
  theme_plant +
  scale_x_log10(name = "Stem Diameter (cm)", breaks = c(1, 3, 10, 30, 100, 300), expand = c(0,                                                                                           0)) +
  scale_y_continuous(labels = scales::percent, name = "Relative Richness", expand = c(0,0)) 


#---------- Bar plot
ggplot(data = obs_richnessbydiameter %>% filter(!fg2 %in% c("all", "unclassified")),
       aes(x = bin_midpoint, y = richness, fill = fg2)) +
  geom_col(position = position_fill(reverse = TRUE)) +
  #geom_col(position = "fill") +
  scale_fill_manual(values = guild_fills_fg) +
  theme_plant +
  scale_x_log10(name = "Stem Diameter (cm)", breaks = c(1, 3, 10, 30, 100, 300), expand = c(0,0)) +
  scale_y_continuous(labels = scales::percent, name = "Relative Richness", expand = c(0,0)) 

ggplot(data = obs_richnessbydiameter %>% filter(!fg2 %in% c("all", "unclassified")),
       aes(x = bin_midpoint, y = richness, fill = fg2)) +
  geom_col(position = position_fill(reverse = TRUE), width = .125) +
  scale_fill_manual(values = guild_fills_fg) +
  theme_plant +
  scale_x_log10(name = "Stem Diameter (cm)", breaks = c(1, 3, 10, 30, 100, 300), expand = c(0,0)) +
  scale_y_continuous(labels = scales::percent, name = "Relative Richness", expand = c(0,0)) 


###################################################################################
################# Biomass growth heat map #####################
###################################################################################
raw_prod <- do.call(rbind, map2(alltreedat, seq(1985,2010,5), function(x, y) cbind(year = y, x %>% filter(!recruit) %>% select(fg, dbh_corr, production))))
raw_prod <- raw_prod %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg', fg), 'unclassified'))

hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000,10000), 
                                             labels = c(1,10,100,1000,10000), limits=c(1,30000))

#--error bars
growth_range = obs_indivprod %>%
  filter(mean_n_individuals >= 20,fg == "all") %>%
  group_by(bin_midpoint) %>%
  summarize(
    min = min(mean, na.rm = T),
    max = max(mean, na.rm = T))

(p <- ggplot() +
    geom_hex(data = raw_prod %>% filter(year == 1995), aes(x = dbh_corr, y = production)) + #+
    geom_point(data = obs_indivprod %>% 
                 filter(year == 1995, mean_n_individuals >= 20,fg == "all"),
               aes(x = bin_midpoint, y = mean), shape = 21, size = 2.5,  color = "black", fill = "black") +
    geom_errorbar(data = growth_range, aes(x = bin_midpoint, ymin = min, ymax = max ), color = "black", width = 0) +
    geom_smooth(data = raw_prod %>% filter(year == 1995, dbh_corr < 156), aes(x = dbh_corr, y = production),
              color = "black", method = "lm", se = T, size = .5, alpha = 0.2) +
    hex_scale_log_colors +
    scale_x_log10(name = "Stem Diameter (cm)", breaks = c(1, 3, 10, 30, 100, 300)) + 
    scale_y_log10(name = expression(paste("Growth (kg yr"^-1, ")")), position = "left",
                  #breaks = c(0.0001, 0.01, 1, 100, 10000), labels = c("0.0001", 0.01, 1, 100, "10,000"),
                  breaks = c(0.001, 0.01, 0.1,1, 10, 100, 1000), labels = c("0.001", 0.01, 0.1, 1, 10,100, "1000"),
                  limits = c(0.0005, 3000)) +
    theme_plant_small() + theme(legend.position = "right", legend.text = element_text(size = 14), legend.title = element_text(size = 14))
)
p_hex <- set_panel_size(p, width=unit(14.3,"cm"), height=unit(14.3,"cm"))
p_hex  <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p_hex)

raw_prod$fg
sum(!is.na(raw_prod$production[raw_prod$year == 1995]))

#pdf("/Volumes/GoogleDrive/ForestLight/Figures/New_main/growth_hex2.pdf")

pdf(file.path(gdrive_path,'Figures/new_main/Final_figs/Fig_3/growth_hex.pdf'))
grid.draw(p_hex)
dev.off()
system2(command = "pdfcrop", 
        args  = c(file.path(gdrive_path2,'Figures/new_main/Final_figs/Fig_3/growth_hex.pdf'), 
                  file.path(gdrive_path2,'Figures/new_main/Final_figs/Fig_3/growth_hex.pdf')) 
)




###################################################################################
################# Abundance, Richness, Production   #####################
###################################################################################

#----------------------------------------------------------------------------------
#--------------------------   Abundance          ----------------------------------
#----------------------------------------------------------------------------------

obs_dens <- obs_dens %>%
  filter(bin_count >= 20) 

#--error bars
fg_mid = c("fg", "bin_midpoint")
dens_range <- obs_dens %>%
  filter(fg != "unclassified",bin_count >= 20) %>%
  group_by(across(fg_mid )) %>%
  summarize(
    min = min(bin_value),
    max = max(bin_value)
  )
dens_range

plot_dens <- function(year_to_plot = 1995, 
                        fg_names = c("fg1", "fg2", "fg3", "fg4", "fg5", "all"),
                        model_fit = 1, 
                        x_limits, 
                        x_breaks = c(1, 3, 10, 30, 100, 300),
                        y_limits, 
                        y_breaks, 
                        y_labels, 
                        fill_names = guild_fills_all, #c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                        color_names = guild_colors_all, #c("black","#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),  
                        x_name = "Stem Diameter (cm)",
                        y_name = expression(paste("Density (n ha"^-1, "cm"^-1, ")")), 
                        geom_size = 4, 
                        obsdat = obs_dens, 
                        preddat = pred_dens, 
                        plot_abline = TRUE, 
                        abline_slope = -2, 
                        abline_intercept = -10,
                        dodge_width = 0.07) 
{
  pos <-  ggplot2::position_dodge(width = dodge_width)
  obsdat <- obsdat %>% dplyr::filter(fg %in% fg_names, year == year_to_plot, bin_count >= 20) %>% 
    dplyr::filter(bin_value > 0) %>% 
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1')))
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% 
    dplyr::summarize(min_obs = min(bin_midpoint), 
                     max_obs = max(bin_midpoint))
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% 
    dplyr::filter(dens_model %in% model_fit, fg %in% fg_names, 
                  year == year_to_plot) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% 
    dplyr::filter(dbh >=  min_obs & dbh <= max_obs) %>%  
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1')))
  p <- ggplot2::ggplot() + 
    ggplot2::geom_ribbon(data = preddat, 
                         ggplot2::aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), 
                         alpha = 0.4)
  if (plot_abline) {
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", linetype = "dashed", 
                                  size = 0.75)
  }
  p + ggplot2::geom_line(data = preddat, 
                         ggplot2::aes(x = dbh, y = q50, group = fg, color = fg)) + 
    ggplot2::geom_line(data = preddat[preddat$fg == "fg5", ], 
                       ggplot2::aes(x = dbh, y = q50), color = "gray") + 
    ggplot2::geom_point(data = obsdat, position = pos,
                        ggplot2::aes(x = bin_midpoint, y = bin_value, group = fg, fill = fg), 
                        size = geom_size, shape = 21, color = "black") + 
    ggplot2::geom_errorbar(data = dens_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), width = 0) +
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) + 
    ggplot2::scale_color_manual(values = color_names) + 
    ggplot2::scale_fill_manual(values = fill_names) + 
    theme_plant() + theme_no_x()
}

(p <- plot_dens(year_to_plot = 1995,
                 fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                 model_fit = DENS,
                 dodge_width = 0.0,
                 x_limits = c(.9, 160),
                 y_limits = c(0.007, 20000),
                 x_breaks = c(1, 10, 100),
                 y_labels = c(0.001, 0.1, 10,1000),
                 y_breaks = c(0.001, 0.1,  10, 1000)))
#p <- p +annotation_custom(grob_text_b)

p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/new_main/abundance/abundance.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/new_main/abundance/abundance.pdf'), 
                    file.path(gdrive_path2,'Figures/new_main/abundance/abundance.pdf')) 
)
#----------------------------------------------------------------------------------
#--------------------------   Richness  -------------------------------------------
#----------------------------------------------------------------------------------

# filter fitted sizes by binned sizes
max_dbh_fg <- obs_richnessbydiameter  %>%
  filter(n_individuals >= 20) %>%
  group_by(fg) %>%
  summarize(max_dbh = max(bin_midpoint) +2 )

min_dbh_fg <-obs_richnessbydiameter  %>%
  filter(n_individuals >= 20) %>%
  group_by(fg) %>%
  summarize(min_dbh = min(bin_midpoint) -0.1)

fitted_richnessbydiameter_filtered <- fitted_richnessbydiameter %>%
  left_join(max_dbh_fg) %>%
  left_join(min_dbh_fg) %>%
  group_by(fg) %>%
  filter(dbh <= max_dbh, dbh >= min_dbh)

unique(fitted_richnessbydiameter_filtered$year)

fg_mid = c("fg", "bin_midpoint")
rich_range <- binned_data %>%
  filter(fg != "unclassified", abundance >= 20) %>%
  group_by(across(all_of(fg_mid_groups ))) %>%
  summarize(
    min = min(richness_by_bin_width),
    max = max(richness_by_bin_width)
  )
ggplot2::geom_errorbar(data = dens_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), width = 0) +
  
#---- Plot Richness ------

(p_rich_cm <- ggplot() + 
   theme_plant() +
  scale_x_log10(name = 'Stem Diameter (cm)',
                 limit = c(.9, 160)) + 
   scale_y_log10(labels = signif,
                limit = c( .1, 1500), 
                 position = "left",
                 name = expression(paste("Richness (50 ha"^-1," cm"^-1, " )"))) +
   scale_fill_manual(values = guild_fills_all) +
   scale_color_manual(values = guild_colors_all) +
   geom_ribbon(data = fitted_richnessbydiameter_filtered  %>% 
                 arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))),
               aes(x = dbh, ymin = q025, ymax = q975,
                   group = fg, fill = fg), alpha = 0.4) +
   geom_line(data = fitted_richnessbydiameter_filtered  %>% 
               arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))),
             aes(x = dbh, y = q50, group = fg, color = fg)) +
   geom_point(data = binned_data %>% 
                 arrange(desc(fg)) %>%
                 filter(!fg %in% 'unclassified', abundance >= 20) %>% #n_individuals >= 20
                 filter(year == "1995") %>%
                 arrange(desc(fg)), 
               aes(x = bin_midpoint, y = richness_by_bin_width, fill = fg, color = fg),
               shape = 21, size = 4, color = "black") + theme_no_x() +
  geom_errorbar(data = rich_range %>%  arrange(desc(fg)), aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), width = 0)) 
  
 
p1 <- set_panel_size(p_rich_cm, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/new_main/Final_figs/Fig_4/richness/Richness.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args  = c(file.path(gdrive_path2,'Figures/new_main/Final_figs/Fig_4/richness/Richness.pdf'), 
                  file.path(gdrive_path2,'Figures/new_main/Final_figs/Fig_4/richness/Richness.pdf')) 
)

#----------------------------------------------------------------------------------
#--------------------------   Productivity  ---------------------------------------
#----------------------------------------------------------------------------------

obs_totalprod <- obs_totalprod %>%
  filter(bin_count >= 20)
grob_text <- grobTree(textGrob("Energy Equivalence: Slope = 0", x = 0.17, y = 0.88, hjust = 0,
                               gp = gpar(col = "gray52", fontsize = 20))) 
geom_size = 3.5

fg_mid = c("fg", "bin_midpoint")

prod_range <- obs_totalprod %>%
  filter(fg != "unclassified",bin_count >= 20) %>%
  group_by(across(all_of(fg_mid))) %>%
  summarize(
    min = min(bin_value),
    max = max(bin_value)
  )

prod_range <- binned_data %>%
  filter(fg != "unclassified", abundance >= 20) %>%
  group_by(across(all_of(fg_mid))) %>%
  summarize(
    min = min(production)/42.84,
    max = max(production)/42.84,
  )

#ggplot2::geom_errorbar(data = dens_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), width = 0) +
  
# Modifications:  changed plotting order, geom_size, geom colors
plot_totalprod <-function(year_to_plot = 1995, 
                           fg_names = c("fg1", "fg2", "fg3","fg4", "fg5", "all"), 
                           model_fit_density = 1, 
                           model_fit_production = 1, 
                           x_limits, 
                           x_breaks = c(1, 3, 10, 30, 100, 300), 
                           y_limits = c(0.03, 100), 
                           y_breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                           y_labels, 
                           fill_names = guild_fills_all, # c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                           color_names = guild_colors_all, #c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray"), 
                           x_name = "Diameter (cm)", 
                           y_name = expression(paste("Productivity (kg cm"^-1, " ha"^-1, "  yr"^-1, ")")),
                           geom_size = 4, 
                           obsdat = obs_totalprod, 
                           preddat = fitted_totalprod, 
                           plot_abline = TRUE, 
                           abline_slope = 0, 
                           dodge_width = 0.0,
                           abline_intercept = 20) 
{
  pos <-  ggplot2::position_dodge(width = dodge_width)
  obsdat <- obsdat %>% 
    dplyr::filter(fg %in% fg_names, year == year_to_plot, bin_count >= 20) %>% 
    dplyr::filter(bin_value > 0) %>% 
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1')))
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% 
    dplyr::summarize(min_obs = min(bin_midpoint), 
                     max_obs = max(bin_midpoint)) 
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% 
    dplyr::filter(dens_model %in% model_fit_density, 
                  prod_model %in% model_fit_production, 
                  fg %in% fg_names, year == year_to_plot) %>% 
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% 
    dplyr::filter(dbh >=  min_obs & dbh <= max_obs)
  
  p <- ggplot2::ggplot() + 
    ggplot2::geom_ribbon(data = preddat, 
                         aes(x = dbh, ymin = q025, ymax = q975, 
                             group = fg, fill = fg), alpha = 0.4, show.legend = F) + 
    ggplot2::geom_line(data = preddat, show.legend = F,
                       aes(x = dbh, y = q50, group = fg, color = fg)) + 
    ggplot2::geom_point(data = obsdat, position = pos,
                        aes(x = bin_midpoint, y = bin_value, group = fg, fill = fg), 
                        size = geom_size, color = "black", shape = 21) + 
    ggplot2::geom_errorbar(data = prod_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg), width = 0) +
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, 
                           limits = y_limits, breaks = y_breaks, labels = y_labels,  position = "left") + 
    ggplot2::scale_color_manual(values = guild_colors_all) + 
    ggplot2::scale_fill_manual(values = guild_fills_all) + 
    theme_plant() + #theme_no_x() +
    theme(aspect.ratio = 1)
  if (plot_abline) 
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", 
                                  linetype = "dashed", size = 0.75)
  p
} 
p <- plot_totalprod(year_to_plot = 1995,
                     fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                     model_fit_density = DENS, 
                     model_fit_production = PROD,
                     x_limits = c(0.9,160),
                     y_limits = c(0.3, 200),
                     y_breaks = c(0.1, 1, 10, 100),
                     y_labels = c(0.1, 1, 10, 100),
                     #dodge_width = 0.0,
                     preddat = fitted_totalprod)


#p1 <- set_panel_size(p, width=unit(14.3,"cm"), height=unit(14.3,"cm"))
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/New_main/production/Total_Production.pdf'))
grid.draw(p1)
dev.off()
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/production/Total_Production.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/production/Total_Production.pdf')) 
)

########################################################################
######################   Ratios     ####################################
########################################################################

#----------------------------------------------------------------------------------
#--------------------------   Richness Ratio -------------------------------------------
#----------------------------------------------------------------------------------

# Load 1995 bin data
load(file.path(gdrive_path, 'data/data_binned/bin_object_singleyear.RData'))

prod_ratio_diam <- breeder_stats_bydiam_byyear %>% 
  mutate(ID = 'Short-Tall') %>%
  rename(production_ratio = breeder_production_ratio, 
         density_ratio = breeder_density_ratio) %>%
  rbind(fastslow_stats_bydiam_byyear %>% 
          mutate(ID = 'Fast-Slow') %>%
          rename(production_ratio = fastslow_production_ratio, 
                 density_ratio = fastslow_density_ratio)) %>%
  filter(year == 1995)

prod_ratio_light <- breeder_stats_bylight_byyear %>% 
  mutate(ID = 'Short-Tall') %>%
  rename(production_ratio = breeder_production_ratio, 
         density_ratio = breeder_density_ratio) %>%
  rbind(fastslow_stats_bylight_byyear %>% 
          mutate(ID = 'Fast-Slow') %>%
          rename(production_ratio = fastslow_production_ratio, 
                 density_ratio = fastslow_density_ratio)) %>%
  filter(year == 1995)

# Load fitted values
ratio_fitted_diam <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/ratio_fittedvalues.csv'))
ratio_fitted_lightarea <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/ratio_fittedvalues_lightarea.csv'))

# Limited fitted lines to plotted data bins (min 20 pairs)
ratio_fitted_lightarea_prod <- ratio_fitted_lightarea %>%
  filter(variable == 'total production', light_area > 7, light_area < 197) %>%
  mutate(ratio = if_else(ratio == 'fast:slow', 'Fast-Slow', 'Short-Tall'))

ratio_fitted_diam_density <- ratio_fitted_diam %>%
  filter(variable == 'density', (ratio == 'fast:slow' & dbh < 66) | 
           (ratio == 'pioneer:breeder') & dbh < 16) %>%
  mutate(ratio = if_else(ratio == 'fast:slow', 'Fast-Slow', 'Short-Tall'))


obs_richnessbydiameter_ratio$tall_slow_ratio = obs_richnessbydiameter_ratio$richness_fg2/obs_richnessbydiameter_ratio$richness_fg3


# new ratios 
ratios <- read_csv(file.path(gdrive_path, 'data/data_binned/additional_bins_ratio_year.csv'))

#------------------ Plot Richness tall slow  diameter -----------
ratio_pionslow_range <- ratios %>%
  filter(min_n_individuals_pioneerslow >= 20) %>%
  group_by(bin_midpoint) %>%
  summarize(
    min = min(richness_ratio_pioneerslow),
    max = max(richness_ratio_pioneerslow)
  )

(rich_ratio_tallslow  <- #ggplot(data = obs_richnessbydiameter_ratio %>% 
                          #        filter(n_individuals_fg3 >= 20,n_individuals_fg2 >= 20),
                           #     aes(x =  bin_midpoint, y = tall_slow_ratio, fill = tall_slow_ratio)) + # exclude largest short:tall ratio
   ggplot() + 
    stat_smooth(data = ratios %>% 
                  filter(min_n_individuals_pioneerslow >= 20, year == "1995"),
                aes(x =  bin_midpoint, y = richness_ratio_pioneerslow, fill = richness_ratio_pioneerslow),
  method = "lm", color = "black", alpha = 0.2 ) +
    geom_point(data = ratios %>% 
                 filter(min_n_individuals_pioneerslow >= 20, year == "1995"),
               aes(x =  bin_midpoint, y = richness_ratio_pioneerslow, fill = richness_ratio_pioneerslow),
               shape = 21, stroke = 0.5, 
               size = 4, 
               #  size = 4, 
               color = "black") +
    scale_tall_slow  +
    scale_x_log10(limits = c(.9, 160), breaks = c(1, 10, 100), 
                  position = "bottom", 
                #  position = "top", 
                  expression(paste('Stem Diameter (cm)'))) + 
    scale_y_log10(labels = signif, breaks = c( 0.1, 0.3, 1, 3, 10), 
                  limits=c(0.3, 4),
                  #position = "left",
                  position = "right",
                  name = expression("Richness Ratio")) + 
    
    theme_plant() +theme_no_x()  +
    geom_errorbar(data = ratio_pionslow_range, aes(x = bin_midpoint, ymin = min, ymax = max), width = 0, color = "black") 

  
)


p1 <- set_panel_size(rich_ratio_tallslow , width=unit(10.25,"cm"), height=unit(8,"cm"))
#p1 <- set_panel_size(rich_ratio_tallslow , width=unit(7.6,"cm"), height=unit(7.6,"cm"))

grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path, 'Figures/New_main/ratios/richness_ratio_diam_tallslow.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/ratios/richness_ratio_diam_tallslow.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/ratios/richness_ratio_diam_tallslow.pdf')) 
)

#------------------ Plot Fast Slow  Richness ~ diameter -----------

(rich_ratio_fastslow  <- ggplot(data = obs_richnessbydiameter_ratio %>% 
                                  filter(n_individuals_fg3 >= 20,n_individuals_fg2 >= 20),
                                aes(x =  bin_midpoint, y = richness_ratio_fastslow, fill = richness_ratio_fastslow)) + # exclude largest short:tall ratio
    stat_smooth(method = "lm", color = "black", alpha = 0.2 ) +
    geom_point(shape = 21, stroke = 0.5, 
               size = 4, 
               #  size = 4, 
               color = "black") +
    scale_fast_slow3  +
    scale_x_log10(limits = c(.9, 160), breaks = c(1, 3, 10, 30, 100), 
                  #position = "bottom", 
                  position = "top", 
                  expression(paste('Stem Diameter (cm)'))) + 
    scale_y_log10(labels = signif, breaks = c( 0.1, 0.3, 1, 3, 10), 
                  limits=c(0.3, 4),
                  #position = "left",
                  position = "left",
                  name = expression("Richness Ratio")) + 
    
    theme_plant() + theme_no_x() +theme_no_y()
  
)


p1 <- set_panel_size(rich_ratio_fastslow, width=unit(10.25,"cm"), height=unit(8,"cm"))
#p1 <- set_panel_size(rich_ratio_fastslow  , width=unit(8,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path, '/Figures/New_main/ratios/richness_ratio_diam_fastslow.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/ratios/richness_ratio_diam_fastslow.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/ratios/richness_ratio_diam_fastslow.pdf')) 
)

#----------------------------------------------------------------------------------
#--------------------------  Abundance Ratio -------------------------------------------
#----------------------------------------------------------------------------------

obs_richnessbydiameter_ratio$tall_slow_abun_ratio = obs_richnessbydiameter_ratio$n_individuals_fg2/obs_richnessbydiameter_ratio$n_individuals_fg3
obs_richnessbydiameter_ratio$fast_slow_abun_ratio = obs_richnessbydiameter_ratio$n_individuals_fg1/obs_richnessbydiameter_ratio$n_individuals_fg3

#------------------------ Abundance tall slow ----------------
(abun_ratio_tallslow  <- ggplot(data = obs_richnessbydiameter_ratio %>% 
                                  filter(n_individuals_fg3 >= 20,n_individuals_fg2 >= 20),
                                aes(x =  bin_midpoint, y = tall_slow_abun_ratio, fill = tall_slow_abun_ratio)) + # exclude largest short:tall ratio
   #geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
   stat_smooth(method = "lm", color = "black", alpha = 0.2 ) +
   geom_point(shape = 21, stroke = 0.5, 
              size = 4, 
              color = "black") +
   scale_tall_slow  +
   scale_x_log10(limits = c(.9, 160), breaks = c(1, 10,  100), 
                 position = "bottom", 
                 expression(paste('Stem Diameter (cm)'))) + 
   scale_y_log10(labels = signif, 
                 breaks = c(0.03, 0.1, 1, 0.3, 1, 3), 
                 limits=c(0.03, 4),
                 position = "right",
                 name = expression("Abundance Ratio")) + 
   theme_plant() + theme_no_x()
 
)

p1 <- set_panel_size(abun_ratio_tallslow , width=unit(10.25,"cm"), height=unit(8,"cm"))
#p1 <- set_panel_size(abun_ratio_tallslow , width=unit(8,"cm"), height=unit(8,"cm"))

grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path, 'Figures/New_main/ratios/abun_ratio_diam_tall_slow_right.pdf'))
#pdf(file.path(gdrive_path, '/Figures/New_main/ratios/abun_ratio_diam_tall_slow_square.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/ratios/abun_ratio_diam_tall_slow_right.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/ratios/abun_ratio_diam_tall_slow_right.pdf')) 
)

#------------------------ Abundance Fast slow ----------------

(abun_ratio_fastslow  <- ggplot(data = obs_richnessbydiameter_ratio %>% 
                                  filter(n_individuals_fg3 >= 20, n_individuals_fg1 >= 20),
                                aes(x =  bin_midpoint, y = fast_slow_abun_ratio, fill = fast_slow_abun_ratio)) + # exclude largest short:tall ratio
    stat_smooth(method = "lm", color = "black", alpha = 0.2 ) +
    geom_point(shape = 21, stroke = 0.5, 
               size = 4, 
               #  size = 4.5, 
               color = "black") +
    scale_fast_slow3  +
    scale_x_log10(limits = c(.9, 160), breaks = c(1, 10,  100), 
                  #position = "bottom", 
                  position = "top", 
                  expression(paste('Stem Diameter (cm)'))) + 
    scale_y_log10(labels = signif,  breaks = c( 0.03, 0.3, 3), 
                  limits=c(0.03, 4),
                  position = "left",
                  #position = "right",
                  name = expression("Abundance Ratio")) + 
    
    theme_plant() + theme_no_x() +theme_no_y()
  
)


p1 <- set_panel_size(abun_ratio_fastslow, width=unit(10.25,"cm"), height=unit(8,"cm"))
#p1 <- set_panel_size(abun_ratio_fastslow, width=unit(8,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


pdf(file.path(gdrive_path, '/Figures/New_main/ratios/abun_ratio_diam_fast_slow.pdf'))
#pdf(file.path(gdrive_path, '/Figures/New_main/ratios/abun_ratio_diam_fast_slow_square.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/ratios/abun_ratio_diam_fast_slow.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/ratios/abun_ratio_diam_fast_slow.pdf')) 
)


#----------------------------------------------------------------------------------
#--------------------------  Productivity  Ratio -------------------------------------------
#----------------------------------------------------------------------------------

prod1995 = obs_totalprod %>%
  filter(year == 1995, bin_count >= 20) #%>%
#mutate(fast_slow_prod = NA, tall_slow_prod = NA)
fast_slow_prod = prod1995$bin_value[prod1995$fg == "fg1"]/prod1995$bin_value[prod1995$fg == "fg3"]
tall_slow_prod = prod1995$bin_value[prod1995$fg == "fg2"][1:15]/prod1995$bin_value[prod1995$fg == "fg3"]
prod_ratio = data.frame("bin_midpoint" = prod1995$bin_midpoint[1:15], "fast_slow_prod" = fast_slow_prod, "tall_slow_prod" = tall_slow_prod)   



(prod_fastslow  <- ggplot(data = prod_ratio, 
                          aes(x =  bin_midpoint, y = fast_slow_prod , fill = fast_slow_prod )) + # exclude largest short:tall ratio
    stat_smooth(method = "lm", color = "black", alpha = 0.2 ) +
    geom_point(shape = 21, stroke = 0.5, 
               size = 4, 
               color = "black") +
    scale_fast_slow3  +
    scale_x_log10(
      limits = c(.9, 160),
      position = "bottom", 
      expression(paste('Stem Diameter (cm)'))) + 
    scale_y_log10(labels = signif, breaks = c(0.03,  0.3, 3), 
                  limits=c(0.03, 4),
                  position = "right",
                  name = expression("Productivity Ratio")) + 
    
    theme_plant() #+theme_no_y()#+ theme_no_x() 
)


p1 <- set_panel_size(prod_fastslow  , width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)
#pdf(file.path(gdrive_path, '/Figures/New_main/ratios/production_ratio_diam_fastslow_square.pdf'))
pdf(file.path(gdrive_path, 'Figures/New_main/ratios/production_ratio_diam_fastslow.pdf'))
grid.draw(p1)
dev.off()
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/ratios/production_ratio_diam_fastslow.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/ratios/production_ratio_diam_fastslow.pdf')) 
)

#---------- productivity Tall vs Slow --------------
(prod_tallslow <- ggplot(data = prod_ratio, 
                         aes(x =  bin_midpoint, y = tall_slow_prod , fill = tall_slow_prod )) + # exclude largest short:tall ratio
    stat_smooth(method = "lm", color = "black", alpha = 0.2 ) +
    geom_point(shape = 21, stroke = 0.5, 
               size = 4, 
               #  size = 4, 
               color = "black") +
    scale_tall_slow  +
    scale_x_log10(#limits = c(1, 100), breaks = c(1, 3, 10, 30, 100), 
      limits = c(.9, 160),
      #position = NULL, 
      expression(paste('Stem Diameter (cm)'))) + 
    scale_y_log10(labels = signif, breaks = c( 0.03, 0.1, 0.3, 1, 3, 10), 
                  limits=c(0.03, 4),
                  #position = "left",
                  position = "right",
                  name = expression("Productivity Ratio")) + 
    theme_plant() + theme_no_x() +theme_no_y()
)

p1 <- set_panel_size(prod_tallslow , width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)

#pdf(file.path(gdrive_path, '/Figures/New_main/ratios/production_ratio_diam_tallslow_square.pdf'))
#pdf(file.path(gdrive_path, '/Figures/New_main/ratios/production_ratio_diam_tallslow.pdf'))
pdf(file.path(gdrive_path, 'Figures/New_main/ratios/production_ratio_diam_tallslow_right.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/ratios/production_ratio_diam_tallslow_right.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/ratios/production_ratio_diam_tallslow_right.pdf')) 
)

######################################################################################################
#########################         Mechanisms             ############################################
######################################################################################################


#---------------------------------------------------------------------------------------------------
#--------------------------     Diameter Growth         -------------------------------------------
#---------------------------------------------------------------------------------------------------
#used modified mass growth function


diam_growth_range <- obs_indivdiamgrowth %>%
  filter(!fg %in% c("unclassified", "all"), mean_n_individuals >= 20) %>%
  group_by(across(groups)) %>%
  summarize(
    min = min(mean),
    max = max(mean)
  )
diam_growth_range

plot_diam <- function (year_to_plot = 1995, 
                        fg_names = c("Fast", "Tall", "Slow",  "Short", "Medium"), 
                        model_fit = 1, 
                        x_limits, 
                        x_breaks = c(1, 3, 10, 30, 100, 300), 
                        y_limits, 
                        y_labels, 
                        y_breaks, 
                        fill_names = guild_fills_fg, # c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                        color_names = guild_colors_fg, #c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),
                        x_name = "Diameter (cm)", 
                        y_name = expression(paste("Diameter Growth (cm yr"^-1, ")")),
                        average = "mean", 
                        plot_errorbar = FALSE, 
                        error_min = "ci_min",
                        error_max = "ci_max", 
                        error_bar_width = 0.1,
                        error_bar_thickness = 0.5, 
                        dodge_width = 0.07, 
                        dodge_errorbar = TRUE, 
                        geom_size = 4, 
                        obsdat = obs_indivprod, 
                        preddat = fitted_indivprod, 
                        plot_abline = TRUE, 
                        abline_slope = 2, 
                        abline_intercept = -1.25) {
  pos <- if (dodge_errorbar) 
    ggplot2::position_dodge(width = dodge_width)
  else "identity"
  obsdat <- obsdat %>% 
    dplyr::filter(fg %in% fg_names, year ==  year_to_plot, 
                  !is.na(mean), mean_n_individuals >= 20) %>% 
    dplyr::group_by(bin_midpoint) %>% 
    dplyr::mutate(width = error_bar_width *  dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    arrange(desc(fg))
  obs_limits <- obsdat %>% 
    dplyr::group_by(fg) %>% 
    dplyr::summarize(min_obs = min(bin_midpoint), 
                     max_obs = max(bin_midpoint))
  preddat <- preddat %>% 
    dplyr::left_join(obs_limits) %>% 
    dplyr::filter(prod_model %in% model_fit, 
                  fg %in% fg_names, year == year_to_plot) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% 
    dplyr::filter(dbh >=  min_obs & dbh <= max_obs) %>% 
    arrange(desc(fg))
  p <- ggplot2::ggplot() + 
    ggplot2::geom_ribbon(data = preddat %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))), 
                         ggplot2::aes(x = dbh, ymin = q025, ymax = q975, 
                                      group = fg, fill = fg), alpha = 0.4) + 
    ggplot2::geom_line(data = preddat %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))), 
                       ggplot2::aes(x = dbh, y = q50, group = fg, color = fg))
  if (plot_errorbar) {
    p <- p + ggplot2::geom_errorbar(data = obsdat %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))), 
                                    ggplot2::aes_string(x = "bin_midpoint", 
                                                        ymin = error_min, 
                                                        ymax = error_max, 
                                                        group = "fg", 
                                                        color = "fg", 
                                                        width = "width"), 
                                    position = pos, 
                                    size = error_bar_thickness)
  }
  p <- p + ggplot2::geom_line(data = preddat[preddat$fg == "Medium", ], 
                              ggplot2::aes(x = dbh, y = q50), color = "gray") + 
    ggplot2::geom_point(data = obsdat %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))),
                        ggplot2::aes_string(x = "bin_midpoint", 
                                            y = average, group = "fg", fill = "fg"), 
                        size = geom_size, color = "black", shape = 21, position = pos) + 
    ggplot2::geom_errorbar(data = diam_growth_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg), width = 0) +
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) + 
    #theme_no_x() + 
    ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")) + 
    ggplot2::scale_color_manual(values = color_names) + 
    ggplot2::scale_fill_manual(values = fill_names) + 
    theme_plant() +#+ theme_no_x() +
    if (plot_abline) {
      p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                    slope = abline_slope, color = "gray72", linetype = "dashed", 
                                    size = 0.75)
    }
  p
}
p
#change to stem diameter growth
p <- plot_diam(year_to_plot = 1995,
                fg_names = c('fg1','fg2','fg3','fg4','fg5'),
                model_fit = PROD,
                x_limits = c(.9, 160),
                x_breaks = c(1, 10, 100),
                y_limits = c(0.02, 1.3),
                y_breaks = c(0.03, 0.1, 0.3, 1),
                y_labels = c(0.03, 0.1, 0.3, 1),
                error_bar_width = 0,
                dodge_width = 0.00,
                obsdat = obs_indivdiamgrowth, #diameter growth
                preddat = fitted_indivdiamgrowth,
                plot_abline = FALSE,
                x_name = 'Stem Diameter (cm)',
                y_name = expression(paste('Diameter Growth (cm yr'^-1,')')))
p
(p1 <- p + theme(axis.text.x = element_text(), axis.ticks.x = element_line()) + 
    annotation_custom(grob_fast) + annotation_custom(grob_tall) + annotation_custom(grob_medium) + 
    annotation_custom(grob_slow) + annotation_custom(grob_short) +
    labs(x = 'Stem Diameter (cm)') + theme(plot.margin = grid::unit(c(0,1,0,0), "cm")) +theme_plant_small())


p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p2)
pdf(file.path(gdrive_path,'Figures/New_main/growth_diam/Diam_growth_x.pdf'))
grid.draw(p2)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/growth_diam/Diam_growth_x.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/growth_diam/Diam_growth_x.pdf')) 
)

##################################################
# new quentin piecewise
# Script for additional mortality versus diameter fit
# Hinged piecewise functional form
# Each nonlinear parameter (except hinge smoothness parameter) has functional group as fixed effect
# 21 Sep 2022


# Setup -------------------------------------------------------------------

library(dplyr)
library(readr)
library(brms)
library(tidybayes)
library(ggplot2)
library(tidyr)
library(forestscaling)
library(purrr)

options(mc.cores = 4, brms.backend = 'cmdstanr', brms.file_refit = 'on_change')

# Load data
load('data/rawdataobj1995.RData')

diam_data <- alltreedat[[3]] %>%
  filter(!is.na(fg), !recruit) %>%
  select(fg, dbh_corr, diam_growth_rate) %>%
  mutate(fg = paste0('fg', fg))

# Fit model ---------------------------------------------------------------

diam_hinge_fixef_fit <- brm(
  bf(
    log(diam_growth_rate) ~ beta0 + beta1low * (log(dbh_corr) - log(x0)) + (beta1high - beta1low) * delta * log(1 + exp((log(dbh_corr) - log(x0)) / delta)),
    beta0 ~ 0 + fg,
    beta1low ~ 0 + fg,
    beta1high ~ 0 + fg,
    x0 ~ 0 + fg,
    delta ~ 1,
    nl = TRUE
  ),
  data = diam_data, family = gaussian(link = 'identity'),
  prior = c(
    prior(normal(0, 2), nlpar = beta0),
    prior(lognormal(1, 1), nlpar = beta1low, lb = 0),
    prior(lognormal(1, 1), nlpar = beta1high, lb = 0),
    prior(lognormal(1, 1), nlpar = x0, lb = 0),
    prior(exponential(10), nlpar = delta, lb = 0)
  ),
  chains = 4, iter = 3000, warmup = 2000, seed = 27603,
  file = '~/temp/forestlight/diam_hinge_fixef_brmfit'
)

# Postprocessing and plotting  -------------------------------------

# Set number of bins
numbins <- 20

# Remove trees not belonging to any functional group
alltreedat_classified <- map(alltreedat, ~ filter(., !is.na(fg)))

# Bin classified trees. (log binning of density)
allyeardbh_classified <- map(alltreedat_classified[-1], ~ pull(., dbh_corr)) %>% unlist
dbhbin_allclassified <- logbin(x = allyeardbh_classified, y = NULL, n = numbins)
bin_edges <- c(dbhbin_allclassified$bin_min,dbhbin_allclassified$bin_max[numbins])

qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
diam_bins <- diam_data %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = bin_edges, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  summarize(p = qprobs, q = quantile(diam_growth_rate, probs = qprobs), n = n()) %>%
  pivot_wider(names_from = p, values_from = q) %>%
  setNames(c('fg', 'dbh_bin', 'n', 'q025', 'q25', 'median', 'q75', 'q975')) %>%
  mutate(bin_midpoint = dbhbin_allclassified$bin_midpoint[as.numeric(dbh_bin)])


# Prediction grid: dbh x fg
# Limit ranges to the observed data points and bins with > 20 individuals
diam_max <- diam_bins %>%
  filter(n >= 20) %>%
  group_by(fg) %>%
  filter(bin_midpoint == max(bin_midpoint))

pred_dat <- diam_max %>%
  group_by(fg) %>%
  group_modify(~ data.frame(dbh_corr = exp(seq(log(1), log(.$bin_midpoint), length.out = 101))))

# Get expected values of posterior means (epred) for each combination of diameter and functional group
diamhinge_pred <- pred_dat %>%
  add_epred_draws(diam_hinge_fixef_fit) %>%
  mutate(.epred = exp(.epred))

### Quick diag. plot
ggplot(diam_bins %>% filter(n >= 20), aes(x=bin_midpoint, y = median)) +
  stat_lineribbon(aes(y = .epred, x = dbh_corr), data = diamhinge_pred, size = 0.5) +
  geom_pointrange(aes(ymin = q25, ymax = q75), color = 'gray40') +
  facet_wrap(~ fg) + scale_x_log10(name = 'diameter (cm)') + scale_y_log10(name = 'diameter growth (cm/y)') +
  scale_fill_brewer(palette = 'Blues') +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2), strip.background = element_blank())


# Generate parameter tables -----------------------------------------------

diam_hinge_params <- diam_hinge_fixef_fit %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  separate(.variable, into = c('b', 'parameter', 'fg')) %>%
  select(-b) %>%
  mutate(fg = as.integer(substr(fg, 5, 5))) 

diam_hinge_quant <- diam_hinge_params %>%
  group_by(parameter, fg) %>%
  median_qi(.width = c(0.5, 0.9, 0.95)) %>%
  pivot_wider(id_cols = c(parameter, fg, .value), names_from = .width, values_from = c(.lower, .upper)) %>%
  select(fg, parameter, .value, .lower_0.95, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9, .upper_0.95) %>%
  setNames(c('fg', 'parameter', 'median', 'q025', 'q05', 'q25', 'q75', 'q95', 'q975')) %>%
  mutate(parameter = factor(parameter, levels = c('beta0', 'beta1low', 'beta1high', 'x0', 'delta'))) %>%
  arrange(parameter, fg) %>%
  mutate(fg = c(rep(c('fast', 'large pioneer', 'slow', 'small breeder', 'medium'), 4), '(none)'),
         parameter = c(rep(c('beta0 (intercept)', 'beta1low (slope at small size)', 'beta1high (slope at large size)', 'x0 (cutoff between small and large sizes)'), each = 5), 'delta (smoothing parameter for hinge)'))

# Generate plotting data --------------------------------------------------

diam_hinge_plotting_data <- diamhinge_pred %>%
  group_by(fg, dbh_corr) %>%
  median_qi(.epred, .width = c(0.5, 0.9, 0.95)) %>%
  pivot_wider(id_cols = c(fg, dbh_corr, .epred), names_from = .width, values_from = c(.lower, .upper)) %>%
  select(fg, dbh_corr, .lower_0.95, .lower_0.9, .lower_0.5, .epred, .upper_0.5, .upper_0.9, .upper_0.95) %>%
  setNames(c('fg', 'dbh', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975')) 

#write_csv(diam_hinge_quant, 'data/clean_summary_tables/clean_diamgrowthbydiameter_parameters.csv')
#write_csv(diam_hinge_plotting_data, 'data/data_forplotting/fitted_diamgrowthbydiameter.csv')

#---------------------------------------------------------------------------------------------------
#--------------------------   Mass Growth  --------------------------------------------------------
#---------------------------------------------------------------------------------------------------

mass_growth_range <- obs_indivprod %>%
  filter(!fg %in% c("unclassified", "all"), mean_n_individuals >= 20) %>%
  group_by(across(groups )) %>%
  summarize(
    min = min(mean),
    max = max(mean)
  )
mass_growth_range
ggplot2::geom_errorbar(data = dens_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), width = 0) +
  
plot_growth <- 
  function (year_to_plot = 1995, 
            fg_names = c("Fast", "Tall", "Slow",  "Short", "Medium", "All"), 
            model_fit = 1, 
            x_limits, 
            x_breaks = c(1, 3, 10, 30, 100, 300), 
            y_limits, 
            y_labels, 
            y_breaks, 
            fill_names = guild_fills, # c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
            color_names = guild_colors, #c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),
            x_name = "Diameter (cm)", 
            y_name = expression(paste("Mass Growth (kg yr"^-1, ")")),
            average = "mean", 
            position = "left",
            plot_errorbar = FALSE, 
            error_min = "ci_min",
            error_max = "ci_max", 
            error_bar_width = 0.1,
            error_bar_thickness = 0.5, 
            dodge_width = 0.07, 
            dodge_errorbar = TRUE, 
            geom_size = 4, 
            obsdat = obs_indivprod, 
            preddat = fitted_indivprod, 
            plot_abline = TRUE, 
            abline_slope = 2, 
            abline_intercept = -1.25) {
    pos <- if (dodge_errorbar) 
      ggplot2::position_dodge(width = dodge_width)
    else "identity"
    obsdat <- obsdat %>% 
      dplyr::filter(fg %in% fg_names, year ==  year_to_plot, 
                    !is.na(mean), mean_n_individuals >= 20) %>% 
      dplyr::group_by(bin_midpoint) %>% 
      dplyr::mutate(width = error_bar_width *  dplyr::n()) %>% 
      dplyr::ungroup() %>% 
      arrange(desc(fg))
    obs_limits <- obsdat %>% 
      dplyr::group_by(fg) %>% 
      dplyr::summarize(min_obs = min(bin_midpoint), 
                       max_obs = max(bin_midpoint))
    preddat <- preddat %>% 
      dplyr::left_join(obs_limits) %>% 
      dplyr::filter(prod_model %in% model_fit, 
                    fg %in% fg_names, year == year_to_plot) %>% 
      dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                       dplyr::all_vars(. > min(y_limits))) %>% 
      dplyr::filter(dbh >=  min_obs & dbh <= max_obs) %>% 
      arrange(desc(fg))
    p <- ggplot2::ggplot() + 
      ggplot2::geom_ribbon(data = preddat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), 
                           ggplot2::aes(x = dbh, ymin = q025, ymax = q975, 
                                        group = fg, fill = fg), alpha = 0.4) + 
      ggplot2::geom_line(data = preddat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), 
                         ggplot2::aes(x = dbh, y = q50, group = fg, color = fg))
    if (plot_errorbar) {
      p <- p + ggplot2::geom_errorbar(data = obsdat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), 
                                      ggplot2::aes_string(x = "bin_midpoint", 
                                                          ymin = error_min, 
                                                          ymax = error_max, 
                                                          group = "fg", 
                                                          color = "fg", 
                                                          width = "width"), 
                                      position = pos, 
                                      size = error_bar_thickness)
    }
    p <- p + ggplot2::geom_line(data = preddat[preddat$fg == "Medium", ], 
                                ggplot2::aes(x = dbh, y = q50), color = "gray") + 
      ggplot2::geom_point(data = obsdat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))),
                          ggplot2::aes_string(x = "bin_midpoint", 
                                              y = average, group = "fg", fill = "fg"), 
                          size = geom_size, color = "black", shape = 21, position = pos) + 
      ggplot2::geom_errorbar(data = mass_growth_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), width = 0) +
      ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
      ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels, position = position) + 
      #theme_no_x() + 
      ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")) + 
      ggplot2::scale_color_manual(values = color_names) + 
      ggplot2::scale_fill_manual(values = fill_names) + 
      theme_plant() #+ theme_no_x()
    if (plot_abline) {
      p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                    slope = abline_slope, color = "gray72", linetype = "dashed", 
                                    size = 0.75)
    }
    p
  }

obs_indivprod <- obs_indivprod %>%
  filter(mean_n_individuals >= 20)

p <- plot_growth(year_to_plot = 1995,
                fg_names = c('fg1','fg2','fg3','fg4','fg5'),
                model_fit = PROD,
                x_limits = c(0.9, 160),
                y_limits = c(0.003, 160),
                y_breaks = c(0.01, 0.1, 1, 10, 100),
                plot_errorbar = F,
                error_min = 'q25',
                error_max = 'q75',
                error_bar_width = 0,
                y_labels = c( 0.01, 0.1, 1, 10, 100),
                abline_slope = 2, 
                abline_intercept = 10, #-1.4
                dodge_width = 0.0) #0.07


p0 <- p #+ annotation_custom(grob_text_a) 
p1 <- set_panel_size(p0, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/New_main/growth_mass/mass_growth.pdf'))
grid.draw(p1)
dev.off()
system2(command = "pdfcrop", 
        args  = c(file.path(gdrive_path2,'Figures/New_main/growth_mass/mass_growth.pdf'), 
                  file.path(gdrive_path2,'Figures/New_main/growth_mass/mass_growth.pdf')) 
)


#---------------------------------------------------------------------------------------------------
#--------------------------   Mortalty   --------------------------------------------------------
#---------------------------------------------------------------------------------------------------


bin_mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_mortalitybins.csv')) # Load binned data
mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_mortalityindividuals.csv')) # Load raw data
fitted_mort <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/mortality_ci_by_fg.csv')) # Load fitted values
#growth_diam <-  read_csv(file.path(gdrive_path, 'data/clean_summary_tables/clean_parameters_individualdiametergrowth.csv')) 

bin_mort <- bin_mort %>% arrange(factor(fg, levels = c('all', 'fg5','fg4', 'fg3', 'fg2','fg1')))
fitted_mort <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/mortality_ci_by_fg.csv')) # Load fitted values

#---  Plot mortality 
obs_range_mort <- bin_mort %>% 
  arrange(factor(fg, levels = c('fg5','fg4', 'fg3', 'fg2','fg1'))) %>%
  filter(variable %in% 'dbh', lived + died >= 20) %>%
  group_by(fg) %>%
  dplyr::summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))

fitted_mort_trunc <- fitted_mort %>%
  arrange(factor(fg, levels = c('fg5','fg4', 'fg3', 'fg2','fg1'))) %>%
  left_join(obs_range_mort) %>%
  left_join(guild_lookup) #%>%
#filter(dbh >= min_obs & dbh <= max_obs)
unique(fitted_mort_trunc$fg)

(p <- ggplot() +#data = fitted_mort_trunc %>% mutate(fg = factor(fg, labels = fg_labels))) +
    #geom_ribbon(aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
    #geom_line(aes(x = dbh, y = q50, group = fg, color = fg)) +
    
     stat_smooth(span = 2, data = bin_mort %>%  # to single out gray and make CI darker
                filter(variable == 'dbh', fg == "fg5", 
                     (lived+died) >= 20),
            aes(x = bin_midpoint, y = mortality, se = T),
               color = "gray", fill = "gray", alpha = 0.2) +
    stat_smooth(span = 2, data = bin_mort %>% 
                  filter(variable == 'dbh', !fg %in% c('all','unclassified'), 
                         (lived+died) >= 20)  %>% 
                  mutate(fg = factor(fg, labels = fg_labels)),
                aes(x = bin_midpoint, y = mortality, color = fg, fill = fg), alpha = 0.25) +
  
    scale_y_continuous(trans = "logit", position = "right", 
                       breaks = c(0.03, 0.1, 0.3, 0.6), 
                       labels =  c(0.03, 0.1, 0.3, 0.6), 
                       limits = c(0.02, 0.8),
                       name = expression(paste("Mortality (5 yr"^-1,")"))) +
    geom_point(data = bin_mort %>% 
                 filter(variable == 'dbh', !fg %in% c('all','unclassified'), 
                        (lived+died) >= 20)  %>% 
                 mutate(fg = factor(fg, labels = fg_labels)),
               aes(x = bin_midpoint, y = mortality, fill = fg),
               shape = 21, size = 4) +
    scale_x_log10(name = parse(text = 'Stem~Diameter~(cm)'), 
                  # breaks = c(3, 30, 300), 
                  limits = c(.9, 160)
    ) +
    scale_color_manual(values = guild_colors) +
    scale_fill_manual(values = guild_fills) +
    theme_plant() )#+ theme_no_x())


p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p2)  #notice that LLP and Slow show steep slow of mortality with light!
pdf(file.path(gdrive_path,'Figures/New_main/mortality/mortality.pdf'))
grid.newpage()
grid.draw(p2)  
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/mortality/mortality.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/mortality/mortality.pdf')) 
)


#---------------------------------------------------------------------------------------------------
#--------------------------   Richness ~ Abundance   --------------------------------------------------------
#---------------------------------------------------------------------------------------------------

#------------- Richness vs abundance

max_size_fg <- obs_richnessbydiameter %>%
  filter(n_individuals > 0) %>%
  group_by(fg) %>%
  dplyr::summarize(max_size = max(abundance_by_bin_width)) 

min_size_fg <- obs_richnessbydiameter %>%
  filter(n_individuals > 0) %>%
  group_by(fg) %>%
  dplyr::summarize(min_size = min(abundance_by_bin_width)) 

fitted_richnessvsabundance_filt <- fitted_richnessvsabundance %>%
  left_join(max_size_fg, by = "fg") %>%
  left_join(min_size_fg, by = "fg") %>%
  group_by(fg) %>% 
  filter(abundance_by_bin_width >= min_size) %>%
  filter(abundance_by_bin_width <= max_size)


# per ha
plot_area = 50
(rich_abun <- ggplot() + 
    theme_plant() +
    scale_fill_manual(values = guild_fills) +
    scale_color_manual(values = guild_colors) +
    geom_jitter(data = obs_richnessbydiameter %>%   arrange(desc(fg)) %>% 
                  filter(n_individuals > 0) %>%
                  filter(!fg %in% c('unclassified', 'all') & richness > 0),  #  & n_individuals >= 20
                aes(x = abundance_by_bin_width/plot_area, y = richness_by_bin_width/plot_area,
                    fill = fg, color = fg),
                shape = 21, size = 4, color = "black", width = 0)  +
    theme(axis.title.y = element_text(vjust = -3)) + 
    scale_x_log10(name = expression(paste("Abundance (cm"^-1," ha"^-1,")")),
                  breaks = c(0.01,  1,  100),
                  labels = c("0.01",  "1", "100"),
                  #labels = signif,
                  limit = c(0.001, 300)
    ) + 
    scale_y_log10(labels = signif,
                  limit = c(0.0003, 3), 
                  position = "right",
                  name = expression(paste("Richness (cm"^-1," ha"^-1,")"))) +
    geom_ribbon(data = fitted_richnessvsabundance_filt %>% 
                  filter(!fg %in% c('unclassified', 'all')) %>%
                  arrange(desc(fg)),
                aes(x = abundance_by_bin_width/plot_area, ymin = q025/plot_area, ymax = q975/plot_area, fill = fg, color = NA), alpha = 0.2, col = NA) + 
    geom_line(data = fitted_richnessvsabundance_filt %>% 
                filter(!fg %in% c('unclassified', 'all')) %>%
                arrange(desc(fg)),
              aes(x = abundance_by_bin_width/plot_area, y = q50/plot_area, color = fg), size = 0.5) )

#+#fill = c("fg5" = "gray93", "all" = "black", "fg4" =  "#87Cefa", "fg3" = "#27408b", "fg2" = "#267038", "fg1" = "#BFE046"),
#color = c("fg5" = "gray", "all" = "black", "fg4" =  "#87Cefa", "fg3" = "#27408b", "fg2" = "#267038", "fg1" = "#BFE046")),


p1 <- set_panel_size(rich_abun, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/New_main/rich_abun/rich_abun_per_ha.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/New_main/rich_abun/rich_abun_per_ha.pdf'), 
                    file.path(gdrive_path2,'Figures/New_main/rich_abun/rich_abun_per_ha.pdf')) 
)


