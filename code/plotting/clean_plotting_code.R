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

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users/jgradym/Google_Drive/ForestLight'))
github_path <- ifelse(Sys.info()['user'] == 'qread', '~/documents/github/', file.path('/Users/jgradym/Documents/Github'))

library(broom)
library(forestscaling) # Packaged all the functions and ggplot2 themes here!
library(tidyverse)
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
#install_github('qdread/forestscaling')

# Define color schemes and labels
guild_fills <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")
guild_fills2 <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")
guild_colors <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")
guild_colors2 <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray")

fg_labels <- c("fg1", "fg2","fg3", "fg4", "fg5")
fg_labels2 <- c("all", "fg1", "fg2","fg3", "fg4", "fg5")
fg_labels3 <- c("all", "fg1", "fg2","fg3", "fg4", "fg5", "unclassified")

guild_names <- paste('fg', 1:5, sep = '')
guild_labels <- c('Fast','Tall', 'Slow', 'Short', 'Medium')
guild_labels2 <- c('All', 'Fast','Tall', 'Slow', 'Short', 'Medium')
guild_labels3 <- c('All', 'Fast','Tall', 'Slow', 'Short', 'Medium',  "Unclassified")

guild_lookup <- data.frame(fg = c('fg1','fg2','fg3','fg4','fg5','all','unclassified'), 
                           fg_name = c('Fast','Tall','Slow','Short','Medium','All','Unclassified'))

year_to_plot = 1995
geom_size <- 4

guide <- guides(fill = guide_legend(title = NULL), override.aes = list(fill = NA), color = F)
guide2 <- guides(color = guide_legend(title = NULL))


grob_a <- grobTree(textGrob("a", x = 0.04, y = 0.93,  hjust = 0,
                            gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob_b <- grobTree(textGrob("b", x = 0.04, y = 0.93,  hjust = 0,
                            gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 

#guide <- guides(color = guide_legend(title=NULL), fill = guide_legend(title=NULL)) #, color = F, override.aes=list(fill=NA))

# To add back the legend
theme_plant2 <- theme_plant() + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
theme_plant <- theme_plant() + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position = "right", legend.text = element_text(size = 14 ), 
        legend.key = element_blank())


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
grob_slow3 <- grobTree(textGrob("Slow", x = 0.04, y = 0.81,  hjust = 0,
                               gp = gpar(col = "#27408b", fontsize = 15, fontface = "italic"))) 
grob_short3 <- grobTree(textGrob("Short", x = 0.04, y = 0.74,  hjust = 0,
                                gp = gpar(col = "#87Cefa", fontsize = 15, fontface = "italic"))) 

# Quentin, I wanted to add annotation_grob() these below outside of the main plot in lieu of a legend but ggplot won't let you.
# found something at https://stackoverflow.com/questions/12409960/ggplot2-annotate-outside-of-plot but it screws with the border and proportions
#if you know a better fix I'd love to hear it

grob_fast2 <- grobTree(textGrob("Fast", x = 1.1, y = 0.95,  hjust = 0,
                               gp = gpar(col = "#BFE046", fontsize = 13, fontface = "italic"))) 
grob_tall2 <- grobTree(textGrob("Tall", x = 1.1, y = 0.88,  hjust = 0,
                               gp = gpar(col = "#267038", fontsize = 13, fontface = "italic"))) 
grob_medium2 <- grobTree(textGrob("Medium", x = 1.1, y = 0.81,  hjust = 0,
                                 gp = gpar(col = "gray70", fontsize = 13, fontface = "italic"))) 
grob_slow2 <- grobTree(textGrob("Slow", x = 1.1, y = 0.74,  hjust = 0,
                               gp = gpar(col = "#27408b", fontsize = 13, fontface = "italic"))) 
grob_short2 <- grobTree(textGrob("Short", x = 1.1, y = 0.67,  hjust = 0,
                                gp = gpar(col = "#87Cefa", fontsize = 13, fontface = "italic"))) 
grob_all2 <- grobTree(textGrob("All", x = 1.1, y = 0.60,  hjust = 0,
                              gp = gpar(col = "black", fontsize = 13, fontface = "italic"))) 

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

# source the extra extraction functions that aren't in the package
source(file.path(github_path, 'forestscalingworkflow/R_functions/model_output_extraction_functions.r'))
source(file.path(github_path, 'forestlight/stan/get_ratio_slopes_fromfit.R'))


################################################################################################
# ------------------------------ Fig 1: Light Interception ---------------------------------
################################################################################################

lightperareacloudbin_fg <- read.csv(file.path(fp, 'lightperareacloudbin_fg.csv'), stringsAsFactors = FALSE)
lightpervolcloudbin_fg <- read.csv(file.path(fp, 'lightpervolcloudbin_fg.csv'), stringsAsFactors = FALSE)
unscaledlightbydbhcloudbin_fg <- read.csv(file.path(fp, 'unscaledlightbydbhcloudbin_fg.csv'), stringsAsFactors = FALSE)


exl <- expression(atop('Light per Crown Area', paste('(W m'^-2, ')')))
exv <- expression(atop('Light per Crown Volume', paste('(W m'^-3, ')')))
exd <- 'Diameter (cm)'


#----------------------   Fig 1a: Light per crown area by diameter -----------------------------

p <- ggplot() + geom_point(alpha = 0.01, data = alltree_light_95, 
                           aes(x = dbh_corr, y = light_received_byarea), color = 'chartreuse3') +
  geom_pointrange(data = lightperareacloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = mean, ymin = q25, ymax = q75)) +
  scale_x_log10(limits = c(0.8, 300), name = exd) +
  scale_y_log10(name = exl, limits = c(0.8, 1000)) +
  theme_plant_small() 


# --- to add secondary height axis, first mask theme_no_x() above
area <- p + scale_x_log10(name = 'Diameter (cm)',
                        breaks = c(1,10,100), limits = c(0.8, 300),
                        sec.axis = sec_axis(~ gMM(., a = 57.17, b = 0.7278, k = 21.57),
                                            name = "Height (m)", breaks = c(3, 10, 30))) 
area

p1 <- set_panel_size(area, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path,'Figures/Light_Individual/light_area.pdf'))
grid.draw(p1)
dev.off()

#code to crop out white space - may only work on unix/macs
system2(command = "pdfcrop", 
        args  = c(file.path(gdrive_path,'Figures/Light_Individual/light_area.pdf'), 
                  file.path(gdrive_path,'Figures/Light_Individual/light_area.pdf')) 
)
lm1 <- lm(log(light_received_byarea) ~ log(dbh_corr), data = alltree_light_95)
summary(lm1)
confint(lm1)

# using linear fit to estimate increase up to 100 cm dbh
exp(0.987*log(100)) #94.2
#using ratio of smallest and largest bin
412.197295/9.477105 #43
#----------------------   Fig 1b: Light per crown volume by diameter -----------------------------

vol <- ggplot() +
  geom_point(alpha = 0.01, data = alltree_light_95, 
             aes(x = dbh_corr, y = light_received_byvolume), color = 'chartreuse3') +
  geom_pointrange(data = lightpervolcloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = mean, ymin = q25, ymax = q75)) +
  scale_x_log10(limits = c(0.8, 300), name = exd) +
  #stat_smooth(method = "lm", color = "black", 
   #           data = alltree_light_95, aes(x = dbh_corr, y = light_received_byvolume)) +
  scale_y_log10(name = exv, breaks = c(1, 10, 100), limits = c(0.8, 200)) +
  theme_plant_small()

vol
p1 <- set_panel_size(vol, width=unit(10.25,"cm"), height=unit(7,"cm"))

grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path,'Figures/Light_Individual/light_volume.pdf'))
grid.draw(p1)
dev.off()
system2(command = "pdfcrop", 
        args  = c(file.path(gdrive_path,'Figures/Light_Individual/light_volume.pdf'), 
                  file.path(gdrive_path,'Figures/Light_Individual/light_volume.pdf')) 
)

lm1 <- lm(log(light_received_byvolume) ~ log(dbh_corr), data = alltree_light_95)
summary(lm1)
confint(lm1)
lm2 <- lm(log(light_received_byvolume) ~ log(dbh_corr), data = alltree_light_95)
summary(lm2)
confint(lm2)
lm3 <- lm(log(light_received_byvolume) ~ log(dbh_corr), data = alltree_light_95)
summary(lm3)
confint(lm3)
exp(0.448*log(100)) #7.9
# or nonlinear
nlm1 <- nls(log(light_received_byvolume) ~ a * exp(b*log(dbh_corr)), 
            start = list(a = 1, b = 0.1), data = alltree_light_95) 
summary(nlm1)

#ratio of largest points for volume
35.6/8.5 #4.2
#combine

g_area  <- ggplotGrob(area)
g_vol<- ggplotGrob(vol)
g1 <- rbind(g_area , g_vol , size = "first")
g1$widths <- unit.pmax(g_area$widths, g_vol$widths)
grid.newpage()
grid.draw(g1)
ggsave(g1, height = 8, width = 11, filename = file.path(gdrive_path,'Figures/Light_Individual/light_combo.pdf'))

system2(command = "pdfcrop", 
        args  = c(file.path(gdrive_path,'Figures/Light_Individual/light_combo.pdf'), 
                  file.path(gdrive_path,'Figures/Light_Individual/light_combo.pdf')) 
)


################################################################################################
# ---------------------------- Fig 2: Scaling Schematic -------------------------
################################################################################################


################################################################################################
# ---------------------------- Fig 3: Life Histories Classification and Rates -------------------------
################################################################################################


# Load Nadja's data (new functional groups 25 June)
# fg5 is the new column (we originally used fg from the older df)
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
fgbci[fgbci$genus == "Luehea",]
fgbci[fgbci$genus == "Cecropia",]
#anacardium excelsum  = slow
#Luehea seemannii = slow


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
      #lab_x <- expression(paste(italic('Slow'), 'to', italic('Fast')))


####--------- Fig 3a, PCA
Fig_3a <- ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, fill = factor(fg5))) +
  geom_point(shape = 21, size = geom_size, color = "black") + 
  labs(x = 'Survivorship–Growth Tradeoff', y = 'Stature—Recruitment Tradeoff') +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group')+
  scale_fill_manual(values = guild_fills) + theme_plant()
Fig_3a 

p2  <- set_panel_size(Fig_3a , width=unit(14.3,"cm"), height=unit(14.3,"cm"))
grid.newpage()
grid.draw(p2)
pdf(file.path(gdrive_path, 'Figures/Life_History/LH_a.pdf'))
grid.draw(p2)
dev.off()


###-------  Fig 3b Growth by Light

title_x <- expression(paste('Light per Crown Area (W m'^-2,')',sep=''))
title_y <- expression(atop('Growth per Crown Area', paste('(kg yr'^-1, ' m'^-2,')', sep='')))
scale_y_log10(name =  expression(atop('Growth per Crown Area',
                                      paste('(kg y'^-1, ' m'^-2,')'))))


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

dodge_width <- 0.07
error_bar_width <- 0.04


# Do some additional computation to correct the error bar width for the number of groups in each bin
obs_light_binned_plotdata <- obs_light_binned %>% 
  arrange(factor(fg, levels = c('all', 'fg5','fg4','fg2', 'fg3', 'fg1'))) %>%
  filter(year == year_to_plot, mean_n_individuals >= 20, !fg %in% c('alltree', 'unclassified')) %>%
  group_by(bin_midpoint, year) %>% 
  mutate(width = sum(c('fg1','fg2','fg3','fg4','fg5') %in% fg)) %>% 
  ungroup

#----- Plot Growth vs Light
Fig_3b  <- ggplot(obs_light_binned_plotdata) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
              aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.4) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot),
            aes(x = light_area, y = q50, color = fg)) +
  geom_errorbar(aes(x = bin_midpoint, ymin = q25, ymax = q75, 
                    group = fg, color = fg, width = 0), #width = error_bar_width * width), 
                position = position_dodge(width = dodge_width)) + 
  geom_point(aes(x = bin_midpoint, y = mean, group = fg, fill = fg),
             size = 4, shape = 21, position = position_dodge(width = dodge_width)) +
  scale_x_log10(name = title_x, limits = c(1.5, 412)) + 
  scale_y_log10(name = title_y, position = "right", breaks = c(0.01, 0.03, 0.1, 0.3), 
                labels = c( 0.01, 0.03, 0.1, 0.3)) +
  scale_color_manual(name = 'Functional group', values = guild_fills, labels = fg_labels) +
  scale_fill_manual(values = guild_fills, labels = fg_labels, guide = FALSE) +
  theme_no_x() + theme_plant2
Fig_3b 

p2 <- set_panel_size(Fig_3b , width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p2)

pdf(file.path(gdrive_path, "Figures/Life_History/LH_b.pdf"))
grid.draw(p2)
dev.off()


#----------------------- Fig 3c: Mortality vs Light -----

bin_mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_mortalitybins.csv')) # Load binned data
mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_mortalityindividuals.csv')) # Load raw data
fitted_mort <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/mortality_ci_by_fg.csv')) # Load fitted values
growth_diam <-  read_csv(file.path(gdrive_path, 'data/clean_summary_tables/clean_parameters_individualdiametergrowth.csv')) 

# Truncate fitted mortality lines to not exceed the range of the observed data, using 20 individuals as the cutoff.
bin_mort <- bin_mort %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg2', 'fg3', 'fg1')))
obs_range_mort <- bin_mort %>% 
  arrange(factor(fg, levels = c('fg5','fg4','fg2', 'fg3', 'fg1'))) %>%
  filter(variable %in% 'light_per_area', lived + died >= 20) %>%
  group_by(fg) %>%
  summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))

fitted_mort_trunc <- fitted_mort %>%
  arrange(factor(fg, levels = c('fg5','fg4','fg2', 'fg3', 'fg1'))) %>%
  left_join(obs_range_mort) %>%
  left_join(guild_lookup) %>%
  filter(light_per_area >= min_obs & light_per_area <= max_obs)
unique(fitted_mort_trunc$fg)

#---  Plot mortality vs light
p <- ggplot(data = fitted_mort_trunc %>% mutate(fg = factor(fg, labels = fg_labels))) +
  geom_ribbon(aes(x = light_per_area, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.4) +
  geom_line(aes(x = light_per_area, y = q50, group = fg, color = fg)) +
  geom_point(data = bin_mort %>% 
               filter(variable == 'light_per_area', !fg %in% c('all','unclassified'), 
                      (lived+died) >= 20)  %>% 
               mutate(fg = factor(fg, labels = fg_labels)),
             aes(x = bin_midpoint, y = mortality, fill = fg),
             shape = 21, size = geom_size) +
  scale_x_log10(name = parse(text = 'Light~per~Crown~Area~(W~m^-2)'), breaks = c(3, 30, 300), limits = c(1.5, 412)) +
  scale_y_continuous(trans = "logit", position = "right", breaks = c(0.03, 0.1, 0.3, 0.6), 
                     labels = c(0.03, 0.1, 0.3, 0.6), limits = c(0.02, 0.65),
                name = expression(paste("Mortality (5 yr"^-1,")"))) +
  scale_color_manual(values = guild_colors) +
  scale_fill_manual(values = guild_fills) +
  theme_plant()
p

p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p2)  #notice that LLP and Slow show steep slow of mortality with light!

pdf(file.path(gdrive_path, "Figures/Life History/LH_c.pdf"))
grid.draw(p2)
dev.off()


################################################################################################
# ------------------------ Fig 4: Main Scaling Plots --------------------------------
################################################################################################

# Plot of slopes in different segments by different functional groups.


# Create plots.
#Model fit 1 = pareto, 1 segment
#Model Fit 2  = 2 segments, etc

grob_text_a <- grobTree(textGrob("a", x = 0.06, y = 0.92, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
grob_text_b <- grobTree(textGrob("b", x = 0.06, y = 0.88, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
grob_text_c <- grobTree(textGrob("c", x = 0.05, y = 0.95, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
geom_size = 2


# Specify dodging with a certain width of error bar
# Model fit 1 = power law
# Model fit 2 = power law exp


# ------------------- Fig 4a: Individual Growth ~ Diameter --------------------------

plot_prod2 <- 
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
            y_name = expression(paste("Growth (kg yr"^-1, ")")),
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
            abline_intercept = -1.3) {
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
      ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
      ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) + 
      #theme_no_x() + 
      ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")) + 
      ggplot2::scale_color_manual(values = color_names) + 
      ggplot2::scale_fill_manual(values = fill_names) + 
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
                x_limits = c(0.8, 200),
                y_limits = c(0.001, 2000),
                y_breaks = c(0.001,0.1, 10, 1000),
                plot_errorbar = T,
                error_min = 'q25',
                error_max = 'q75',
                error_bar_width = 0,
                y_labels = c(0.001,0.1,10,1000),
                dodge_width = 0.07)
p

p0 <- p #+ annotation_custom(grob_text_a) 
p1 <- set_panel_size(p0, width=unit(10.25,"cm"), height=unit(7,"cm"))

grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path,'Figures/Main_Scaling/Growth.pdf'))
grid.draw(p1)
dev.off()

# --- to add secondary height axis, first mask theme_no_x() above
p0 <- p + scale_x_log10(name = 'Diameter (cm)',
                        breaks = c(1,10,100), limits = c(0.8, 230),
                        sec.axis = sec_axis(~ gMM(., a = 57.17, b = 0.7278, k = 21.57),
                                            name = "Height (m)", breaks = c(3, 10, 30))) #+
 
p0


p1 <- set_panel_size(p0, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path,'Figures/Main_Scaling/Growth_Height.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
                  args  = c(file.path(gdrive_path,'Figures/Main_Scaling/Growth_Height.pdf'), 
                            file.path(gdrive_path,'Figures/Main_Scaling/Growth_Height.pdf')) 
           )
# ------------------- Fig 4b: Density ~ Diameter --------------------------

obs_dens <- obs_dens %>%
  filter(bin_count >= 20)
plot_dens2 <- function (year_to_plot = 1995, 
                        fg_names = c("fg1", "fg2", "fg3", "fg4", "fg5", "all"),
                        model_fit = 1, 
                        x_limits, 
                        x_breaks = c(1, 3, 10, 30, 100, 300),
                        y_limits, 
                        y_breaks, 
                        y_labels, 
                        fill_names = guild_fills2, #c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                        color_names = guild_colors2, #c("black","#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),  
                        x_name = "Diameter (cm)",
                        y_name = expression(paste("Density (n ha"^-1, "cm"^-1, ")")), 
                        geom_size = 4, 
                        obsdat = obs_dens, 
                        preddat = pred_dens, 
                        plot_abline = TRUE, 
                        abline_slope = -2, 
                        abline_intercept = 4,
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
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) + 
    ggplot2::scale_color_manual(values = color_names) + 
    ggplot2::scale_fill_manual(values = fill_names) + 
    theme_plant()
  }
p <- plot_dens2(year_to_plot = 1995,
          fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
          model_fit = DENS,
          dodge_width = 0.0,
          x_limits = c(.8, 230),
          y_limits = c(0.007, 20000),
          x_breaks = c(1, 10, 100),
          y_labels = c(0.001, 0.1, 10,1000),
          y_breaks = c(0.001, 0.1,  10, 1000))
#p <- p +annotation_custom(grob_text_b)

p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)
pdf(file.path(gdrive_path,'Figures/Main_Scaling/Density.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Main_Scaling/Density.pdf'), 
                    file.path(gdrive_path,'Figures/Main_Scaling/Density.pdf')) 
)

#------------------- Fig 4C: Total Production ~ Diameter --------------------------
obs_totalprod <- obs_totalprod %>%
  filter(bin_count >= 20)
grob_text <- grobTree(textGrob("Energy Equivalence: Slope = 0", x = 0.17, y = 0.88, hjust = 0,
                               gp = gpar(col = "gray52", fontsize = 20))) 

# Modifications: added jitter, changed plotting order, geom_size, geom colors
plot_totalprod2 <-function(year_to_plot = 1995, 
                           fg_names = c("fg1", "fg2", "fg3","fg4", "fg5", "all"), 
                           model_fit_density = 1, 
                           model_fit_production = 1, 
                           x_limits, 
                           x_breaks = c(1, 3, 10, 30, 100, 300), 
                           y_limits = c(0.03, 100), 
                           y_breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                           y_labels, 
                           fill_names = guild_fills2, # c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                           color_names = guild_colors2, #c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray"), 
                           x_name = "Diameter (cm)", 
                           y_name = expression(paste("Production (kg cm"^-1, " ha"^-1, "  yr"^-1, ")")),
                           geom_size = 4.5, 
                           obsdat = obs_totalprod, 
                           preddat = fitted_totalprod, 
                           plot_abline = TRUE, 
                           abline_slope = 0, 
                           dodge_width = 0.0,
                           abline_intercept = 2) 
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
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, 
                           limits = y_limits, breaks = y_breaks, labels = y_labels,  position = "right") + 
    ggplot2::scale_color_manual(values = guild_colors2) + 
    ggplot2::scale_fill_manual(values = guild_fills2) + 
    theme_plant() + 
    theme(aspect.ratio = 1)
  if (plot_abline) 
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", 
                                  linetype = "dashed", size = 0.75)
  p
} 
p <- plot_totalprod2(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
               model_fit_density = DENS, 
               model_fit_production = PROD,
               x_limits = c(0.7,230),
               y_limits = c(0.5, 200),
               y_breaks = c(0.1, 1, 10, 100),
               y_labels = c(0.1, 1, 10, 100),
               #dodge_width = 0.0,
               preddat = fitted_totalprod)
p

p1 <- set_panel_size(p, width=unit(14.3,"cm"), height=unit(14.3,"cm"))

grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/Main_Scaling/Total_Production.pdf'))
grid.draw(p1)
dev.off()

# ---------------- Add height sec axis

#---------- Alternative to main plot
p0 <- p + scale_x_log10(name = 'Diameter (cm)',
                        breaks = c(1, 3, 10, 30,100), limits = c(0.7, 230),
                        sec.axis = sec_axis(~ gMM(., a = 57.17, b = 0.7278, k = 21.57),
                                            name = "Height (m)", 
                                            breaks = c(2, 3, 5, 10, 20, 30, 40))) +
  scale_y_log10(position = "right", limits = c(0.5, 200), 
                breaks = c(0.1, 1, 10, 100),
                labels = c(0.1, 1, 10, 100), 
                name = expression(paste("Production (kg cm"^-1, " ha"^-1, "  yr"^-1, ")"))) +
  annotation_custom(grob_text)
p0
p1 <- set_panel_size(p0, width=unit(14.3,"cm"), height=unit(14.3,"cm"))
grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/Main_Scaling/Total_Production_Height.pdf'))
grid.draw(p1)
dev.off()


# ------- Supplemental Height Secondary axis
p0 <- p + scale_x_log10(name = 'Diameter (cm)',
                        breaks = c(1,10,100), limits = c(0.7, 230),
                        sec.axis = sec_axis(~ gMM(., a = 57.17, b = 0.7278, k = 21.57),
                                            name = "Height (m)", 
                                            breaks = c(2, 3, 5, 10, 20, 30, 40))) +
  scale_y_log10(position = "left", limits = c(0.5, 200), 
                breaks = c(0.1, 1, 10, 100),
                labels = c(0.1, 1, 10, 100),
                name =expression(atop('Production', paste('(kg yr'^-1,' cm'^-1,' ha'^-1,')')))) +
  theme(aspect.ratio = 0.8) + geom_point(size = 1)

p0

p1 <- set_panel_size(p0, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/Supplementals/Total_Production_Height/Total_prod_height.pdf'))
grid.draw(p1)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Total_Production_Height/Total_prod_height.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Total_Production_Height/Total_prod_height.pdf')) 
)



########################################################################################
# ------------------------------- Fig 5 Symmetry Plots ------------------------------------
########################################################################################


fitted_indivlight <- read.csv(file.path(fp, 'fitted_indivlight.csv'), stringsAsFactors = FALSE)
fitted_totallight <- read.csv(file.path(fp, 'fitted_totallight.csv'), stringsAsFactors = FALSE)
indivlightbins_fg <- read.csv(file.path(fp, 'obs_indivlight.csv'), stringsAsFactors = FALSE)
totallightbins_fg <- read.csv(file.path(fp, 'obs_totallight.csv'), stringsAsFactors = FALSE)

# Fig 5a      
totallightbins_fg <- totallightbins_fg %>%
  filter(bin_count >= 20)
p <- plot_totalprod(year_to_plot = 1995,
                         fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                         model_fit_density = DENS, 
                         model_fit_production = PROD,
                         x_limits = c(0.9,150),
                         y_limits = c(100, 200000),
                         geom_size = 4,
                         y_breaks = c(100, 1000, 10000, 100000),
                         y_labels = c("0.1", "1", "10", "100"),
                         y_name = expression(paste('Total Light Intercepted (kW cm'^-1,' ha'^-1,')')),
                         preddat = fitted_totallight,
                         obsdat = totallightbins_fg,
                         plot_abline = FALSE)
p2 <- p + theme_plant()
p2
# ---------------Fig 5b--------------

params <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/piecewise_paramci_by_fg.csv'), stringsAsFactors = FALSE)

xvalues <- params %>% 
  filter(variable == 'density', model == DENS) %>%
  group_by(fg) %>%
  summarize(mid_cutoff = (mean[parameter == 'tau_high'] + mean[parameter == 'tau_low'])/2)


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

fg_full_names <- c('Fast', 'Tall', 'Slow', 'Short', 'Medium', 'All', 'Unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

allslopes <- rbind(growth_slopes_atmiddle, light_slopes_atmiddle) %>%
  ungroup %>%
  mutate(fg = factor(fg, levels = fgs, labels = fg_full_names))
grob0 <- grobTree(textGrob("b", x = 0.04, y = 0.9,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob1 <- grobTree(textGrob("Solar", x = 0.68, y = 0.94, hjust = 0,
                           gp = gpar(col = "gold3", fontsize = 18))) 
grob2 <- grobTree(textGrob("Production", x = 0.68, y = 0.86, hjust = 0,
                           gp = gpar(col = "darkgreen", fontsize = 18)))
grob3 <- grobTree(textGrob("Energy Equivalence", x = 0.28, y = 0.51, hjust = 0,
                           gp = gpar(col = "black", fontsize = 18))) 


# Plot
grob_text <- grobTree(textGrob("Solar Equivalence", x = 0.28, y = 0.87, hjust = 0,
                               gp = gpar(col = "gold3", fontsize = 18))) 

grob_text2 <- grobTree(textGrob("a", x = 0.06, y = 0.91, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
totallightbins_fg <- totallightbins_fg %>%
  filter(bin_count >= 20)
tot_light <- plot_totalprod(year_to_plot = 1995,
                            fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                            model_fit_density = DENS, 
                            model_fit_production = PROD,
                            x_limits = c(0.9,150),
                            y_limits = c(100, 200000),
                            geom_size = 3.5,
                            y_breaks = c(100, 1000, 10000, 100000),
                            y_labels = c("0.1", "1", "10", "100"),
                            y_name = expression(paste('Total Light Intercepted (kW cm'^-1,' ha'^-1,')')),
                            preddat = fitted_totallight,
                            obsdat = totallightbins_fg,
                            plot_abline = FALSE)


tot_light2 <- tot_light  + 
  scale_y_continuous(position = "left", trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("0.1", "1", "10", "100"), 
                     limits = c(200, 200000),
                     name = expression(atop('Total Light Intercepted',paste('(kW cm'^-1,' ha'^-1,')'))))  +
  theme(aspect.ratio = 0.75) + 
  geom_abline(intercept = log10(70000), slope = 0, color ="#C9B074",
              linetype="dashed", size=.75) +
  annotation_custom(grob_text) + annotation_custom(grob_text2)
plot(tot_light2)
p_tot_light <- tot_light2


g_tot_light <- ggplotGrob(p_tot_light)

#---------------------Fig 5b--------------
# compare slopes
slopes <- ggplot(allslopes %>% filter(!fg %in% 'Unclassified'), 
                 aes(x = fg, y = q50, ymin = q025, ymax = q975, fill =  variable, color =variable)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = .75) +
  geom_errorbar(position = position_dodge(width = 0.6), size = 0.75, width = 0) +
  geom_point( position = position_dodge(width = 0.6), shape = 21, size = 3.5, color = "black", stroke = 0.5) +
  geom_errorbar(data = allslopes %>% filter(fg %in% c('Tall', 'Slow', 'All')), 
                position = position_dodge(width = 0.6), size = 0.75, width = 0) +
  labs( x = NULL, y = 'Scaling Slope') +
  scale_y_continuous(limits = c(-1.05, 1.3)) +
  scale_fill_manual(values = c('gold1', 'darkolivegreen3')) +
  scale_color_manual(values = c('gold3', 'darkgreen')) +
  theme_plant() + theme(axis.text.x = element_text(angle = 25,  vjust = .7, face = "italic", size = 18)) +
  annotation_custom(grob1) + annotation_custom(grob2) + annotation_custom(grob3) +annotation_custom(grob0)
slopes

g_slopes <- ggplotGrob(slopes)

combo <- rbind(g_tot_light, g_slopes, size = "first")
combo$widths <- unit.pmax(g_tot_light$widths,g_slopes$widths)
grid.newpage()
grid.draw(combo)
ggsave(combo, height = 8.6, width = 6, filename = file.path(gdrive_path,'Figures/Symmetry/combo.pdf'))

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Symmetry/combo.pdf'), 
                    file.path(gdrive_path,'Figures/Symmetry/combo.pdf')) 
)


########################################################################################
# ------------------------------- Fig 6, Ratio Scaling ------------------------------------
########################################################################################
# ------------------------------- Fig 6 Relative Abundance & Production ------------------------------------

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

# ---------------------------------   Fig 6a  Production by light ------------------------------------

prod_ratio <- prod_ratio_light   %>%
  filter(n_individuals >= 20) %>%
  ggplot() + #aes(x = bin_midpoint, y = production_ratio, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_ribbon(aes(x = light_area, ymin = q025, ymax = q975, group = ratio, fill = ratio),
              alpha = 0.3, data = ratio_fitted_lightarea_prod) +
  geom_line(aes(x = light_area, y = q50, group = ratio,  color = ratio), 
            data = ratio_fitted_lightarea_prod) +
  geom_point(aes(x = bin_midpoint, y = production_ratio, fill = ID), 
             shape = 21, size = 4,  stroke = .5, color = "black") +
  scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey"))+
  scale_color_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey50")) +
  theme_plant_small() + 
  scale_x_log10(name = expression(paste('Light per Crown Area (W m'^-2,')')), 
                limits = c(2,300), breaks = c(3, 30, 300)) +
  scale_y_log10(labels = signif, breaks = c(0.01,0.1, 1,10,100,1000), 
                limits = c(0.01,400), position = "left",
                name = "Ratio") 
prod_ratio

prod_ratio2 <- set_panel_size(prod_ratio, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(prod_ratio2)

pdf(file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/Prod_ratio.pdf'))
grid.draw(prod_ratio2)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/Prod_ratio.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/Prod_ratio.pdf')) 
)
#----------------------------- Fig 6B Abundance by Diameter ----------------------------
dens_ratio <- prod_ratio_diam %>% 
  filter(density_ratio > 0) %>%
  filter(n_individuals >= 20) %>%
  ggplot() +
  geom_ribbon(aes(x = dbh, ymin = q025, ymax = q975, group = ratio, fill = ratio), 
              alpha = 0.4, data = ratio_fitted_diam_density) +
  geom_line(aes(x = dbh, y = q50, group = ratio, color = ratio), 
            data = ratio_fitted_diam_density) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_point(aes(x = bin_midpoint, y = density_ratio, fill = ID),
             shape = 21, size = 4,  stroke = .5,  color = "black") +
  scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey")) +
  scale_color_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey50")) +
  scale_x_log10(limits = c(1,100), breaks = c(1,10, 100), 
                name = expression(paste('Diameter (cm)'))) + 
  scale_y_log10(labels = signif, breaks = c(0.01,0.1, 1,10,100,1000), 
                limits=c(0.01,200),
                name = expression("Density Ratio")) + 
  theme_plant_small() +theme_no_y()
dens_ratio 

# ---------- combine

g_dens  <- ggplotGrob(dens_ratio )
g_prod <- ggplotGrob(prod_ratio)
g6 <- cbind(g_prod, g_dens , size = "first")
g6$heights <- unit.pmax(g_dens$heights, g_prod $heights)
grid.newpage()
grid.draw(g6)

ggsave(g6, height = 8, width = 11, filename = file.path(gdrive_path,'Figures/Ratios/Ratio_no_line.pdf'))
ggsave(g6, height = 8, width = 11, filename = file.path(gdrive_path,'Figures/Ratios/Ratio2.pdf'))

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Ratios/Ratio.pdf'), 
                    file.path(gdrive_path,'Figures/Ratios/Ratio.pdf')) 
)

#for visual comparison
g_dens  <- ggplotGrob(dens_ratio )
g_prod <- ggplotGrob(prod_ratio)
g6 <- rbind(g_prod, g_dens , size = "first")
g6$widths <- unit.pmax(g_dens$widths, g_prod $widths)
grid.newpage()
grid.draw(g6)

ggsave(g6, height = 8, width = 11, filename = file.path(gdrive_path,'Figures/Supplementals/Ratios_PCA/Main_comparison.pdf'))

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Ratios_PCA/Main_comparison.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Ratios_PCA/Main_comparison.pdf')) 
)



#---------------------------------------------------------------------------------------------
# ------------------------ Supplements  ---------------------------- 
#---------------------------------------------------------------------------------------------


#---------------Fig S1:  LAI, Tranmittance and Crown Depth
lai_0 <- read_csv('/Users/jgradym/Google_Drive/ForestLight/data/data_forplotting/LAI_Depth.csv')
pfd_0 <- read_csv('/Users/jgradym/Google_Drive/ForestLight/data/data_forplotting/PFD_LAI.csv')
lai <- lai_0 %>% filter(Species != "Cecropia")
pfd <- pfd_0 %>% filter(Species != "Cecropia")

fill_sp <- c("darkolivegreen4", "cadetblue2",  "brown3", "gray22") 
color_sp <- c("darkolivegreen", "cadetblue3",  "brown4", "black") 

#------------Fig S1A: T LAI vs Crown Depth
depth <- ggplot(data = lai, 
            aes(x = Depth, y = LAI, fill = Species, fill = Species, color = Species )) +
  geom_smooth(method = "lm", alpha = 0.2, 
              aes(fill = Species, color = Species),
              show.legend = FALSE) +
  geom_point(aes(fill = Species), size = 4.5, shape = 21, color = "black") + 
    theme_plant  +
  scale_y_log10(limits = c(0.6, 10),labels = signif, breaks = c(1, 3, 10)) +
  scale_x_log10(limits = c(0.2, 15), breaks = c( 0.3, 1, 3, 10), 
                labels = signif, name = "Crown Depth (m)") +
  scale_color_manual(values = color_sp) +
  scale_fill_manual(values = fill_sp) + guide +
  annotation_custom(grob_a) 
depth


pdf(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_depth.pdf'))
grid.draw(depth)
dev.off()

depth2 <- set_panel_size(depth, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(depth2)
pdf(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_depth.pdf'))
grid.draw(depth2)
dev.off()
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_depth.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/LAI/lai_depth.pdf')) 
)

# Mixed model slopes
lmer1 <- lmer(log(LAI) ~ log(Depth)+ (log(Depth) |Species), data = lai)
summary(lmer1) # 0.33
performance::r2(lmer1)

#try to run bayesian style
# weird error
b_lmer <- stan_lmer(log(LAI) ~ log(Depth) + (log(Depth) |Species), 
                    data = lai, seed = 1000)

blm <- stan_glm(log(LAI) ~ log(Depth) + Species, data = lai)
library(mlmRev)
library(lme4)
b_lmer <- stan_lmer(formula = LAI ~ Depth + (Depth |Species), 
                    data = lai, seed = 100)

data(Gcsemv, package = "mlmRev")
GCSE <- subset(x = Gcsemv, 
               select = c(school, student, gender, course))

M1_stanlmer <- stan_lmer(formula = course ~ 1 + (1 | school), 
                         data = GCSE,
                         seed = 1000)

model_code <- 'parameters {real y;} model {y ~ normal(0,1);}'
fit <- stan(model_code=model_code)  # this will take a minute

# calculated separately
lm2 <- lm(log(LAI) ~ log(Depth), data = lai)
summary(lm2)
lm_each <-  lai %>%
  nest(-Species) %>% #group variable
  mutate(
    fit = map(data, ~ lm(log(LAI) ~ log(Depth), data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)%>%
  filter(term != '(Intercept)')
lm_each
#average
mean(lm_each$estimate) #0.37


#------------Fig S1B: Transmittance vs LAI 

trans <- ggplot(data = pfd %>% filter(Species != "Cecropia"),
            aes( x = LAI, y = PFD, fill = Species, color = Species )) +
  geom_smooth(method = "lm", alpha = 0.3, 
              aes(fill = Species, color = Species), show.legend = F) +
  geom_point(aes(fill = Species), size = 4.5, shape = 21, color = "black") + 
  theme_plant  +
  scale_y_log10(limits = c(2.5, 200), name = "% PFD Transmittance") +
  scale_x_continuous(limits = c(-0.3, 8)) +
  scale_color_manual(values = color_sp) +
  scale_fill_manual(values = fill_sp) + guide +
  annotation_custom(grob_b) +
  theme(axis.title = element_text(margin = margin(t = 1, r = 0, b = 0, l = 0, unit = "cm")))
trans

pdf(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_pfd.pdf'))
trans
dev.off()

trans2 <- set_panel_size(trans, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(trans2 )
pdf(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_pfd.pdf'))
grid.draw(trans2 )
dev.off()
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_pfd.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/LAI/lai_pfd.pdf')) 
)

g_depth  <- ggplotGrob(depth)
g_trans <- ggplotGrob(trans)
g1 <- rbind(g_depth, g_trans, size = "first")
g1$widths <- unit.pmax(g_depth $widths, g_trans$widths)
grid.newpage()
grid.draw(g1)
ggsave(g1, height = 8, width = 11, filename = file.path(gdrive_path,'Figures/Light_Individual/light_combo.pdf'))


# all 4 species
mlm_pdf <- lmer(log(PFD) ~ LAI + (LAI |Species), pfd) 
mlm_pdf # 0.49 sloep



# pfd all
lm_pfd <- lm(log(PFD) ~ LAI, pfd)
lm_pfd
lm2_pfd <- lm(log(PFD) ~ LAI + Species, pfd)
lm2_pfd


lm_each_pfd <-  pfd %>%
  nest(-Species) %>% #group variable
  mutate(
    fit = map(data, ~ lm(log(PFD) ~ LAI, data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)%>%
  filter(term != '(Intercept)')
lm_each_pfd
mean(lm_each_pfd$estimate) #0.49



# ------------------------ Fig S2: Adapted from Fig 2a and Fig 1a in Ruger 2018 --------------------


# ------------------------ Fig S3: Diameter growth Scaling -------------------------

p <- plot_prod2(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3','fg4','fg5'),
               model_fit = PROD,
               x_limits = c(1, 230),
               x_breaks = c(1, 10, 100),
               y_limits = c(0.02, 1.3),
               y_breaks = c(0.03, 0.1, 0.3, 1),
               y_labels = c(0.03, 0.1, 0.3, 1),
               error_bar_width = 0,
               dodge_width = 0.07,
               obsdat = obs_indivdiamgrowth,
               preddat = fitted_indivdiamgrowth,
               plot_abline = FALSE,
               x_name = 'Diameter (cm)',
               y_name = expression(paste('Diameter growth (cm yr'^-1,')')))
               
p1 <- p + theme(axis.text.x = element_text(), axis.ticks.x = element_line()) + 
  annotation_custom(grob_fast) + annotation_custom(grob_tall) + annotation_custom(grob_medium) + 
  annotation_custom(grob_slow) + annotation_custom(grob_short) +
  labs(x = 'Diameter (cm)') + theme(plot.margin = grid::unit(c(0,1,0,0), "cm")) +theme_plant_small()


p2 <- set_panel_size(p1, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p2)
pdf(file.path(gdrive_path,'Figures/Supplementals/Diameter_Growth/Diam_growth.pdf'))
grid.draw(p2)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Diameter_Growth/Diam_growth.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Diameter_Growth/Diam_growth.pdf')) 
)


#-------------------------------------- ------------------- -------------------  
#-------------------   Fig S4: Growth Light Heat Map    -------------------
#------------------- ------------------- ------------------- ------------------- 

fg_labeler <- c("fg1" = "Fast", "fg2" = "Tall", "fg3" = "Slow", "fg4" = "Short", "fg5" = "Medium")
theme_facet2 <- function () 
{
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = NA), 
                 panel.border = ggplot2::element_rect(color = "black", 
                                                      fill = NA, size = 0.75),
                 legend.position = "none", 
                 panel.background = ggplot2::element_blank(), 
                 strip.text.x = ggplot2::element_text(size = 13), 
                 axis.text = ggplot2::element_text(size = 13, color = "black"), 
                 axis.ticks.length = ggplot2::unit(0.2, "cm"),
                 axis.title = ggplot2::element_text(size = 13))
}
# ------------------- growth Light Heat Map ------------------- 


#obs_light_raw$fg <- factor(obs_light_raw$fg, levels = c(paste0('fg', 1:5), 'unclassified'),
 #                          labels = c("Fast", "Pioneer", "Slow", "Breeder", "Medium", "Unclassified"))
#pred_light_5groups$fg <- factor(pred_light_5groups$fg, levels = c(paste0('fg', 1:5)),
 #                               labels = c("Fast", "Pioneer", "Slow", "Breeder", "Medium"))
hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000), 
                                             labels = c(1,10,100,1000), limits=c(1,5000))

unique(obs_light_raw$fg)
growth_light_hex<- ggplot(obs_light_raw %>% filter(year == year_to_plot, fg %nin% c('alltree','unclassified'))) +
  facet_wrap(~ fg, ncol = 2, labeller = as_labeller(fg_labeler)) +
  geom_hex(aes(x = light_area, y = production_area)) +
  #geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
  #          aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.3) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), 
            aes(x = light_area, y = q50), size = 1, color = 'black') +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, labels=signif) +
  hex_scale_log_colors + 
  theme_plant_small() +
  theme(strip.text = element_text(size=14))+
  guides(color = FALSE) +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA),
        legend.position = c(0.7, 0.15),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=15))
growth_light_hex


pdf(file.path(gdrive_path, "Figures/Supplementals/Growth_light/growth_light_hex.pdf"))
growth_light_hex
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Growth_light/growth_light_hex.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Growth_light/growth_light_hex.pdf')) 
)
#-----------------------------------------------------------------------------
#------------------------ Fig S5: Mean Growth with Light
#------------------------------------------------------------------------------

growth_light_mean <- ggplot(obs_light_binned %>% 
                              filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, ncol = 2, labeller = as_labeller(fg_labeler)) +
  geom_ribbon(data = pred_light_5groups %>% 
                filter(year == year_to_plot), 
              aes(x = light_area, ymin = q025, ymax = q975, group = fg, color=NA, fill = fg),
              alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% 
              filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), 
            aes(x = light_area, y = q50, group = fg, color = fg), size = 0.5) +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q25, yend = q75), size = 0.3) +
  geom_segment(data = cast_pars %>% 
                 filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), 
               aes(x = x_max_q50 * 0.5, xend = x_max_q50 * 2, y = y_max_q50 * 0.5, yend = y_max_q50 * 2), 
               color = 'brown1', size = .5) +
  geom_point(shape=21, aes(x = bin_midpoint, y = mean)) +
  scale_color_manual(values = guild_fills ) +
  scale_fill_manual(values = guild_colors) +
  scale_x_log10(name = title_x, breaks = c(1,10,100)) + 
  scale_y_log10(name = title_y, labels = signif, breaks = c(0.01,0.1,1,10)) +
  theme_plant_small() +
  theme(strip.text = element_text(size=14))+
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA))

growth_light_mean 

pdf(file.path(gdrive_path, "Figures/Supplementals/Growth_light/mean_growth_light_max.pdf"))
p_mean_panels
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Growth_light/mean_growth_light_max.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Growth_light/mean_growth_light_max.pdf')) 
)


#---------------------------- Raw growth by light + fg --------------------------------
p_raw_panels <- ggplot(obs_light_raw %>% 
                         filter(year == year_to_plot, !fg %in% c('alltree','unclassified'))) +
  facet_wrap(~ fg, ncol = 2, labeller = as_labeller(fg_labeler)) +
  geom_point(shape = 21,  alpha = 0.1, aes(x = light_area, y = production_area)) + #alpha = fg
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
              aes(x = light_area, ymin = q025, ymax = q975, group=fg,color=NA,fill=fg), alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), 
            aes(x = light_area, y = q50, group = fg, color = fg), size= 1) +
  #scale_alpha_manual(values = c(0.15, 0.15, 0.05, 0.008, 0.008)) +
  scale_alpha_manual(values = c(0.15, 0.15, 0.15, 0.15, 0.15)) +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, breaks=c(0.001,0.01,1),labels=signif) +
  scale_color_manual(values = guild_colors) +
  scale_fill_manual(values = guild_colors) +
  theme_plant_small() + theme_facet2()

p_raw_panels
pdf(file.path(gdrive_path, "Figures/Supplementals/Growth_light/growth_light_fg_raw_alph0.1.pdf"))
p_raw_panels
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Growth_light/growth_light_fg_raw_alph0.1.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Growth_light/growth_light_fg_raw_alph0.1.pdf')) 
)

#-------------------------------------------------------------------------
#---------------------- Fig S6 Max Growth Responsiveness with Light ------
#-------------------------------------------------------------------------

# 1. Plot of maximum slope by functional group

# Remove all tree and unclassified groups
param_ci$fg <- factor(param_ci$fg ,levels = c("fg1", "fg2", "fg5", "fg4", "fg3"))
guild_labels2_ <- c("Fast", "Tall", "Medium", "Short", "Slow")
growth_slope <- ggplot(param_ci %>% filter(fg != 'NA', year == year_to_plot, parameter %in% 'log_slope', !fg %in% c('alltree','unclassified')),
       aes(x = fg, y = mean, ymin = q025, ymax = q975)) + 
  geom_errorbar(width = 0) + geom_point(size = 4) +
  theme(axis.text.x = element_text(angle = 25,  vjust = 0.7))+
  scale_x_discrete(name = 'Life History Guild', labels = guild_labels2_) +
  scale_y_continuous(expression(atop('Max. Growth Responsiveness',paste('to Light (kg yr'^-1, ' W'^-1,')'))), 
                     limits = c(0.65, 1.1),
                     breaks = seq(0, 1.1, 0.1), labels = seq(0, 1.1, 0.1)) +
  theme_plant_small() + theme(aspect.ratio = 0.75) + annotation_custom(grob_b)
growth_slope 
pdf(file.path(gdrive_path, "Figures/Supplementals/Growth_light/max_g_light_slope.pdf"))
p
dev.off()
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Growth_light/max_g_light_slope.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Growth_light/max_g_light_slope.pdf')) 
)


#-------------------------------------------------------------------------
#---------------------- Fig S6b Mortality Slope with Light ------
#-------------------------------------------------------------------------
mortal <- read_csv(file.path(gdrive_path, 'data/clean_summary_tables/clean_parameters_mortality.csv'))
mortal <- mortal %>% filter(fg != "--")
# 1. Plot of maximum slope by functional group

# Remove all tree and unclassified groups
mortal$fg <- factor(mortal$fg ,levels = c("fast", "large pioneer", "medium", "small breeder", "slow"))
guild_labels2_ <- c("Fast", "Tall", "Medium", "Short", "Slow")
mort_slope <- ggplot(mortal %>% filter(parameter %in% 'slope'),
            aes(x = fg, y = mean, ymin = q025, ymax = q975)) + 
  geom_errorbar(width = 0) + geom_point(size = 4) +
  theme(axis.text.x = element_text(angle = 25,  vjust = 0.7))+
  scale_x_discrete(name = 'Life History Guild', labels = guild_labels2_) +
  scale_y_continuous(limits = c(-1.3, -0.25), breaks = c(-1.2, -0.9, -0.6, -0.3), 
                     expression(atop('Mortality Responsiveness', paste('to Light (yr'^-1,' W'^-1,' m'^-2,')')))) +
  theme_plant_small() + theme(aspect.ratio = 0.75) + theme_no_x() + annotation_custom(grob_a)
mort_slope


#combine
g_growth<- ggplotGrob(growth_slope)
g_mort <- ggplotGrob(mort_slope)

combo <- rbind(g_mort, g_growth, size = "first")
combo$widths <- unit.pmax(g_mort$widths, g_growth$widths)
grid.newpage()
grid.draw(combo)
ggsave(combo, height = 8.6, width = 6, filename = file.path(gdrive_path,'Figures/Supplementals/Max_growth_mort_light/box_growth_mort_light.pdf'))

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Max_growth_mort_light/box_growth_mort_light.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Max_growth_mort_light/box_growth_mort_light.pdf')) 
)

#-------------------------------------------------------------------------
#-------------------   Fig S7 Piecewise Scaling Slopes  -------------------

piece <- file.path(gdrive_path, 'data/data_piecewisefits')
ics <- read.csv(file.path(piece , 'piecewise_ics_by_fg.csv'), stringsAsFactors = FALSE)
ics$fg <- factor(ics$fg , labels = c("All", "Fast", "Tall", "Slow", "Short", "Medium", "Unclassified"))

slopes <- read.csv(file.path(piece , 'piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
slopes$fg <- factor(slopes$fg , labels = c("All", "Fast", "Tall", "Slow", "Short", "Medium", "Unclassified"))
slopes$variable <- factor(slopes$variable, labels = c("Density", "Individual Growth", "Production"))
colors <- c("sienna4", "yellowgreen", "springgreen4")

# Using 3 segment density and 1 segment production
p <- ggplot(slopes %>% filter((dens_model == 3 & is.na(prod_model)) |
                                (is.na(dens_model) & prod_model == 1) | 
                                (dens_model == 3 & prod_model == 1), 
                              !fg %in% 'Unclassified'), 
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
  scale_y_continuous(breaks = c(-15, -8, -5, -3, -2, 0, 2)) +
  theme_bw() + 
  theme(strip.background = element_blank(), panel.grid = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = 'bottom', 
        legend.spacing.x=unit(.2, "cm"), legend.title=element_blank(),
        legend.text=element_text(size = 12),axis.title = element_text(size = 15),
        axis.text = element_text(color = "black", size = 11)) +
  labs(y = 'Slope') +
  ggtitle('Fitted Slopes: \n 3 Segment Density & 1 Segment Growth Models')+
  theme(plot.title = element_text(hjust=0.5)) #+ theme(legend.title=element_blank())
p 
pdf(file.path(gdrive_path, "Figures/Supplementals/Piecewise_Slopes/slopes_1_seg_growth.pdf"))
p
dev.off()

# Extra: Using 3 segment density and 2 segment production
p <- ggplot(slopes %>% filter((dens_model == 3 & is.na(prod_model)) |
                                (is.na(dens_model) & prod_model == 2) | 
                                (dens_model == 3 & prod_model == 2), !fg %in% 'Unclassified'), 
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
pdf(file.path(gdrive_path, "Figures/Supplementals/Piecewise_Slopes/slopes_2_seg_growth.pdf"))
p
dev.off()

# --------------------------- Fig S8: Heat map of growth Scaling 
raw_prod <- do.call(rbind, map2(alltreedat, seq(1985,2010,5), function(x, y) cbind(year = y, x %>% filter(!recruit) %>% select(fg, dbh_corr, production))))

raw_prod <- raw_prod %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg', fg), 'unclassified'))


hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000,10000), 
                                             labels = c(1,10,100,1000,10000), limits=c(1,10000))


plot_prod_withrawdata2 <- function (year_to_plot = 1995, 
                                    fg_names = c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified"), 
                                    full_names = c("Fast", "Slow", "Pioneer", "Breeder", "Medium", "Unclassified"), 
                                    func_names = c("power law"), #, "2-segment power law"),
                                    x_limits = c(1, 300), x_breaks = c(1, 10, 100), 
                                    y_limits, 
                                    y_breaks, 
                                    x_name = "Diameter (cm)", 
                                    y_name = expression(paste("Growth (kg yr"^-1,")")),
                                    line_types = c("solid","dashed"),
                                    aspect_ratio = 0.75, 
                                    hex_scale = ggplot2::scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9,"RdYlBu")), bias = 1)(50), 
                                                                              trans = "log", name = "Individuals", 
                                                                              breaks = c(1, 10, 100, 1000, 10000), 
                                                                              labels = c(1, 10, 100, 1000, 10000), 
                                                                              limits = c(1, 10000)), 
                                    obsdat = raw_prod, 
                                    preddat = fitted_indivprod,
                                    plot_abline = TRUE, 
                                    abline_slope = 2, 
                                    abline_intercept = -2.1, 
                                    plot_fits = FALSE) 
{
  obsdat <- obsdat %>% dplyr::filter(fg %in% fg_names, year ==  year_to_plot) %>% 
    arrange(desc(fg))
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% 
    dplyr::summarize(min_obs = min(production), 
                     max_obs = max(production)) %>% 
    arrange(desc(fg))
  preddat <- preddat %>% filter(prod_model == 1) %>%
    dplyr::left_join(obs_limits) %>% 
    dplyr::filter(fg %in% fg_names, year == year_to_plot) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% 
    dplyr::filter(dbh >=  min_obs & dbh <= max_obs) %>% 
    dplyr::mutate(prod_model = factor(prod_model, 
                                      labels = func_names)) %>%
    arrange(desc(fg))
  labels <- setNames(full_names, fg_names)
  p <- ggplot2::ggplot() + 
    ggplot2::geom_hex(data = obsdat, ggplot2::aes(x = dbh_corr, y = production)) + 
    ggplot2::facet_wrap(~fg, ncol = 2, labeller = ggplot2::labeller(fg = labels)) + 
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks) + 
    ggplot2::scale_linetype_manual(name = "Growth fit", values = line_types) + 
    hex_scale + 
    theme_plant_small() + 
    ggplot2::coord_fixed(ratio = aspect_ratio) + 
    ggplot2::theme(legend.position = "right", 
                   strip.background = ggplot2::element_blank(), 
                   strip.text = ggplot2::element_text(size = 13), 
                   legend.key = ggplot2::element_rect(fill = NA))
  if (plot_abline) 
    p <- p + ggplot2::geom_abline(slope = abline_slope, intercept = abline_intercept, linetype = "dashed") + 
    ggplot2::guides(linetype = "none")
  if (plot_fits) 
    p <- p + ggplot2::geom_line(data = preddat, 
                                ggplot2::aes(x = dbh, y = q50), size = 0.5) #, group = prod_model, linetype = prod_model
  return(p)
}

p <- plot_prod_withrawdata2(year_to_plot = 1995,
                            fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
                            full_names = c('Fast', 'Tall', 'Slow', 'Short', 'Medium', 'Unclassified'),
                            x_limits = c(1, 316),
                            x_breaks = c(1,10, 100),
                            y_limits = c(0.001, 1000),
                            y_breaks = c(.001, .1, 10, 1000),
                            line_types = c('solid', 'dashed'),
                            hex_scale = hex_scale_log_colors,
                            plot_abline = FALSE,
                            plot_fits = TRUE)

p <- p + theme(legend.position = 'right', legend.text = element_text(size = 13), 
               legend.title = element_text(size = 15)) + 
  scale_y_log10(labels = c(0.01, 1, 100), 
                breaks = c(0.01, 1, 100), 
                name = expression(paste('Growth (kg yr'^-1,')'))) 

p 
pdf(file.path(gdrive_path, 'Figures/Supplementals/Growth_Hex/growth_hex.pdf'))
p
dev.off()


#-------------------------------------------------------------------------------
# ------------------------   Fig S9: Growth WAIC of Piecewise Models  -----------
#-------------------------------------------------------------------------------

# This section was edited by QDR, 20 Jun 2019, for the updated model fits.

# Growth model
base_size <- 11
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
pdf(file.path(gdrive_path, 'Figures/Supplementals/WAIC/WAIC_growth.pdf'))
grid.draw(p1)
dev.off()

#-------------------------------------------------------------------------------
# ------------------------   Fig S10: Growth WAIC of Piecewise Models  -----------------------------------
#-------------------------------------------------------------------------------

# Density model

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
pdf(file.path(gdrive_path, 'Figures/Supplementals/WAIC/WAIC_density.pdf'))
grid.draw(p1)
dev.off()


#---------------------------------------------------------------------------------------------
#----------------------- Fig. S11: total production from 1985 - 2010 --------------------------
#---------------------------------------------------------------------------------------------
minmax_prod_bycensus <- obs_totalprod %>%
  filter(bin_value > 0, !fg %in% 'unclassified') %>%
  group_by(fg, bin_midpoint) %>%
  summarize(range_min = min(bin_value), range_max = max(bin_value))
minmax_prod_bycensus $fg <- factor(minmax_prod_bycensus $fg , labels = c("All", "Fast", "Tall", "Slow", "Short", "Medium"))

p <- ggplot(minmax_prod_bycensus, 
            aes(x = bin_midpoint, ymin = range_min, ymax = range_max, color = fg)) +
  geom_errorbar(size = 1) +
  scale_x_log10(name = 'Diameter (cm)', breaks = c(1,3,10,30,100,300)) + 
  scale_y_log10(expression(paste('Production (kg ha'^-1,' cm'^-1,' yr'^-1,')')),
                breaks = 10^(-2:3), labels = as.character(10^(-2:3)), limits = c(0.1, 200)) +
  scale_color_manual(values = guild_colors2, labels = guild_labels2) +
  theme_plant_small(legend = T) + guide2
p

p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p1)

pdf(file.path(gdrive_path,'Figures/Supplementals/Total_Prod_Range/Total_Prod_Range.pdf'))
grid.draw(p)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Total_Prod_Range/Total_Prod_Range.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Total_Prod_Range/Total_Prod_Range.pdf')) 
)

#-------------------------------------------------------------------------------------------
#----------------------- Fig S12: individual light capture scaling --------------------------
#-------------------------------------------------------------------------------------------

alpha_value <- 1
hexfill2 <- scale_fill_gradient(low = 'forestgreen', high = 'navy', trans = 'log', breaks = c(1,3,10,30,100,300))
exl <- expression(atop('Intercepted Light', paste('per Individual (W)')))
exd <- 'Diameter (cm)'

labels = trans_format("log10", math_format(10^.x))

p <- ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received)) +
  geom_pointrange(data = unscaledlightbydbhcloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = mean, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, breaks = c(1,100,10000, 1000000), limits = c(1,1000000), 
                labels = trans_format("log10", math_format(10^.x))) +
  theme_plant() + theme(legend.position = "right", legend.text = element_text(size = 15), legend.title = element_text(size = 16))+
  hex_scale_log_colors +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))# +
p
pdf(file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/indiv_light.pdf'))
p
dev.off()


system2(command = "pdfcrop", 
        args  = c(file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/indiv_light.pdf'), 
                  file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/indiv_light.pdf')) 
)



# ------------------------------- Fig S13: Total Volume Vollume Scaling ------------------------------------

# Section added by QDR 20 June 2019: new plots of total light scalings and total volume scalings, including fits and CIs

# Fitted values for individual light, total light, and total volume

fitted_totalvol <- read.csv(file.path(fp, 'fitted_totalvol.csv'), stringsAsFactors = FALSE)
totalvolbins_fg <- read.csv(file.path(fp, 'obs_totalvol.csv'), stringsAsFactors = FALSE)


#rough_slope_all <- (log(1650.623080)  -	log(844.331539)) / (log(36.878615) - log(4.919149))	

totalvolbins_fg <- totalvolbins_fg %>%
  filter(bin_count >= 20)
p <- plot_totalprod2(year_to_plot = 1995,
                     fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                     model_fit_density = DENS, 
                     model_fit_production = PROD,
                     x_limits = c(0.9, 150),
                     x_breaks = c(1, 10, 100),
                     y_limits = c(5, 5500),
                     y_breaks = c(1, 10, 100, 1000),
                     y_labels = c(1, 10, 100, 1000),
                     y_name = expression(paste('Total Crown Volume (m'^3, ' cm'^-1, ' ha'^-1,')')), 
                     preddat = fitted_totalvol,
                     obsdat = totalvolbins_fg, 
                     plot_abline = FALSE,
                     geom_size = 4)
p
p0 <- p + scale_y_continuous(position = "left", trans = "log10", limits = c(9, 5000),
                             name = expression(atop('Total Crown Volume',paste('(m'^3, ' cm'^-1,' ha'^-1,')')))) +
  theme_plant_small(legend = TRUE)  + guide +
  scale_fill_manual(values = guild_fills2, labels = guild_labels2) 
plot(p0)

p_tot_vol<- set_panel_size(p0, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p_tot_vol)
pdf(file.path(gdrive_path,'Figures/Supplementals/Total_Vol/total_vol.pdf'))
grid.draw(p_tot_vol)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Total_Vol/total_vol.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Total_Vol/total_vol.pdf')) 
)

## -------------------- ------Fig S14: PCA score   -----------------------------

grob_p <- grobTree(textGrob("Short-Tall Continuum", x = 0.13, y = 0.93,  hjust = 0,
                            gp = gpar(col = "black", fontface = "italic", fontface = "bold", fontsize = 15))) 

grob_f <- grobTree(textGrob("Slow-Fast Continuum", x = 0.13, y = 0.83,  hjust = 0,
                            gp = gpar(col = "gray43",  fontface = "italic",fontsize =15))) 


fastslowscore_bin_bylight <- fastslowscore_bin_bylight_byyear %>%
  left_join(fastslow_stats_bylight_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Fast-Slow")

breederscore_bin_bylight <- breederscore_bin_bylight_byyear  %>%
  left_join(breeder_stats_bylight_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Short-Tall")

PCA_score_by_light <- rbind(fastslowscore_bin_bylight,breederscore_bin_bylight)

error_bar_width <- .15

PCA_light <- PCA_score_by_light %>%
  filter(year == 1995, n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, fill = ID)) +
  geom_errorbar(width = error_bar_width) +theme_plant() +
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey"))+
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+ annotation_custom(grob_a) +
  annotation_custom(grob_p) + annotation_custom(grob_f) +
  scale_x_log10(limits=c(1.8,450), breaks = c(3, 30, 300),
                name = expression(paste('Light per Crown Area (W m'^-2,')'))) + 
  scale_y_continuous(limits=c(-1.5,1.25),breaks=c(-1,0,1),name = 'PCA Score') + theme_plant_small()
PCA_light


#---Fig S14B:  PCA by Diameter

fastslowscore_bin_bydiam <- fastslowscore_bin_bydiam_byyear %>%
  left_join(fastslow_stats_bydiam_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Fast-Slow")

breederscore_bin_bydiam <- breederscore_bin_bydiam_byyear %>%
  left_join(breeder_stats_bydiam_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = 'Short-Tall')

score_bin_bydiam <- rbind(fastslowscore_bin_bydiam, breederscore_bin_bydiam)


grob_a <- grobTree(textGrob("a", x = 0.04, y = 0.93,  hjust = 0,
                            gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob_b <- grobTree(textGrob("b", x = 0.04, y = 0.93,  hjust = 0,
                            gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 

PCA_diam <- score_bin_bydiam %>%
  filter(year == 1995, n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  geom_errorbar(width = error_bar_width) +theme_plant_small() +
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey"))+
  scale_x_log10(name = 'Diameter (cm)', limits = c(.8,100)) + 
  annotation_custom(grob_b) +
  scale_y_continuous(limits=c(-1.5,1.25),breaks=c(-1,0,1),name = "PCA Score") #+

PCA_diam

# -- combine
g_light  <- ggplotGrob(PCA_light)
g_diam <- ggplotGrob(PCA_diam)
gPCA <- rbind(g_light, g_diam, size = "first")
gPCA$widths <- unit.pmax(g_light$widths,g_diam $widths)
grid.newpage()
grid.draw(gPCA)

ggsave(gPCA, height = 8, width = 11, filename = file.path(gdrive_path,'Figures/Supplementals/Ratios_PCA/PCA.pdf'))

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Ratios_PCA/PCA.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Ratios_PCA/PCA.pdf')) 
)
#-------------------------------------------------------------------------------
#----------------------Fig S15: Light Scaling -------------------
#----------------------------------------------------------------------------

# Load plotting data
load(file.path(gdrive_path, 'data/data_forplotting/light_scaling_plotting_data.RData'))
fg_names <- c("Fast", "Tall", "Slow", "Short")
fill_scale <- scale_fill_manual(values = guild_fills[1:4], name = NULL, labels = fg_names, guide = guide_legend(override.aes = list(shape = 21)))
color_scale <- scale_color_manual(values = guild_colors[1:4], name = NULL, labels = fg_names, guide = FALSE)

area_core <- 42.84

obs_indivprod_lightarea$fg <- factor(obs_indivprod_lightarea$fg , labels = c("Fast", "Tall", "Slow", "Short"))


#------------ Fig S15a: Growth with light

l_growth <- ggplot() +
  geom_ribbon(data = prod_pred_dat_lightarea %>% filter(light_area >= 7), 
              aes(x = light_area, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.3, show.legend = F) +
  geom_line(data = prod_pred_dat_lightarea %>% filter(light_area >= 7),  
            aes(x = light_area, y = q50, group = fg, color = fg)) +
  geom_point(data = obs_indivprod_lightarea %>% filter(mean_n_individuals >= 20), 
             aes(x = bin_midpoint, y = mean, group = fg, fill = fg), 
             shape = 21, color = 'black', size = 4, show.legend = FALSE) +
  scale_x_log10(breaks = c(3, 30, 300), name = parse(text = 'Light~per~crown~area~(W~m^-2)'), limits=c(1, 600)) +
  scale_y_log10(labels = signif, limits = c(0.01, 100), name = parse(text = 'Growth~(kg~yr^-1)')) +
  annotation_custom(grob_fast) + annotation_custom(grob_tall) + 
  annotation_custom(grob_slow3) + annotation_custom(grob_short3) +
  theme_plant_small() + theme_no_x() +
  fill_scale +
  color_scale
l_growth

l_growth<- set_panel_size(l_growth , width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(l_growth)

pdf(file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/light_growth_v2.pdf'))
grid.draw(l_growth_low2)
dev.off()
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/light_growth_v2.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/light_growth_v2.pdf')) 
)


#------------ Fig S15b: GAbundance with light

l_abun<- ggplot() +
  geom_ribbon(data = dens_pred_dat_lightarea %>% filter(light_area >= 7), 
              aes(x = light_area, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.3) +
  geom_line(data = dens_pred_dat_lightarea %>% filter(light_area >= 7), 
            aes(x = light_area, y = q50, group = fg, color = fg)) +
  geom_point(data = obs_dens_lightarea %>% filter(bin_count >= 20, bin_value > 0), 
             aes(x = bin_midpoint, y = bin_value/area_core, group = fg, fill = fg),
             shape = 21, color = 'black', size = 4, show.legend = FALSE) +
  scale_x_log10(name = parse(text = 'Light~per~crown~area~(W~m^-2)'), breaks = c(3, 30, 300),limits = c(1,400)) +
  scale_y_log10(labels = signif, limits = c(0.01, 300), name = parse(text = 'Abundance~(ha^-1~cm^-1)')) +
  theme_plant_small() +
  fill_scale +
  color_scale
l_abun

l_abun2 <- set_panel_size(l_abun, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(l_abun2 )

pdf(file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/light_density_v2.pdf'))
grid.draw(l_abun2)
dev.off()
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/light_density_v2.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/light_density_v2.pdf')) 
)


#----Fig S15C:Total prod scaling with light-----------

l_prod<- ggplot() +
  geom_ribbon(data = totalprod_pred_dat_lightarea %>% filter(light_area >= 7), 
              aes(x = light_area, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.3, show.legend = F) +
  geom_line(data = totalprod_pred_dat_lightarea %>% filter(light_area >= 7), 
            aes(x = light_area, y = q50, group = fg, color = fg)) +
  geom_point(data = obs_totalprod_lightarea %>% filter(bin_count >= 20, bin_value > 0), 
             aes(x = bin_midpoint, y = bin_value/area_core, group = fg, fill = fg),
             shape = 21, color = 'black', size = 4, show.legend = FALSE) +
  scale_x_log10(breaks = c(3, 30, 300), name = parse(text = 'Light~per~crown~area~(W~m^-2)'), limits = c(1,600)) +
  scale_y_log10(labels = signif, limits = c(0.01, 10), position = "right",
                name  = expression(atop('Production', paste('(kg yr'^-1,' cm'^-1,' ha'^-1,')')))) +
  theme_plant_small() +theme(aspect.ratio = 1) +# theme_no_x() + 
  fill_scale + 
  color_scale

l_prod
l_prod2 <- set_panel_size(l_prod, width=unit(10.25,"cm"), height=unit(7,"cm"))
l_prod2 <- set_panel_size(l_prod, width=unit(14,"cm"), height=unit(14,"cm"))

grid.newpage()
grid.draw(l_prod2)

pdf(file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/light_totalprod_v4.pdf'))
grid.draw(l_prod_low2)
dev.off()
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/light_totalprod_v3.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Light_Scaling/light_totalprod_v3.pdf')) 
)


#---- Fig S15d: production ratio with light

dens_ratio <- prod_ratio_diam %>% 
  filter(density_ratio > 0) %>%
  filter(n_individuals >= 20) %>%
  ggplot() +
  geom_ribbon(aes(x = dbh, ymin = q025, ymax = q975, group = ratio, fill = ratio), 
              alpha = 0.4, data = ratio_fitted_diam_density) +
  geom_line(aes(x = dbh, y = q50, group = ratio, color = ratio), 
            data = ratio_fitted_diam_density) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_point(aes(x = bin_midpoint, y = density_ratio, fill = ID),
             shape = 21, size = 4,  stroke = .5,  color = "black") +
  scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey")) +
  scale_color_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey50")) +
  scale_x_log10(limits = c(1,100), breaks = c(1,10, 100), 
                name = expression(paste('Diameter (cm)'))) + 
  scale_y_log10(labels = signif, breaks = c(0.01,0.1, 1,10,100,1000), 
                limits=c(0.01,200),
                name = expression("Density Ratio"), position = "right") + 
  theme_plant() #or small
dens_ratio 


# -------------------- ----- Supp: Production ratio by diameter, abundance by light  -----------------------------

p_r_d <- prod_ratio_diam   %>%
  filter(n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  geom_point(shape = 21, size = 4,  stroke = .5, color = "black")+
  scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey")) +
  theme_plant() + annotation_custom(grob_b) + 
  scale_x_log10(name = 'Diameter (cm)', limits=c(1,150), breaks=c(1, 10, 100)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.03,100), 
                name = expression("Production Ratio")) 

p_r_d 

a_r_l <- prod_ratio_light   %>%
  filter(n_individuals >= 20) %>%
  ggplot() + 
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_point(aes(x = bin_midpoint, y = density_ratio, fill = ID), 
             shape = 21, size = 4.5,  stroke = .5, color = "black") +
  scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey"))+
  theme_plant() +  annotation_custom(grob_a) + annotation_custom(grob_f) + annotation_custom(grob_p) +
  scale_x_log10(name = expression(paste('Light per Crown Area (W m'^-2,')')), limits=c(2,330), breaks=c(3,  30,  300)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100), limits=c(0.01,100),
                name = expression("Density Ratio")) +
  theme(plot.title = element_text(size = 14, color = "black")) 

a_r_l 

g_p_r_d  <- ggplotGrob(p_r_d )
g_a_r_l <- ggplotGrob(a_r_l )
g_ratio_sup <- rbind(g_a_r_l ,g_p_r_d , size = "first")
g_ratio_sup$widths <- unit.pmax(g_a_r_l$widths, g_p_r_d $widths)
grid.newpage()
grid.draw(g_ratio_sup )

ggsave(g_ratio_sup , height = 8, width = 11, filename = file.path(gdrive_path,'Figures/Supplementals/Ratios_PCA/Ratio.pdf'))
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/Ratios_PCA/Ratio.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/Ratios_PCA/Ratio.pdf')) 
)



#---------------Into the weeds of the light scaling -------------
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
fg_names <- c('Fast', 'Tall', 'Slow', 'Short')

dens_pred_fg_lightarea_quantiles <- map2(dens_pred_fg_lightarea, map(stan_data_list, 'N'), function(dat, N) {
  do.call(cbind, dat) %>%
    sweep(2, N/area_core, `*`) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    t %>% as.data.frame %>% setNames(c('q025', 'q50', 'q975')) %>%
    mutate(light_area = la_pred)
})

dens_pred_dat_lightarea <- map2_dfr(dens_pred_fg_lightarea_quantiles, fg_names,
                                    ~ data.frame(fg = .y, .x)) %>%
  mutate(fg = factor(fg, levels = fg_names))

prod_pred_fg_lightarea_quantiles <- map(prod_pred_fg_lightarea, function(dat) {
  do.call(cbind, dat) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    t %>% as.data.frame %>% setNames(c('q025', 'q50', 'q975')) %>%
    mutate(light_area = la_pred)
}) 

prod_pred_dat_lightarea <- map2_dfr(prod_pred_fg_lightarea_quantiles, fg_names,
                                    ~ data.frame(fg = .y, .x)) %>%
  mutate(fg = factor(fg, levels = fg_names))

totalprod_pred_fg_lightarea_quantiles <- map2(totalprod_pred_fg_lightarea, map(stan_data_list, 'N'), function(dat, N) {
  dat %>%
    sweep(2, N/area_core, `*`) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    t %>% as.data.frame %>% setNames(c('q025', 'q50', 'q975')) %>%
    mutate(light_area = la_pred)
})

totalprod_pred_dat_lightarea <- map2_dfr(totalprod_pred_fg_lightarea_quantiles, fg_names,
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

obs_dens_lightarea <- data_to_bin %>%
  group_by(fg) %>%
  group_modify(~ logbin_setedges(x = .$light_received_byarea, edges = binedgedat))

obs_totalprod_lightarea <- data_to_bin %>%
  group_by(fg) %>%
  group_modify(~ logbin_setedges(x = .$light_received_byarea, y = .$production, edges = binedgedat))

obs_indivprod_lightarea <- data_to_bin %>%
  group_by(fg) %>%
  group_modify(~ cloudbin_across_years(dat_values = .$production, dat_classes = .$light_received_byarea, edges = binedgedat, n_census = 1))


# Themes

fill_scale <- scale_fill_manual(values = guild_fills[1:4], name = NULL, labels = fg_names, guide = guide_legend(override.aes = list(shape = 21)))
color_scale <- scale_color_manual(values = guild_colors[1:4], name = NULL, labels = fg_names, guide = FALSE)

