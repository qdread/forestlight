# "Extra" plotting code
# Deprecated plots and those not included in either the main MS or supplemental figures

# Set path to data on google drive
gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users',Sys.info()['user'],'Google Drive/ForestLight'))

library(forestscaling) # Packaged all the functions and ggplot2 themes here!

library(tidyverse)
library(egg)
library(scales)
library(RColorBrewer)
library(gtable)
library(grid)
library(reshape2)

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



# Unused variants of hexagon plots ----------------------------------------


# ------------------- Each group
ggplot() +
  geom_point(alpha = 0.05, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownvolume), color = 'chartreuse3') +
  geom_pointrange(data = lightpervolfakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg, ncol = 2, labeller = labeller(fg = labels)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv) +
  theme_facet2()
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
  theme_facet2() +
  hex_scale_log_colors  +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# All together
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownarea)) +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, breaks = c(1,10,100,1000), limits = c(1,1000)) +
  theme_plant() +
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
  theme_facet2() +
  hex_scale_log_colors  +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# All together
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownvolume)) +
  geom_pointrange(data = lightpervolfakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv, breaks = c(1,10,100,1000), limits = c(1,1000)) +
  theme_plant() +
  hex_scale_log_colors  +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))



# More unused hexagon plots -----------------------------------------------

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
  theme_plant() +
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
  geom_ribbon() + geom_line() + theme_facet2()

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



#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(gray.colors(9, start=.9, end=.1), bias=3)(50),
#                                  trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))

#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(c('khaki1', 'gold', 'red3'), bias=3)(50),
#                                  trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))

#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(c('aliceblue', 'forestgreen','khaki1','red3'), bias=3)(50),
#                        trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))
#hex_scale_log <- scale_fill_gradientn(colours = colorRampPalette(c('forestgreen', 'red3'), bias=3)(50),
#                        trans = 'log', name = 'Number of\nindividuals', breaks = c(1,10,100,1000), labels = c(1,10,100,1000))
#hex_scale_3b <- scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(9,'RdYlBu'), bias=0.1)(50)), guide = F)

# Another unused variant of light by area plots --------------------------------------------------

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
  theme_plant() + theme(aspect.ratio = 0.7) #+
# theme(panel.border = element_rect(fill=NA),
#      legend.position = c(0.2, 0.8))


p_median_1panel
p<-set_panel_size(p_median_1panel, width=unit(11.8,"cm"), height=unit(8.26,"cm"))
plot(p)


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
  theme_plant() +
  ggtitle('Symmetry of total growth and total incoming energy patterns', 'Slope of scaling relationship at 10 cm dbh')

###
# Observed values with fitted lines superimposed (update to Figure 4)


# light per volume

# Light Vs Size Plot
# Edited 21 March: include volume in addition to area.

#load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

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
pdf(file.path(gdrive_path, 'Figures/Growth_light/light_crown_area_fg.pdf'))
p
dev.off()
lightperareafakebin_fg$fg
# All together
p <- ggplot() +
  geom_point(alpha = 0.01, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownarea), color = 'chartreuse3') +
  geom_pointrange(data = lightperareafakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, limits = c(1, 1000), breaks = c(1, 10, 100, 1000)) +
  theme_plant() 

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
  theme_facet() + 
  theme(axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 13)) + 
  theme(legend.position="none") +
  theme(strip.background = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = .75),
        strip.text.x = element_text(size = 12, face = "italic"),
        panel.grid = element_blank(), legend.position = 'none') +
  scale_color_manual(values = guild_fills_nb0)
p

pdf(file.path(gdrive_path, 'Figures/Growth_light/light_crown_vol_fg.pdf'))
p
dev.off()

colors <- c("sienna4", "yellowgreen", "springgreen4")




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
  theme_plant() + #geom_errorbar(aes(width = width), position=p_dodge) +
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
  theme_plant() + geom_errorbar(aes(width = width), position = p_dodge) +
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
  theme_plant() + geom_errorbar(aes(width = width), position=p_dodge) +
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



# Light by diameter model -------------------------------------------------

# Faceted individual light plot
# Each group
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received)) +
  geom_pointrange(data = unscaledlightbydbhfakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg, ncol = 2, labeller = labeller(fg = fg_labels)) +theme_plant+
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, breaks = c(1,100,10000, 1000000), limits = c(1,1000000), 
                labels = trans_format("log10", math_format(10^.x))) +
  theme_facet2() +
  hex_scale_log_colors +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))


# Very simple model
loglogregressions <- alltree_light_95 %>%
  group_by(fg) %>%
  do(model = lm(log10(light_received) ~ log10(dbh_corr), data = .))
library(broom)
loglogregressions$model
lapply(loglogregressions$model, summary)
#(Intercept)     0.903027   0.005262   171.6   <2e-16 ***
#log10(dbh_corr) 2.343608   0.007530   311.2   <2e-16 ***


# Distinguish 2 and 5 census ratio plots ----------------------------------

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
  geom_errorbar(width = error_bar_width) + theme_plant() +
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Breeder" = "black", "Fast-Slow" = "grey"))+
  
  scale_x_log10(limits=c(.7,160),breaks=c(1,10, 100), name = expression(paste('Diameter (cm)'))) + 
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.003,100),
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

p <- fastslow_stats_breeder_bydiam_2census %>% 
  filter(density_ratio_mean > 0) %>%
  filter(mean_n_individuals > 10) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, 
             ymin = density_ratio_min, ymax = density_ratio_max, fill = ID)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  geom_errorbar(width = error_bar_width) + theme_plant() +
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black")+
  scale_fill_manual(values = c("Breeder-Pioneer" = "black", "Fast-Slow" = "grey"))+
  
  scale_x_log10(limits=c(.7,330),breaks=c(1,10, 100), name = expression(paste('Diameter (cm)'))) + 
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.006,100),
                name = expression("Ratio")) 


p
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
plot(p1)
pdf(file.path(gdrive_path, 'Figures/Fig_6/Fig_6b_Density_2yr.pdf'))
plot(p1)
dev.off()
