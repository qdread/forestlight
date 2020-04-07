
DENS = 3
PROD = 1


#maybe change theme_plant to include the legend?

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users/jgradym/Google Drive/ForestLight'))

#library(broom)
library(forestscaling) 
library(tidyverse)

#library(egg)
library(scales)
library(RColorBrewer)
library(grid)
#library(gtable)
#library(grid)
library(reshape2)
#library(hexbin)

#
guild_fills <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")
guild_fills2 <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")
guild_colors <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")
guild_colors2 <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray")

fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('Fast','Tall', 'Slow', 'Short', 'Medium')
fg_labels2 <- c('Fast','Tall', 'Slow', 'Short', 'Medium', 'All')
fg_labels3 <- c('Fast','Tall', 'Slow', 'Short', 'Medium', 'All', 'Unclassified')

fg_lookup <- data.frame(fg = c('fg1','fg2','fg3','fg4','fg5','all','unclassified'), 
                        fg_name = c('Fast','Tall','Slow','Short','Medium','All','Unclassified'))

year_to_plot = 1995
geom_size <- 4
#probaly worth adding the legend back in
theme_plant <- theme_plant() + 
  theme(legend.position = "right", legend.text = element_text(size = 14 ), 
                      legend.key = element_blank())
guide <- guides(fill=guide_legend(title=NULL), color = F, override.aes=list(fill=NA))

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
################################################################################################
# ------------------------------ Fig 1: Light Interception ---------------------------------
################################################################################################


lightperareacloudbin_fg <- read.csv(file.path(fp, 'lightperareacloudbin_fg.csv'), stringsAsFactors = FALSE)
lightpervolcloudbin_fg <- read.csv(file.path(fp, 'lightpervolcloudbin_fg.csv'), stringsAsFactors = FALSE)
unscaledlightbydbhcloudbin_fg <- read.csv(file.path(fp, 'unscaledlightbydbhcloudbin_fg.csv'), stringsAsFactors = FALSE)


exl <- expression(atop('Light per Crown Area', paste('(W m'^-2, ')')))
exv <- expression(atop('Light per Crown Volume', paste('(W m'^-3, ')')))
exd <- 'Diameter (cm)'


#----------------------   Fig 1a: Light per crown volume by diameter -----------------------------

Fig_1a <- ggplot() + geom_point(alpha = 0.01, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownarea), color = 'chartreuse3') +
  geom_pointrange(data = lightperareacloudbin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_plant 
Fig_1a



#----------------------   Fig 1b: Light per crown volume by diameter -----------------------------



Fig_1b <- ggplot() +
  geom_point(alpha = 0.01, data = alltree_light_95, 
             aes(x = dbh_corr, y = light_received/crownvolume), color = 'chartreuse3') +
  geom_pointrange(data = lightpervolcloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exv, breaks = c(1, 10, 100)) +
  theme_plant
Fig_1b 
################################################################################################
# ------------------------------ Fig 3: Life Histories ---------------------------------
################################################################################################

# ### Fig 3a, PCA

fgbci <- read.table(file.path(gdrive_path, 'data/Ruger/fgroups_dynamics_new.txt'), stringsAsFactors = FALSE)
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5)) # correct in advance?
fgbci$PC_slow_to_fast <- -fgbci$X1new
fgbci$PC_breeder_to_pioneer <- fgbci$X2new


Fig_3a <- ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, fill = factor(fg5))) +
  geom_point(shape = 21, size = geom_size, color = "black") + 
  labs(x = 'Slow - Fast Tradeoff', y = 'Recruitment - Stature Tradeoff') +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,3))+
  scale_color_manual(values = guild_colors)+
  scale_fill_manual(values = guild_fills, labels = fg_labels) + 
  theme_plant + guide
Fig_3a

#########  Fig 3b
# Axis titles
title_x <- expression(paste('Light per Crown Area (W m'^-2,')',sep=''))
title_y <- expression(atop('Growth per Crown Area', paste('(kg yr'^-1, ' m'^-2,')', sep='')))
scale_y_log10(name =  expression(atop('Growth per Crown Area',
                                      paste('(kg y'^-1, ' m'^-2,')'))))


obs_light_binned <- read.csv(file.path(fp, 'obs_light_binned.csv'), stringsAsFactors = FALSE)
obs_light_raw <- read.csv(file.path(fp, 'obs_light_raw.csv'), stringsAsFactors = FALSE)
pred_light <- read.csv(file.path(fp, 'pred_light.csv'), stringsAsFactors = FALSE)
param_ci <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/lightbyarea_paramci_by_fg.csv'), stringsAsFactors = FALSE)

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

Fig_3b <- ggplot(obs_light_binned_plotdata) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
              aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.4) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot),
            aes(x = light_area, y = q50, color = fg)) +
  geom_errorbar(aes(x = bin_midpoint, ymin = q25, ymax = q75, 
                    group = fg, color = fg, width = 0),
                position = position_dodge(width = dodge_width)) + 
  geom_point(aes(x = bin_midpoint, y = mean, group = fg, fill = fg),
             size = 4, shape = 21, position = position_dodge(width = dodge_width)) +
  
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, position = "right", breaks = c(0.01, 0.03, 0.1, 0.3), 
                labels = c( 0.01, 0.03, 0.1, 0.3)) +
  scale_color_manual(values = guild_colors) +
  scale_fill_manual(values = guild_fills, labels = fg_labels) +
  theme_plant +guide
Fig_3b

### ################## Fig 3c, Mortality

bin_mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_mortalitybins.csv')) # Load binned data
mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting/obs_mortalityindividuals.csv')) # Load raw data
fitted_mort <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/mortality_ci_by_fg.csv')) # Load fitted values

# Truncate fitted mortality lines to not exceed the range of the observed data, using 20 individuals as the cutoff.
obs_range_mort <- bin_mort %>% 
  filter(variable %in% 'light_per_area', lived + died >= 20) %>%
  group_by(fg) %>%
  summarise(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))

fitted_mort_trunc <- fitted_mort %>%
  left_join(obs_range_mort) %>%
  filter(light_per_area >= min_obs & light_per_area <= max_obs)


#Mortality
Fig_3c<- ggplot(data = fitted_mort_trunc %>% mutate(fg = factor(fg, labels = fg_labels))) +
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
  scale_color_manual(values = guild_colors) +
  scale_fill_manual(values = guild_fills) +
  theme_plant + guide
Fig_3c




########################################################################################
# ------------------------------- Fig 4 Scaling Plots ------------------------------------
########################################################################################



# I'm trying to keep the original function, but we should reverse geom order here: %>% arrange(desc(fg))
# maybe change names so that the names and order is All, Medium, Short, Tall, Slow, Fast 
obs_indivprod <- obs_indivprod %>%
  filter(mean_n_individuals >= 20)
Fig_4a <- plot_prod(year_to_plot = 1995,
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
Fig_4a+ theme_plant +guide
Fig_4a 

obs_dens <- obs_dens %>%
  filter(bin_count >= 20)
Fig_4b <- plot_dens(year_to_plot = 1995,
                fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                model_fit = DENS,
                x_limits = c(.8, 230),
                y_limits = c(0.01, 20000),
                x_breaks = c(1, 10, 100),
                y_labels = c(0.001, 0.1, 10,1000),
                y_breaks = c(0.001, 0.1,  10, 1000))
Fig_4b +theme_plant +guide
#p <- p +annotation_custom(grob_text_b)

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
Fig_5a <- plot_totalprod(year_to_plot = 1995,
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
Fig_5a + theme_plant + guide

# Fig 5b

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

fg_full_names <- c('Fast', 'Tall', 'Slow', 'Short', 'Medium', 'All Trees', 'Unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

allslopes <- rbind(growth_slopes_atmiddle, light_slopes_atmiddle) %>%
  ungroup %>%
  mutate(fg = factor(fg, levels = fgs, labels = fg_full_names))
grob0 <- grobTree(textGrob("b", x = 0.04, y = 0.9,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob1 <- grobTree(textGrob("Light Capture", x = 0.75, y = 0.94, hjust = 0,
                           gp = gpar(col = "gold3", fontsize = 18))) 
grob2 <- grobTree(textGrob("Production", x = 0.75, y = 0.86, hjust = 0,
                           gp = gpar(col = "darkgreen", fontsize = 18)))
grob3 <- grobTree(textGrob("Energy Equivalence", x = 0.35, y = 0.51, hjust = 0,
                           gp = gpar(col = "black", fontsize = 18))) 
# Plot

Fig_5b <- ggplot(allslopes %>% filter(!fg %in% 'Unclassified'), aes(x = fg, y = q50, ymin = q025, ymax = q975, fill =  variable, color =variable)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = .75) +
  geom_point(position = position_dodge(width = 0.6), shape = 21, size = 4, color = "black", stroke = 0.5) +
  geom_errorbar(position = position_dodge(width = 0.6), size = 0.75, width = 0) +
  labs( x = NULL, y = 'Scaling Slope') +
  scale_y_continuous(limits = c(-1.05, 1.3)) +
  scale_fill_manual(values = c('gold1', 'darkolivegreen3')) +
  scale_color_manual(values = c('gold3', 'darkgreen')) +
  theme_plant() + theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "italic", size = 18)) +
  annotation_custom(grob1) + annotation_custom(grob2) + annotation_custom(grob3) +annotation_custom(grob0)
Fig_5b



########################################################################################
# ------------------------------- Fig 6, Ratio Scaling ------------------------------------
########################################################################################
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
  mutate(ID = 'Short-Tall') %>%
  rename(production_ratio = breeder_production_ratio, density_ratio = breeder_density_ratio) %>%
  rbind(fastslow_stats_bylight_byyear %>% 
          mutate(ID = 'Fast-Slow') %>%
          rename(production_ratio = fastslow_production_ratio, density_ratio = fastslow_density_ratio)) %>%
  filter(year == 1995)

# ---------------------------------   Fig 6a  Production by light ------------------------------------


prod_ratio <- prod_ratio_light   %>%
  filter(n_individuals >= 20) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio, fill = ID)) +
  geom_point(shape = 21, size = 4.5,  stroke = .5, color = "black")+
  scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey"))+
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  theme_plant + 
  scale_x_log10(name = expression(paste('Light per Crown Area (W m'^-2,')')), limits=c(2,330), breaks=c(3,  30,  300)) +
  scale_y_log10(labels=signif,breaks = c(0.01,0.1, 1,10,100,1000), limits=c(0.01,200),
                name = "Ratio") 
prod_ratio

#----------------------------- Fig 6B Abundance by Diameter ----------------------------
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
  theme_plant 
dens_ratio 



########################################################################################
# ------------------------------- Supplementals ------------------------------------
########################################################################################


# ------------------------ Individual diameter ~ diameter  -------------------------

#note: needs x axis values and ticks
obs_indivdiamgrowth <- obs_indivdiamgrowth %>%
  filter(mean_n_individuals >= 20)
plot_prod(year_to_plot = 1995,
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

# ------------------------ Growth per Crown Area ~ Light per Crown Area -------------------------

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
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), 
            aes(x = light_area, y = q50), color = 'black') +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y, labels = signif) +
  hex_scale_log_colors +
  theme_plant_small() +
  theme(strip.text = element_text(size=14))+
  guides(color = FALSE) +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA),
        legend.position = c(0.7, 0.15),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=15))
p_hex_panels

#---------------------------- Median binned growth by light + fg --------------------------------
fg_colors <- c(fg1 = "#BFE046", fg2 = "#27408b" , fg3 = "#267038", fg4 =  "#87Cefa", fg5 = "gray70") #unclassified = 'brown',#,  alltree = 'black')

#This one isn't working!!!

ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_labels)) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), 
              aes(x = light_area, ymin = q025, ymax = q975, group = fg, color = NA, fill = fg),
              alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), 
            aes(x = light_area, y = q50, group = fg, color = fg), size=0.5) +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q25, yend = q75), size = 0.3) +
  geom_segment(data = cast_pars %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified')), 
               aes(x = x_max_q50 * 0.5, xend = x_max_q50 * 2, y = y_max_q50 * 0.5, yend = y_max_q50 * 2), color = 'brown1', size = .5) +
  geom_point(shape=21, aes(x = bin_midpoint, y = median)) +
  scale_color_manual(values = guild_fills) +
  scale_fill_manual(values = fg_colors )+
  scale_x_log10(name = title_x, breaks=c(1,10,100)) + 
  scale_y_log10(name = title_y, labels=signif, breaks=c(0.01,0.1,1,10)) +
  theme_plant + theme_facet()

# Remove all tree and unclassified groups
param_ci$fg <- factor(param_ci$fg ,levels = c("fg1", "fg2", "fg5", "fg3", "fg4"))
fg_labels2 <- c("Fast", "Tall", "Medium", "Slow", "Short")

#---------------------------- Max Growth Rate with light -----------------------------------
ggplot(param_ci %>% filter(fg != 'NA', year == year_to_plot, parameter %in% 'log_slope', !fg %in% c('alltree','unclassified')),
            aes(x = fg, y = q50, ymin = q025, ymax = q975)) + 
  geom_errorbar(width = 0.4) + geom_point(size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_x_discrete(name = 'Life History Strategy', labels = fg_labels2) +
  scale_y_continuous(expression(paste('Max. Growth Rate (kg yr'^-1, ' m'^-2,')')), 
                     limits = c(0.6, 1.1),
                     breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2)) +
  theme_plant() + theme(aspect.ratio = 0.75)

#------------------- Maximum slope of growth with light --------------------------------


# ------------------------   WAIC of Piecewise Models  -----------------------------------

# This section was edited by QDR, 20 Jun 2019, for the updated model fits.

piece <- file.path(gdrive_path, 'data/data_piecewisefits')
ics <- read.csv(file.path(piece, 'piecewise_ics_by_fg.csv'), stringsAsFactors = FALSE)
ics$fg <- factor(ics$fg , labels = c("All", "Fast", "LL Pioneer", "Slow", "SL Breeder", "Medium", "Unclassified"))

# Density model
base_size <- 11
ggplot(ics %>% filter(criterion == 'WAIC', 
                           variable == 'density', !fg %in% 'Unclassified'), 
            aes(x = factor(dens_model), y = IC_value, ymin = IC_value - IC_stderr, ymax = IC_value + IC_stderr)) +
  facet_wrap(~ fg, labeller = label_value, scales = 'free_y') +
  geom_pointrange() +
  theme_bw(base_size = base_size, base_family = "",
           base_line_size = base_size/22, base_rect_size = base_size/11)+
  theme(strip.background = element_blank(),panel.grid = element_blank(), 
        text = element_text(size = 14),strip.text = element_text(size=12)) +
  labs(x = 'Segments in Density Function', y = "Widely Applicable Information Criterion (WAIC)")



# Growth model

ggplot(ics %>% filter(criterion == 'WAIC', 
                           variable == 'production', !fg %in% 'Unclassified'), 
            aes(x = factor(prod_model), y = IC_value, ymin = IC_value - IC_stderr, ymax = IC_value + IC_stderr)) +
  facet_wrap(~ fg, labeller = label_value,scales = 'free_y') +
  geom_pointrange() + #theme_facet +
  theme_bw(base_size = base_size, base_family = "",
           base_line_size = base_size/22, base_rect_size = base_size/11)+
  theme(strip.background = element_blank(),panel.grid = element_blank(), 
        text = element_text(size = 14),strip.text = element_text(size=12)) +
  labs(x = 'Segments in Growth Function',  y = "Widely Applicable Information Criterion (WAIC)")

#--------------------- Growth by LH Hex ----------------------------
# Get only year, func group, dbh, and production (no more is needed to plot right now)
raw_prod <- do.call(rbind, map2(alltreedat, seq(1985,2010,5), function(x, y) cbind(year = y, x %>% filter(!recruit) %>% select(fg, dbh_corr, production))))

raw_prod <- raw_prod %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg', fg), 'unclassified'))


# Original grayscale

hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias=1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1,10,100,1000,10000), 
                                             labels = c(1,10,100,1000,10000), limits=c(1,10000))

p <- plot_prod_withrawdata(year_to_plot = 1995,
                            fg_names = c('fg1','fg2','fg3','fg4','fg5','unclassified'),
                            full_names = c('Fast', 'Tall', 'Slow', 'Short', 'Medium', 'Unclassified'),
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
# ------------------- Total growth with height -------------------------
p0 <- plot_totalprod(year_to_plot = 1995,
                     fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                     model_fit_density = DENS, 
                     model_fit_production = PROD,
                     x_limits = c(0.5,230),
                     y_limits = c(0.5, 200),
                     y_breaks = c(0.1, 1, 10, 100),
                     y_labels = c(0.1, 1, 10, 100),
                     preddat = fitted_totalprod)
p <- p0 + scale_x_log10(name = 'Diameter (cm)',
                        breaks = c(1,3,10,30,100,300), limits = c(0.9, 230),
                        sec.axis = sec_axis(~ gMM(., a = 57.17, b = 0.7278, k = 21.57),
                                            name = "Height (m)", breaks = c(2, 3, 5, 10, 20, 40))) +
  scale_y_log10(position = "left", limits = c(0.5, 200), breaks = c(0.1, 1, 10, 100),labels = c(0.1, 1, 10, 100),
                name =expression(atop('Total Production', paste('(kg ha'^-1,' cm'^-1,' yr'^-1,')')))) +
  theme(aspect.ratio = 0.8) +geom_point(size = 1)

p


#----------------------- total production range --------------------------
minmax_prod_bycensus <- obs_totalprod %>%
  filter(bin_value > 0, bin_count >= 20, !fg %in% 'unclassified') %>%
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



alpha_value <- 1
hexfill2 <- scale_fill_gradient(low = 'forestgreen', high = 'navy', trans = 'log', breaks = c(1,3,10,30,100,300))
exl <- expression(atop('Intercepted Light', paste('per Individual (W)')))
exd <- 'Diameter (cm)'

labels = trans_format("log10", math_format(10^.x))

#----------------------- individual light capture scaling --------------------------

alpha_value <- 1
hexfill2 <- scale_fill_gradient(low = 'forestgreen', high = 'navy', trans = 'log', breaks = c(1,3,10,30,100,300))
exl <- expression(atop('Intercepted Light', paste('per Individual (W)')))
exd <- 'Diameter (cm)'

labels = trans_format("log10", math_format(10^.x))

p <- ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received)) +
  geom_pointrange(data = unscaledlightbydbhcloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl, breaks = c(1,100,10000, 1000000), limits = c(1,1000000), 
                labels = trans_format("log10", math_format(10^.x))) +
  theme_plant() + theme(legend.position = "right", legend.text = element_text(size = 15), legend.title = element_text(size = 16))+
  hex_scale_log_colors +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))# +
p

#----------------------- Scaling of total crown volume --------------------------


totalvolbins_fg <- read.csv(file.path(fp, 'obs_totalvol.csv'), stringsAsFactors = FALSE)

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
                     y_name = expression(atop('Total Crown Volume',paste('(m'^3, ' cm'^-1,' ha'^-1,')'))), 
                     preddat = fitted_totalvol,
                     obsdat = totalvolbins_fg, 
                     plot_abline = FALSE,
                     geom_size = 4)
p

