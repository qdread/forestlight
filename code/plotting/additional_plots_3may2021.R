# New plots


# Coefficient symmetry plot -----------------------------------------------

# This code was modified from the code under "Fig 5: Symmetry plots" starting appx. line 1050 of clean_plotting_code.R

fitted_indivlight <- read.csv(file.path(fp, 'fitted_indivlight.csv'), stringsAsFactors = FALSE)
fitted_totallight <- read.csv(file.path(fp, 'fitted_totallight.csv'), stringsAsFactors = FALSE)
indivlightbins_fg <- read.csv(file.path(fp, 'obs_indivlight.csv'), stringsAsFactors = FALSE)
totallightbins_fg <- read.csv(file.path(fp, 'obs_totallight.csv'), stringsAsFactors = FALSE)

fitted_indivlightcaptured <- read.csv(file.path(fp, 'fitted_indivlightcaptured.csv'), stringsAsFactors = FALSE)
fitted_totallightcaptured <- read.csv(file.path(fp, 'fitted_totallightcaptured.csv'), stringsAsFactors = FALSE)
indivlightcapturedbins_fg <- read.csv(file.path(fp, 'obs_indivlightcaptured.csv'), stringsAsFactors = FALSE)
totallightcapturedbins_fg <- read.csv(file.path(fp, 'obs_totallightcaptured.csv'), stringsAsFactors = FALSE)


# Fig 5a      

# Plot
grob_text <- grobTree(textGrob("Solar Equivalence", x = 0.27, y = 0.87, hjust = 0,
                               gp = gpar(col = "gold3", fontsize = 18))) 

grob_a <- grobTree(textGrob("a", x = 0.06, y = 0.91, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
totallightcapturedbins_fg <- totallightcapturedbins_fg %>%
  filter(bin_count >= 20)
tot_light <- plot_totalprod(year_to_plot = 1995,
                            fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                            model_fit_density = DENS, 
                            model_fit_production = PROD,
                            x_limits = c(0.9,200),
                            y_limits = c(100, 200000),
                            geom_size = 3.5,
                            y_breaks = c(100, 1000, 10000, 100000),
                            x_breaks = c(1, 10, 100),
                            y_labels = c("0.1", "1", "10", "100"),
                            y_name = expression(paste('Total Light Captured (kW cm'^-1,' ha'^-1,')')),
                            preddat = fitted_totallightcaptured,
                            obsdat = totallightcapturedbins_fg,
                            plot_abline = FALSE)


tot_light2 <- tot_light  + 
  scale_y_continuous(position = "right", trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("0.1", "1", "10", "100"), 
                     limits = c(200, 200000),
                     name = expression(atop('Total Light Captured',paste('(kW cm'^-1,' ha'^-1,')'))))  +
  scale_x_log10(name = "Stem Diameter (cm)", limits = c(0.8, 200), position = "top") +
  theme(aspect.ratio = 0.75) + 
  geom_abline(intercept = log10(70000), slope = 0, color ="#C9B074",
              linetype="dashed", size=.75) +
  annotation_custom(grob_text) #+ annotation_custom(grob_a)
plot(tot_light2)
p_tot_light <- tot_light2
# ---------------Fig 5b--------------

params <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/piecewise_paramci_by_fg.csv'), stringsAsFactors = FALSE)

xvalues <- params %>% 
  filter(variable == 'density', model == DENS) %>%
  group_by(fg) %>%
  summarize(mid_cutoff = (mean[parameter == 'tau_high'] + mean[parameter == 'tau_low'])/2)


growth_slopes <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
light_slopes <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/lightcaptured_piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
growth_slopes_atmiddle <- growth_slopes %>% 
  filter(variable == 'total_production', dens_model == DENS, prod_model == PROD) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

light_slopes_atmiddle <- light_slopes %>% 
  filter(variable == 'total_captured_light', dens_model == DENS, prod_model == PROD) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

fg_full_names <- c('Fast', 'Tall', 'Slow', 'Short', 'Medium', 'All', 'Unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

allslopes <- rbind(growth_slopes_atmiddle, light_slopes_atmiddle) %>%
  ungroup %>%
  mutate(fg = factor(fg, levels = fgs, labels = fg_full_names))
grobb <- grobTree(textGrob("b", x = 0.04, y = 0.9,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob1 <- grobTree(textGrob("Solar", x = 0.68, y = 0.94, hjust = 0,
                           gp = gpar(col = "gold3", fontsize = 18))) 
grob2 <- grobTree(textGrob("Production", x = 0.68, y = 0.86, hjust = 0,
                           gp = gpar(col = "darkgreen", fontsize = 18)))
grob3 <- grobTree(textGrob("Energy Equivalence", x = 0.24, y = 0.51, hjust = 0,
                           gp = gpar(col = "black", fontsize = 18))) 




#---------------------Fig 5b--------------
# compare slopes
slopes <- ggplot(allslopes %>% filter(!fg %in% 'Unclassified'), 
                 aes(x = fg, y = q50, ymin = q025, ymax = q975, fill =  variable, color =variable)) +
  theme_plant() +
  geom_hline(yintercept = 0, linetype = 'dashed', size = .75) +
  geom_errorbar(position = position_dodge(width = 0.6), size = 0.75, width = 0) +
  geom_point(position = position_dodge(width = 0.6), shape = 21, size = 3.5, color = "black", stroke = 0.5) +
  geom_errorbar(data = allslopes %>% filter(fg %in% c('Tall', 'Slow', 'All')), 
                position = position_dodge(width = 0.6), size = 0.75, width = 0) +
  scale_y_continuous(position = "right", ) +
  labs(x = NULL) +
  scale_y_continuous(name = "Scaling Slope", position = "right", limits = c(-1.05, 1.3)) +
  scale_fill_manual(values = c('gold1', 'darkolivegreen3')) +
  scale_color_manual(values = c('gold3', 'darkgreen')) +
  theme(axis.text.x = element_text(angle = 25,  vjust = .7, face = "italic", size = 18)) +
  annotation_custom(grob1) + annotation_custom(grob2) + annotation_custom(grob3) #+annotation_custom(grob0)
slopes

g_slopes <- ggplotGrob(slopes)
g_tot_light <- ggplotGrob(p_tot_light)
combo <- rbind(g_tot_light, g_slopes, size = "first")
combo$widths <- unit.pmax(g_tot_light$widths,g_slopes$widths)
grid.newpage()
grid.draw(combo)


g_hex <- ggplotGrob(light_hex)
g_light <- ggplotGrob(combo)
new <- cbind(light_hex,combo, size = "first")

combo <- cbind(g_tot_light, g_slopes, size = "first")
new$heights <- unit.pmax(new$widths, new$widths)
ggsave(combo, height = 8.6, width = 6, filename = file.path(gdrive_path,'Figures/Symmetry/combo3.pdf'))
grid.draw(new)
dev.off()

system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path2,'Figures/Symmetry/light_combo2.pdf'), 
                    file.path(gdrive_path2,'Figures/Symmetry/light_combo2.pdf')) 
)



# Total leaf area scaling plot --------------------------------------------

# This was modified based on the volume scaling plot at around line 3387 of clean_plotting_code.R

# Fitted values for individual light, total light, and total volume

fitted_totalleafarea <- read.csv(file.path(fp, 'fitted_totalleafarea.csv'), stringsAsFactors = FALSE)
totalleafareabins_fg <- read.csv(file.path(fp, 'obs_totalleafarea.csv'), stringsAsFactors = FALSE)



totalleafareabins_fg <- totalleafareabins_fg %>%
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
                     y_name = expression(paste('Total Leaf Area (m'^2, ' cm'^-1, ' ha'^-1,')')), 
                     preddat = fitted_totalleafarea,
                     obsdat = totalleafareabins_fg, 
                     plot_abline = FALSE,
                     geom_size = 4)
p
p0 <- p + scale_y_continuous(position = "left", trans = "log10", limits = c(9, 5000),
                             name = expression(atop('Total Leaf Area',paste('(m'^2, ' cm'^-1,' ha'^-1,')')))) +
  theme_plant_small(legend = TRUE)  + guide +
  scale_fill_manual(values = guild_fills2, labels = guild_labels2) 
plot(p0)



# Captured light hexbin plot ----------------------------------------------

# This plot was modified from the plot on line 3291 of clean_plotting_code.R, using the captured light data instead of incoming light.

ex_lightcaptured <- expression(atop("Captured Light", paste("per Individual (W)")))

indiv_light_captured <- ggplot() +
  theme_plant() +
  scale_x_log10(name = exd, limits = c(.9, 200)) +
  scale_y_log10(name = ex_lightcaptured, breaks = c(1,100,10000, 1000000), limits = c(1,1000000), 
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_captured)) +
  hex_scale_log_colors +
  geom_pointrange(data = unscaledlightcapturedbydbhcloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = mean, ymin = mean, ymax = mean), size = .5) +
  theme(legend.position = "right", legend.text = element_text(size = 15), legend.title = element_text(size = 16))+
  hex_scale_log_colors +
  geom_smooth(data = alltree_light_95, aes(x = dbh_corr, y = light_captured),
              method = "lm", color = "black", size = .7) +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))
