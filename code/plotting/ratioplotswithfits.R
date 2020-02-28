
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

# Load fitted values
ratio_fitted_diam <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/ratio_fittedvalues.csv'))
ratio_fitted_lightarea <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/ratio_fittedvalues_lightarea.csv'))

# ----------- alternative figure 6a with fit lines -------------------------------------

grob1_1a <- grobTree(textGrob("1:1", x = 0.5, y = 0.52,  hjust = 0,
                              gp = gpar(col = "black", fontsize = 19))) 

ratio_fitted_lightarea_prod <- ratio_fitted_lightarea %>%
  filter(variable == 'total production', light_area > 7) %>%
  mutate(ratio = if_else(ratio == 'fast:slow', 'Fast-Slow', 'Breeder-Pioneer'))

prod_ratio <- prod_ratio_light   %>%
  filter(n_individuals >= 20) %>%
  ggplot() +
  geom_ribbon(aes(x = light_area, ymin = q025, ymax = q975, group = ratio, fill = ratio), alpha = 0.5, data = ratio_fitted_lightarea_prod) +
  geom_line(aes(x = light_area, y = q50, group = ratio), data = ratio_fitted_lightarea_prod) +
  geom_point(aes(x = bin_midpoint, y = production_ratio, fill = ID), shape = 21, size = 4.5,  stroke = .5, color = "black")+
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

#----------------------------- Fig 6B Abundance by Diameter ----------------------------
grob1_1 <- grobTree(textGrob("1:1", x = 0.5, y = 0.54,  hjust = 0,
                             gp = gpar(col = "black", fontsize = 19))) 

grob_r <- grobTree(textGrob("Replacement", x = 0.55, y = 0.35,  hjust = 0,
                            gp = gpar(col = "black", fontsize = 20))) 

ratio_fitted_diam_density <- ratio_fitted_diam %>%
  filter(variable == 'density', (ratio == 'fast:slow' & dbh < 48) | (ratio == 'pioneer:breeder') & dbh < 15) %>%
  mutate(ratio = if_else(ratio == 'fast:slow', 'Fast-Slow', 'Breeder-Pioneer'))

dens_ratio <- prod_ratio_diam %>% 
  filter(density_ratio > 0) %>%
  filter(n_individuals >= 20) %>%
  ggplot() +
  geom_ribbon(aes(x = dbh, ymin = q025, ymax = q975, group = ratio, fill = ratio), alpha = 0.5, data = ratio_fitted_diam_density) +
  geom_line(aes(x = dbh, y = q50, group = ratio), data = ratio_fitted_diam_density) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  geom_point(aes(x = bin_midpoint, y = density_ratio, fill = ID), shape = 21, size = 4.5,  stroke = .5,  color = "black")+
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
