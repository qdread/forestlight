# Plot slopes of total production (energy equivalence slope)
# Edit 15 Aug: path to new data

# Change density model to 'pareto' and/or production model to 'powerlaw' if you want to plot those.
density_model <- 'weibull'
production_model <- 'powerlawexp'
year_to_plot <- 1990

fp <- '~/google_drive/ForestLight/data/data_forplotting_aug2018' ## CHANGE PATH AS NEEDED

fitted_slopes_by_fg <- read.csv(file.path(fp, 'fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)

# Plot total production slope

library(dplyr)
library(cowplot)

colors <- c('black', RColorBrewer::brewer.pal(5, 'Set1'))

# Just show total production slopes for all sizes
fitted_slopes_by_fg %>%
  filter(variable == 'total_production', dens_model == density_model, prod_model == production_model, year == year_to_plot, !fg %in% 'unclassified') %>%
ggplot(aes(x=dbh, y=q50, ymin=q025, ymax=q975, group=fg, color=fg, fill=fg)) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  scale_x_log10() + scale_y_continuous(limits = c(-10, 4)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors)
  

# Just total production slopes, Show intermediate sizes only
fitted_slopes_by_fg %>%
  filter(variable == 'total_production', dens_model == density_model, prod_model == production_model, year == year_to_plot, !fg %in% 'unclassified') %>%
ggplot(aes(x=dbh, y=q50, ymin=q025, ymax=q975, group=fg, color=fg, fill=fg)) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  scale_x_log10(limits = c(5,50), expand=c(0,0), breaks = c(5,10,20,50), name = 'Diameter (cm)') + 
  scale_y_continuous(limits = c(-10, 2), name = 'Total growth slope') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors)

# Slopes of density and growth too ----------------------------------------

fitted_slopes_by_fg %>%
  filter(dens_model == density_model, prod_model == production_model, year == year_to_plot, !fg %in% 'unclassified') %>%
ggplot(aes(x=dbh, y=q50, ymin=q025, ymax=q975, color=fg, fill=fg, group = interaction(fg,variable), linetype = variable)) +
  geom_line() +
  scale_x_log10(limits = c(1.1, 200), expand=c(0,0), breaks = c(2,5,10,20,50,100), name = 'Diameter (cm)') + 
  scale_y_continuous(limits = c(-10, 4), name = 'Slope', expand = c(0,0)) +
  geom_hline(yintercept = c(-2, 0, 2), linetype = 'dashed') +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors)

# Just density slope
fitted_slopes_by_fg %>%
  filter(variable == 'density', dens_model == density_model, prod_model == production_model, year == year_to_plot, !fg %in% 'unclassified') %>%
ggplot(aes(x=dbh, y=q50, ymin=q025, ymax=q975, group=fg, color=fg, fill=fg)) +
  geom_ribbon(alpha = 0.5) +
  geom_line(size = 1) +
  scale_x_log10(expand=c(0,0), breaks = c(5,10,20,50), name = 'Diameter (cm)') + 
  scale_y_continuous(name = 'Slope', expand = c(0,0)) +
  coord_cartesian(xlim = c(5,50), ylim=c(-3.3,-0.5)) +
  geom_hline(yintercept = c(-2), linetype = 'dashed') +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors)
