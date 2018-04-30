
# Change density model to 'pareto' and/or production model to 'powerlaw'
density_model <- 'weibull'
production_model <- 'powerlawexp'

fp <- 'C:/Users/Q/google_drive/ForestLight/data/data_forplotting_12apr2018' ## CHANGE PATH AS NEEDED

# Read all the csvs in directory.
for (i in dir(fp, pattern = '.csv')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}

# Plot total production slope

library(dplyr)
library(cowplot)

# Use formula d log y / d log x = x/y dy/dx to get slope of total production.
slope_totalprod <- pred_totalprod %>%
  filter(year == 1995, !fg %in% 'unclassified', dens_model == density_model, prod_model == production_model) %>%
  group_by(fg) %>%
  mutate_at(vars(starts_with('q')), funs(c(NA, (dbh[-1]/.[-1]) * diff(.)/diff(dbh))))

# Calculate slope of density  
slope_density <- pred_dens %>%
  filter(year == 1995, !fg %in% 'unclassified', dens_model == density_model, prod_model == production_model) %>%
  group_by(fg) %>%
  mutate_at(vars(starts_with('q')), funs(c(NA, (dbh[-1]/.[-1]) * diff(.)/diff(dbh))))

# Calculate slope of individual production
slope_indivprod <- pred_indivprod %>%
  filter(year == 1995, !fg %in% 'unclassified', dens_model == density_model, prod_model == production_model) %>%
  group_by(fg) %>%
  mutate_at(vars(starts_with('q')), funs(c(NA, (dbh[-1]/.[-1]) * diff(.)/diff(dbh))))

colors <- c('black', RColorBrewer::brewer.pal(5, 'Set1'))

# Just show total production slopes for all sizes
ggplot(slope_totalprod, aes(x=dbh, y=q50, ymin=q025, ymax=q975, group=fg, color=fg, fill=fg)) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  scale_x_log10() + scale_y_continuous(limits = c(-10, 4)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors)
  

# Just total production slopes, Show intermediate sizes only
ggplot(slope_totalprod, aes(x=dbh, y=q50, ymin=q025, ymax=q975, group=fg, color=fg, fill=fg)) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  scale_x_log10(limits = c(5,50), expand=c(0,0), breaks = c(5,10,20,50), name = 'Diameter (cm)') + 
  scale_y_continuous(limits = c(-10, 2), name = 'Total growth slope') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors)

# Slopes of density and growth too ----------------------------------------

ggplot(slope_totalprod, aes(x=dbh, y=q50, ymin=q025, ymax=q975, group=fg, color=fg, fill=fg)) +
  geom_line(size = 1.5) +
  geom_line(data = slope_indivprod) +
  geom_line(data = slope_density) +
  scale_x_log10(limits = c(5,50), expand=c(0,0), breaks = c(5,10,20,50), name = 'Diameter (cm)') + 
  scale_y_continuous(limits = c(-10, 4), name = 'Slope', expand = c(0,0)) +
  geom_hline(yintercept = c(-2, 0, 2), linetype = 'dashed') +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors)

