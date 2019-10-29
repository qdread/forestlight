# Plot mortality binned data (binned in the mortality_processing.r script)
# this version uses the correct mixed model fit lines
# QDR / Forestlight / 08 Oct 2019

# Load data ---------------------------------------------------------------

library(tidyverse)

# Check out this slick trick so that we no longer need to comment out any paths - just run this and it sees if it's Quentin or not.
user <- Sys.info()['user']
gdrive_path <- ifelse(user == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))
github_path <- ifelse(user == 'qread', '~/Documents/GitHub/forestlight', file.path('/Users',user,'Documents/GitHub/forestlight'))

bin_mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting_aug2018/obs_mortalitybins.csv')) # Load binned data
mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting_aug2018/obs_mortalityindividuals.csv')) # Load raw data
fitted_mort <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/mortality_ci_by_fg.csv')) # Load fitted values

# Make some simple plots --------------------------------------------------

# Color mapping

theme_plant <- theme(panel.grid = element_blank(), #for Total Production
                     aspect.ratio = .75,
                     axis.text = element_text(size = 15, color = "black"), 
                     axis.ticks.length=unit(0.2,"cm"),
                     axis.title = element_text(size = 15),
                     axis.title.y = element_text(margin = margin(r = 10)),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.x.top = element_text(margin = margin(b = 5)),
                     plot.title = element_text(size = 15, face = "plain", hjust = 10),
                     panel.border = element_rect(color = "black", fill=NA,  size=1),
                     panel.background = element_rect(fill = "transparent",colour = NA),
                     plot.background = element_rect(fill = "transparent",colour = NA),
                     legend.position = "none",
                     rect = element_rect(fill = "transparent"),
                     text = element_text(family = 'Helvetica')) 

fg_labels <- c('Fast','LL Pioneer', 'Slow', 'SL Breeder', 'Medium')
guild_fills_nb <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")
guild_colors_nb <- c("#3B4403", "#02330A", "#031a49", "#02394F", "gray")


(plightarea <- ggplot(data = fitted_mort %>% mutate(fg = factor(fg, labels = fg_labels))) +
  geom_ribbon(aes(x = light_per_area, ymin = q025, ymax = q975, group = fg, fill = fg), alpha = 0.3) +
  geom_line(aes(x = light_per_area, y = q50, group = fg, color = fg)) +
  geom_point(data = bin_mort %>% filter(variable == 'light_per_area', !fg %in% c('all','unclassified'), (lived+died) > 20)  %>% mutate(fg = factor(fg, labels = fg_labels)),
             aes(x = bin_midpoint, y = mortality, fill = fg), #lived + died),
             shape = 21, size = 3) +
  scale_x_log10(name = parse(text = 'Light~per~Crown~Area~(W~m^-2)'), limits = c(1.1, 412)) +
  scale_y_continuous(breaks = c(0.03, 0.3, .1), labels = c(0.03, 0.3, 0.1), limits = c(0.02, .5),
                name = expression(paste("Mortality (5 yr"^-1,")")), trans = 'logit') +
  #scale_size_continuous(name = 'Individuals', trans = 'log', breaks = 10^(0:3)) +
  scale_color_manual(values = guild_colors_nb) +
  scale_fill_manual(values = guild_fills_nb) +
  theme_plant) #+
  #ggtitle('Mortality rate by light per area')
plightarea
p2 <- set_panel_size(plightarea, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p2)

pdf(file.path(gdrive_path, "Figures/Fig_2/Fig_2c.pdf"))
grid.newpage()
grid.draw(p2)
dev.off()
