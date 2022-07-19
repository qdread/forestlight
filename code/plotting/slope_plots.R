# change these if needed.
DENS = 3
PROD = 1

# Set path to data on google drive
#devtools::install_github('qdread/forestscaling')

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Volumes/GoogleDrive/My Drive/ForestLight'))
github_path <- ifelse(Sys.info()['user'] == 'qread', '~/Documents/GitHub/MSU_repos', file.path('/Users/jgradym/Documents/GitHub/'))


gdrive_path2 <-  file.path('/Users/jgradym/Google\\ Drive/ForestLight')
gdrive_path2 <-  file.path('/Volumes/GoogleDrive/My\\ Drive/ForestLight')

library(broom)
library(egg)
library(scales)
library(RColorBrewer)
library(gtable)
library(grid)
library(reshape2)
library(hexbin)
library(forestscaling) # Packaged all the functions and ggplot2 themes here!
library(tidyverse)
#devtools::install_github('qdread/forestscaling')

# Define color schemes and labels
guild_colors<- c("fast" = "#BFE046", "large pioneer" = "#267038", "slow" = "#27408b", "small breeder" = "#87Cefa", medium = "gray50")
guild_fills2 <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")
guild_fills <- c("fast" = "#BFE046", "large pioneer" = "#267038", "slow" = "#27408b", "small breeder" = "#87Cefa", medium = "gray93")

guild_colors <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")
guild_colors2 <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray")


biomass_growth = read_csv('/Volumes/GoogleDrive/.shortcut-targets-by-id/0Bzy2GmZ-I6IcT0JmNk96Sl9iMVU/ForestLight/data/clean_summary_tables/clean_parameters_growth.csv')

# small tree slope
ggplot(data = biomass_growth %>% filter(parameter_description == "slope small trees", fg != "all trees", model == "two segment"),
       aes(x = fg, y = q50)) +
  geom_point(shape = 21, size = 2) + 
  geom_errorbar(aes(ymin = q025, ymax = q975)) +
  theme_plant()

# all tree biomass slope

ggplot(data = biomass_growth %>% filter(parameter_description == "slope", fg != "all trees", model == "one segment"),
       aes(x = fg, y = q50, color = fg, fill = fg)) +
  geom_point(size = 5) + 
  scale_color_manual(values = guild_colors) +
  geom_errorbar(aes(ymin = q025, ymax = q975)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")+
  scale_y_continuous(limits = c(2.1, 2.4), name = "Slope of Biomass Growth", position  = "right") +
  theme_plant()

#------------Diameter Growth ---------------------
diameter_growth = read_csv("/Volumes/GoogleDrive/.shortcut-targets-by-id/0Bzy2GmZ-I6IcT0JmNk96Sl9iMVU/ForestLight/data/clean_summary_tables/clean_parameters_individualdiametergrowth.csv")

ggplot(data = diameter_growth %>% filter(parameter_description == "slope", fg != "all trees", model == "one segment"),
       aes(x = fg, y = q50, color = fg, fill = fg)) +
  geom_point(size = 5) + 
  theme_plant() +
  scale_color_manual(values = guild_colors) +
  geom_errorbar(aes(ymin = q025, ymax = q975)) +
  scale_y_continuous( name = "Slope of Diameter Growth", position = "right") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")+
  scale_x_discrete(name = NULL)
