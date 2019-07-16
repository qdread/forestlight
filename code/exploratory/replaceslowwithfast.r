# Calculate alternative biomass if slow individuals were replaced with fast

library(tidyverse)

agb.allometry <- function(wsg,dbh,H) {0.0509*wsg*((dbh)^2)*H}

gdrive_path <- '~/google_drive/ForestLight'

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

dat1995 <- alltreedat[[3]] %>%
  select(sp, fg, wsg, height_corr, dbh_corr, agb_corr) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

# WSG by individual (show plot)
ggplot(dat1995, aes(x = factor(fg), y = wsg)) + 
  geom_boxplot() + 
  labs(x = 'Functional group', y = 'Wood specific gravity') +
  theme_bw()
  
# WSG by species (show plot)
dat1995 %>%
  group_by(sp) %>%
  slice(1) %>%
  ungroup %>%
  ggplot(aes(x = factor(fg), y = wsg)) +
    geom_boxplot() + 
    geom_jitter() +
    labs(x = 'Functional group', y = 'Wood specific gravity') 

# Replace the slow ones with average fast one (weighted median)
wsg_by_fg <- dat1995 %>%
  group_by(fg) %>%
  summarize(median_wsg = median(wsg))

# Slow have WSG of 0.609, fast have WSG of 0.393
fast_wsg <- wsg_by_fg$median_wsg[wsg_by_fg$fg == 'fg1']

dat1995 <- dat1995 %>%
  mutate(wsg_alternative = if_else(fg == 'fg3', fast_wsg, wsg),
         agb_alternative = agb.allometry(wsg_alternative, dbh_corr, height_corr))

# Create a plot of the old and new biomass

biomass_sums <- dat1995 %>%
  group_by(fg) %>%
  summarize(biomass_true = sum(agb_corr),
            biomass_alternative = sum(agb_alternative)) 

# Rename FG3 to FG1 in alternative scenario
biomass_sums$biomass_alternative[1] <- biomass_sums$biomass_alternative[1] + biomass_sums$biomass_alternative[3]
biomass_sums$biomass_alternative[3] <- 0

# Create plot
fgnames <- c('fast', 'pioneer', 'slow', 'breeder', 'medium', 'unclassified')
biomass_sums %>%
  mutate(fg = factor(fgnames, levels = fgnames)) %>%
  gather(scenario, biomass, -fg) %>%
  ggplot(aes(x = scenario, y = biomass, fill = fg)) +
    geom_col(position = 'stack') +
    theme_bw() +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(5, 'Set1'), 'black'))

sum(biomass_sums$biomass_alternative) / sum(biomass_sums$biomass_true)
