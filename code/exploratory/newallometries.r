library(tidyverse)
library(forestscaling)
library(broom)

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))

allo <- read_csv(file.path(gdrive_path, 'data/BCI_raw/allometry_bci_trees.csv'))

table(allo$species)
length(unique(allo$species)) # 80 species

# Cross reference allometry with FG lookup.
fgbci <- read.table(file.path(gdrive_path, 'data/Ruger/fgroups_dynamics_new.txt'), stringsAsFactors = FALSE)

# Correct functional groups so that: 1 fast, 2 pioneer, 3 slow, 4 breeder, 5 intermediate
# Old 1,2,3,4,5 --> New 2,3,1,4,5
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5))

fg_lookup <- data.frame(species = fgbci$splower, fg = fgbci$fg5)

allo <- allo %>%
  select(-X9) %>%
  left_join(fg_lookup)

# Unique species fg combinations
spp <- unique(allo[,c('species', 'fg')])

table(spp$fg, useNA = 'always')

# Visualize the relationships for each group

guild_colors_nb <- c("#3B4403", "#02330A", "#031a49", "#02394F", "#595A5B")
guild_fills_nb <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray93")

pdf(file.path(gdrive_path, 'figs/allometries_08nov2019.pdf'), height = 6, width = 6)

ggplot(allo %>% filter(!is.na(fg)), aes(x = dbh.cm, y = ht.m, color = factor(fg))) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  scale_x_log10(name = 'Diameter (cm)') + scale_y_log10(name = 'Height (m)') +
  scale_color_manual(values = guild_fills_nb) +
  theme_minimal()

ggplot(allo %>% filter(!is.na(fg)), aes(x = dbh.cm, y = crowndepth.m, color = factor(fg))) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  scale_x_log10(name = 'Diameter (cm)') + scale_y_log10(name = 'Crown depth (m)') +
  scale_color_manual(values = guild_fills_nb) +
  theme_minimal()

ggplot(allo %>% filter(!is.na(fg)), aes(x = dbh.cm, y = crownarea.m2, color = factor(fg))) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  scale_x_log10(name = 'Diameter (cm)') + scale_y_log10(name = 'Crown area (m2)') +
  scale_color_manual(values = guild_fills_nb) +
  theme_minimal()

dev.off()

# Calculate allometries for each variable for each species.
allo_coefs <- allo %>% 
  pivot_longer(ht.m:crownarea.m2, names_to = 'variable') %>%
  group_by(fg, species, variable) %>%
  group_modify(~ tidy(lm(log(value) ~ log(dbh.cm), data = .), conf.int = TRUE), keep = TRUE)

# Plot coefficients

pslopes <- ggplot(allo_coefs %>% filter(!term %in% "(Intercept)"), aes(x = interaction(fg, species, drop = TRUE, lex.order = TRUE), y = estimate, ymin = conf.low, ymax = conf.high, color = factor(fg))) +
  facet_wrap(~ variable, scales = 'free_y') +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_color_manual(values = guild_fills_nb) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  ggtitle('Slopes')

pintercepts <- ggplot(allo_coefs %>% filter(term %in% "(Intercept)"), aes(x = interaction(fg, species, drop = TRUE, lex.order = TRUE), y = estimate, ymin = conf.low, ymax = conf.high, color = factor(fg))) +
  facet_wrap(~ variable, scales = 'free_y') +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_color_manual(values = guild_fills_nb) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  ggtitle('Intercepts')

pdf(file.path(gdrive_path, 'figs/allometries_coefs_08nov2019.pdf'), height = 6, width = 9)
pslopes; pintercepts
dev.off()
