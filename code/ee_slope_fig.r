# Figure showing the energy equivalence slopes (sum of density and production power law slopes)
# Using only trees 5 to 50 cm in size

# Will add the model comparisons soon, as well as plots showing this.

param_cis_midsize <- read.csv('C:/Users/Q/google_drive/ForestLight/data/summarytables_12apr2018/paramci_eeslopes_midsizetrees.csv', stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)

param_cis_midsize %>%
  filter(year == 1995, !fg %in% 'unclassified', parameter %in% c('alpha','beta1','ee_slope')) %>%
  mutate_at(vars(starts_with('q')), funs(if_else(parameter=='alpha', -(.+1), .))) %>%
ggplot(aes(x = fg, y = q50, ymin = q025, ymax = q975, color = parameter, group = parameter)) +
  geom_errorbar(color='black', width = 0.1) +
  geom_point() +
  geom_hline(yintercept = c(-2,0,2), linetype = 'dotted') +
  theme_classic() + theme(panel.border = element_rect(fill=NA))
