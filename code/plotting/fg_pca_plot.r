library(cowplot)
library(dplyr)

# Load Nadja's data (new functional groups 25 June)
# fg5 is the new column (we originally used fg from the older df)
fgbci <- read.table('~/google_drive/ForestLight/data/Ruger/fgroups_dynamics_new.txt', stringsAsFactors = FALSE)

# Correct functional groups so that: 1 fast, 2 pioneer, 3 slow, 4 breeder, 5 intermediate
# Old 1,2,3,4,5 --> New 2,3,1,4,5
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5))

# Currently X1new is high for slow species and low for fast species
# Currently X2new is high for pioneer species and low for breeder species
# Correct these
fgbci$PC_slow_to_fast <- -fgbci$X1new
fgbci$PC_breeder_to_pioneer <- fgbci$X2new


guild_colors <- RColorBrewer::brewer.pal(5, 'Set1')
fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('fast','long-lived pioneer', 'slow', 'short-lived breeder', 'intermediate')


# Plot functional groups
ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, color = factor(fg5))) +
  geom_point() +
  labs(x = 'X1 slow to fast', y = 'X2 breeders to pioneers') +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group')

# ggsave('C:/Users/Q/google_drive/ForestLight/figs/figures_june2018_newcluster/fg5plot (newest groups).pdf')