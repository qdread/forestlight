h_bohl <- function(dbh) 10^(.438 + .595 * log10(dbh))    # Height
library(forestscaling)
library(tidyverse)

ggplot(data.frame(dbh=c(1,300)), aes(x = dbh)) +
  stat_function(fun = h_bohl) +
  stat_function(fun = taper.H, color = 'red') +
  stat_function(fun = notaper.H, color = 'blue') +
  scale_x_log10() + scale_y_log10()
  