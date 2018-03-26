# Production by functional group by year

library(purrr)
library(dplyr)

prod_by_fg <- map(alltreedat[-1], function(x) x %>%
      group_by(fg) %>%
      summarize(sum_prod = sum(production, na.rm = TRUE)))

prod_by_fg <- cbind(year = rep(seq(1990,2010,5), each=6), do.call(rbind, prod_by_fg))
prod_by_fg$fg <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified')
prod_by_fg$prod_per_ha <- prod_by_fg$sum_prod/42.84
write.csv(prod_by_fg, file = 'C:/Users/Q/google_drive/ForestLight/data/production_by_fg.csv', row.names = FALSE)
