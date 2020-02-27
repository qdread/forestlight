
library(tidyverse)

tprod <- read_csv('~/google_drive/ForestLight/data/data_binned/totalproductionbin_byyear.csv')

prodsums <- tprod %>%
  group_by(year, fg) %>%
  summarize(total_production = sum(bin_value),
            n_individuals = sum(bin_count))

write_csv(prodsums, '~/google_drive/ForestLight/data/data_binned/sum_totalproduction_1995.csv')
