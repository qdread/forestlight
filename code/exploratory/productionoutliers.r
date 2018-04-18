# Check big tree production discrepancies.

bigprod <- alltreedat[[6]][alltreedat[[6]]$production > 1000, ]
bigprodid <- bigprod$treeID

bigtree2005 <- alltreedat[[5]][alltreedat[[5]]$treeID %in% bigprodid, c('treeID', 'sp', 'dbh_corr', 'agb_corr', 'production')]
bigtree2010 <- alltreedat[[6]][alltreedat[[6]]$treeID %in% bigprodid, c('treeID', 'sp', 'dbh_corr', 'agb_corr', 'production')]

library(dplyr)

bigtrees <- full_join(bigtree2005, bigtree2010, by = c('treeID','sp'))
bigtrees$dbh_corr.y - bigtrees$dbh_corr.x
table(bigtrees$dbh_corr.y - bigtrees$dbh_corr.x > 20)

# Plot all diameter increments
alltree2005 <- alltreedat[[5]][, c('treeID', 'sp', 'dbh_corr', 'agb_corr', 'production')]
alltree2010 <- alltreedat[[6]][, c('treeID', 'sp', 'dbh_corr', 'agb_corr', 'production')]

alltreebothyears <- full_join(alltree2005, alltree2010, by = c('treeID','sp'))

ggplot(alltreebothyears, aes(x = dbh_corr.y - dbh_corr.x)) +
  geom_histogram() + scale_x_log10()

# plot diameter vs diameter
ggplot(alltreebothyears, aes(x = dbh_corr.x, y = dbh_corr.y - dbh_corr.x)) +
  geom_point() + scale_x_log10(breaks=c(1,2,5,10,20,50,100))+scale_y_log10(breaks=c(1,2,5,10,20,50))
  