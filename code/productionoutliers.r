# Check big tree production discrepancies.

bigprod <- alltreedat[[6]][alltreedat[[6]]$production > 1000, ]
bigprodid <- bigprod$treeID

bigtree2005 <- alltreedat[[5]][alltreedat[[5]]$treeID %in% bigprodid, c('treeID', 'sp', 'dbh_corr', 'agb_corr', 'production')]
bigtree2010 <- alltreedat[[6]][alltreedat[[6]]$treeID %in% bigprodid, c('treeID', 'sp', 'dbh_corr', 'agb_corr', 'production')]

library(dplyr)

bigtrees <- full_join(bigtree2005, bigtree2010, by = c('treeID','sp'))
bigtrees$dbh_corr.y - bigtrees$dbh_corr.x
