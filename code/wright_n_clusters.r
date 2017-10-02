# Determine optimal number of clusters

library(dplyr)

wright <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/wright2010.csv', stringsAsFactors = FALSE)
wright[wright == -99] <- NA # Unknown values were given a value of -99
wright$SPECIES.[109:110] <- c('simplex_var1', 'simplex_var2') # Correct duplicate named subspecies.
wright$Taxon <- with(wright, paste(GENUS., SPECIES.))

wright_df <- with(wright, data.frame(Taxon, mrt = MRT25SAP/100, rgr = RGR95SAP, stringsAsFactors = FALSE))
wright_df <- subset(wright_df, !is.na(mrt))

wright_pca <- with(wright_df, prcomp(data.frame(qlogis(mrt), log10(rgr)), scale=TRUE, center=TRUE)) # 90% of variation on the growth-mortality single axis. Nice.
pca_scores <- wright_pca$x[,1]

# Use nbClust package
library(NbClust)
wright_df <- transform(wright_df, logit_mrt = qlogis(mrt), log_rgr = log10(rgr))

nb_wright <- NbClust(wright_df[,c('logit_mrt','log_rgr')], distance = 'euclidean', min.nc = 2, max.nc = 5, method = 'kmeans', index = 'all', alphaBeale = 0.1)

nclusts <- nb_wright$Best.nc[1,]
nclusts <- nclusts[nclusts > 0 & !is.na(nclusts)]

ggplot(data.frame(n = nclusts), aes(x=n)) + geom_dotplot()

hist(nclusts, breaks = max(nclusts))

set.seed(4930116)
k_bothtrans <- kmeans(x = wright_df[,c('logit_mrt','log_rgr')], centers=2, nstart=25)
wright_df$cluster <- k_bothtrans$cluster
wright_df$pca <- pca_scores

ggplot(wright_df, aes(x = mrt, y = rgr, color = factor(cluster))) + geom_point()
