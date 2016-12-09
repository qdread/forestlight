# Clustering of Wright et al. growth rate and mortality data from BCI.

fp <- 'C:/Users/Q/Google Drive/ForestLight'

library(XLConnect)
wright <- readWorksheetFromFile(file = file.path(fp, 'data/Shade Tolerance/Demographic/Wright et al 2010, growth mortality tradeoffs.xlsx'), sheet = 1, startRow = 26) # Get rid of the lines above header.
wright[wright == -99] <- NA # Unknown values were given a value of -99

library(ggplot2)

ptradeoff <- ggplot(wright, aes(x = MRT25SAP/100, y = RGR95SAP)) + geom_point() + theme_minimal()

ptradeoff
ptradeoff + scale_y_log10()

library(scales)
ptradeoff + coord_trans(x = 'logit', y = 'log10')
ptradeoff + coord_trans(x = 'probit', y = 'log10')

# Complementary log-log
cloglog_trans <- function() trans_new('cloglog', function(x) log(-log(1 - x)), function(x) 1 - exp(-exp(x)))
ptradeoff + coord_trans(x = 'cloglog', y = 'log10')

# k-means clustering

dat <- wright[,c('MRT25SAP', 'RGR95SAP')]
dat <- dat[complete.cases(dat), ]

set.seed(37917)
k_raw <- kmeans(x = dat, centers = 3, nstart = 25)
k_logrgr <- kmeans(x = transform(dat, RGR95SAP = log10(RGR95SAP)), centers=3, nstart=25)
k_bothtrans <- kmeans(x = transform(dat, MRT25SAP = qlogis(MRT25SAP/100), RGR95SAP = log10(RGR95SAP)), centers=3, nstart=25)

# Swap 1 and 3 in k_bothtrans
tmp <- k_bothtrans$cluster
tmp[k_bothtrans$cluster == 1] <- 3
tmp[k_bothtrans$cluster == 3] <- 1

clust_all <- cbind(dat, clust_raw = k_raw$cluster, clust_logrgr = k_logrgr$cluster, clust_bothtrans = tmp)

pclustraw <- ggplot(clust_all, aes(x = MRT25SAP/100, y = RGR95SAP, color = factor(clust_raw))) + geom_point() + theme_minimal()
pclustlog <- ggplot(clust_all, aes(x = MRT25SAP/100, y = RGR95SAP, color = factor(clust_logrgr))) + geom_point() + theme_minimal()
pclustboth <- ggplot(clust_all, aes(x = MRT25SAP/100, y = RGR95SAP, color = factor(clust_bothtrans))) + geom_point() + theme_minimal()

library(gridExtra)
grid.arrange(pclustraw,pclustlog, pclustboth, nrow=1)

# Notes: the results appear to be the same, regardless of the initial random seed. The log transformation does not appear to make any difference, but transforming the mortality with logit does make a difference. It makes the intermediate cluster larger, cutting into the shade-tolerant cluster, but the pioneer cluster also takes over some of the intermediates in that case.
