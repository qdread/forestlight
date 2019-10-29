# Clustering of Wright et al. growth rate and mortality data from BCI.

fp <- 'C:/Users/Q/Google Drive/ForestLight'

library(XLConnect)
wright <- readWorksheetFromFile(file = file.path(fp, 'data/Shade Tolerance/Demographic/Wright et al 2010, growth mortality tradeoffs.xlsx'), sheet = 1, startRow = 26) # Get rid of the lines above header.
wright[wright == -99] <- NA # Unknown values were given a value of -99
wright$SPECIES.[109:110] <- c('simplex_var1', 'simplex_var2') # Correct duplicate named subspecies.
wright$taxon <- with(wright, paste(GENUS., SPECIES.))


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


# Run with pca

library(dplyr)

wright_complete <- wright %>% 
  select(taxon, MRT25SAP, RGR95SAP) %>% 
  filter(!is.na(MRT25SAP), !is.na(RGR95SAP)) %>% 
  mutate(MRT25SAP = MRT25SAP/100,
         log_rgr = log10(RGR95SAP),
         logit_mort = qlogis(MRT25SAP))

set.seed(48823)

pca_raw <- prcomp(x = wright_complete[, c('MRT25SAP', 'RGR95SAP')], scale. = TRUE, center = TRUE)
k_pca_raw <- kmeans(x = pca_raw$x[,'PC1'], centers = 3, nstart = 25)

pca_logrgr <- prcomp(x = wright_complete[, c('MRT25SAP', 'log_rgr')], scale. = TRUE, center = TRUE)
k_pca_logrgr <- kmeans(x = pca_logrgr$x[,'PC1'], centers = 3, nstart = 25)

pca_bothtrans <-  prcomp(x = wright_complete[, c('logit_mort', 'log_rgr')], scale. = TRUE, center = TRUE)
k_pca_bothtrans <- kmeans(x = pca_bothtrans$x[,'PC1'], centers = 3, nstart = 25)

# Swap clusters around so they match

tmpraw <- k_pca_raw$cluster
tmpraw[k_pca_raw$cluster == 2] <- 1
tmpraw[k_pca_raw$cluster == 3] <- 2
tmpraw[k_pca_raw$cluster == 1] <- 3


tmplog <- k_pca_logrgr$cluster
tmplog[k_pca_logrgr$cluster == 1] <- 2
tmplog[k_pca_logrgr$cluster == 2] <- 1


tmpboth <- k_pca_bothtrans$cluster
tmpboth[k_pca_bothtrans$cluster == 2] <- 3
tmpboth[k_pca_bothtrans$cluster == 3] <- 2


wright_complete <- mutate(wright_complete,
                          PC1raw = pca_raw$x[,'PC1'],
                          PC1logrgr = pca_logrgr$x[,'PC1'],
                          PC1bothtrans = pca_bothtrans$x[,'PC1'],
                          clust_raw = tmpraw,
                          clust_logrgr = tmplog,
                          clust_bothtrans = tmpboth)

save(wright_complete, file = 'kmeansoutput.r')

t1 <- theme_minimal() + theme(legend.position = 'bottom')

ppcaclustraw <- ggplot(wright_complete, aes(x = MRT25SAP/100, y = RGR95SAP, color = factor(clust_raw))) + geom_point() + t1
ppcaclustlog <- ggplot(wright_complete, aes(x = MRT25SAP/100, y = RGR95SAP, color = factor(clust_logrgr))) + geom_point() + t1
ppcaclustboth <- ggplot(wright_complete, aes(x = MRT25SAP/100, y = RGR95SAP, color = factor(clust_bothtrans))) + geom_point() + t1

library(gridExtra)
png('C:/Users/Q/Google Drive/ForestLight/figs/threeclusters.png', height=5, width=9, res=300, units='in')
grid.arrange(ppcaclustraw,ppcaclustlog, ppcaclustboth, nrow=1)
dev.off()