# Load John's macrosystem data and explore/visualize it.

harvard <- read.csv('data/harvard.csv', stringsAsFactors = FALSE)
bci <- read.csv('data/bci.csv', stringsAsFactors = FALSE)

library(dplyr)
harvard <- rename(harvard, CII=CI)
bci <- filter(bci, Site == 'Barro Colorado Island')

# Match the two data frames and combine.
names(harvard) %in% names(bci)
names(bci) %in% names(harvard)
commonnames <- intersect(names(harvard), names(bci))

macrodata <- rbind(harvard[, commonnames], bci[, commonnames])

# Correct typos
macrodata$CII[macrodata$CII == 1.2] <- 1.5
macrodata$CII[macrodata$CII == 25] <- 2.5

# Taxon by light availability
with(macrodata, table(Taxon, CII))

macrodata <- mutate(macrodata, 
                    dbhYear2_RGR = log(Year2_DBH) - log(Year1_DBH), 
                    dbhYear3_RGR = log(Year3_DBH) - log(Year2_DBH),
                    dghYear2_RGR = log(Year2_DGH) - log(Year1_DGH_1),
                    dghYear3_RGR = log(Year3_DGH) - log(Year2_DGH),
                    Year2_RGR = pmax(dbhYear2_RGR, dghYear2_RGR, na.rm = TRUE),
                    Year3_RGR = pmax(dbhYear3_RGR, dghYear3_RGR, na.rm = TRUE)
                    )

# Exclude negative growth rates
macrodata$Year2_RGR[macrodata$Year2_RGR < 0] <- NA
macrodata$Year3_RGR[macrodata$Year3_RGR < 0] <- NA


# Sums of species

nwithcii <- macrodata %>% group_by(Taxon) %>% summarize(ncii = sum(!is.na(CII)))
macrodata <- left_join(macrodata, nwithcii)

library(ggplot2)

simpletheme <- theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank())

# To plot number of observations
n_above <- function(x) {
  ypos <- ifelse(sum(!is.na(x)) == 1, 1.2*max(x, na.rm=TRUE), 1.2*(mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))))
  data.frame(y = ypos, label = paste('n =', sum(!is.na(x))))
}
n_below <- function(x) data.frame(y = -0.1, label = paste('n =', sum(!is.na(x))))

harvardplot <- ggplot(macrodata %>% filter(ncii > 10, grepl('Harvard',Site), !is.na(CII), !is.na(Taxon), Taxon!=''), aes(x = factor(CII), y = Year3_RGR)) + 
  stat_summary(fun.y = 'mean', geom = 'bar') +
  stat_summary(geom = 'errorbar', width=0) +
  stat_summary(geom = 'text', fun.data = n_above, size = 2) +
  facet_wrap(~ Taxon, scales = 'free_y') +
  simpletheme +
  labs(x = 'CII') +
  ggtitle('Harvard Forest: Year 3 RGR')

# Multiple pages for Barro Colorado

ggplot(macrodata %>% filter(ncii > 10, grepl('Barro',Site), !is.na(CII), !is.na(Taxon), Taxon!=''), aes(x = factor(CII), y = Year3_RGR)) + 
  stat_summary(fun.y = 'mean', geom = 'bar') +
  stat_summary(geom = 'errorbar', width=0) +
  stat_summary(geom = 'text', fun.data = n_above) +
  facet_wrap(~ Taxon, scales = 'free_y') +
  simpletheme +
  labs(x = 'CII') +
  ggtitle('Barro Colorado Island: Year 3 RGR')

bciplotdat <- macrodata %>% filter(ncii > 10, grepl('Barro',Site), !is.na(CII), !is.na(Taxon), Taxon!='')
taxonlist <- sort(unique(bciplotdat$Taxon))

xx<-seq(1,51,by=9)
rowstouse <- cbind(xx, c((xx-1)[-1],51))

pdf('figs/rgrbarplots.pdf', height=9, width=9)

harvardplot

for (i in 1:nrow(rowstouse)) {
  bciplot_i <- ggplot(bciplotdat %>% filter(Taxon %in% taxonlist[(rowstouse[i,1]):(rowstouse[i,2])]), aes(x = factor(CII), y = Year3_RGR)) + 
    stat_summary(fun.y = 'mean', geom = 'bar') +
    stat_summary(geom = 'errorbar', width=0) +
    stat_summary(geom = 'text', fun.data = n_above, size = 2) +
    facet_wrap(~ Taxon, scales = 'free_y') +
    simpletheme +
    labs(x = 'CII') +
    ggtitle('Barro Colorado Island: Year 3 RGR')
  print(bciplot_i)
}

dev.off()

# Determine how many individuals there are in each diameter bin for each species, then split it by CII.

hist(macrodata$Year3_DBH)
macrodata <- mutate(macrodata, y3dbhclass = cut(Year3_DBH, breaks = c(0,2.5,5,10,100)))

threewaytable <- with(macrodata, table(y3dbhclass, CII, Taxon))

