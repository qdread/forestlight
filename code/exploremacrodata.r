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

ggplot(macrodata %>% filter(ncii > 10, grepl('Harvard',Site), !is.na(CII), !is.na(Taxon), Taxon!=''), aes(x = factor(CII), y = Year3_RGR)) + 
  stat_summary(fun.y = 'mean', geom = 'bar') +
  stat_summary(geom = 'errorbar', width=0) +
  stat_summary(geom = 'text', )
  facet_wrap(~ Taxon, scales = 'free_y') +
  simpletheme +
  ggtitle('Harvard Forest')

ggplot(macrodata %>% filter(ncii > 10, grepl('Barro',Site), !is.na(CII), !is.na(Taxon), Taxon!=''), aes(x = factor(CII), y = Year3_RGR)) + 
  stat_summary(fun.y = 'mean', geom = 'bar') +
  stat_summary(geom = 'errorbar', width=0) +
  facet_wrap(~ Taxon, scales = 'free_y') +
  simpletheme +
  ggtitle('Barro Colorado Island')
