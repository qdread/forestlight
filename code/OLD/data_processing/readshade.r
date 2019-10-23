# Read the shade tolerance tables and merge them
# QDR/Forestlight/11.11.16/12.14.16

library(XLConnect)
comita <- readWorksheetFromFile('C:/Users/Q/Google Drive/ForestLight/data/ComitaAppendix1.xlsx', sheet=1)
laselva <- read.csv('C:/Users/Q/Google Drive/ForestLight/data/shade_tolerance_la_selva.csv', stringsAsFactors = FALSE)

# Clean comita df
library(dplyr)

comita[comita == '-'] <- NA
comita <- comita %>% rename(maxsaplinggrowth = Max.sapling.growth.mm.y.1) %>% mutate(maxsaplinggrowth = as.numeric(maxsaplinggrowth))

laselva <- mutate(laselva, sppcode = gsub('_', ' ', sppcode))

table(comita$Species %in% laselva$sppcode)
table(laselva$sppcode %in% comita$Species) # 27 species overlap.

comita$Shade.tolerance.guild[comita$Shade.tolerance.guild == 'G'] <- 'intolerant'
comita$Shade.tolerance.guild[comita$Shade.tolerance.guild == 'I'] <- 'intermediate'
comita$Shade.tolerance.guild[comita$Shade.tolerance.guild == 'S'] <- 'tolerant'

allshadetol <- full_join(comita, laselva %>% rename(Species = sppcode), by = 'Species')

x <- filter(allshadetol, !is.na(Shade.tolerance.guild) & !is.na(shadetol)) # Quite a few of these do not match!
table(x$Shade.tolerance.guild == x$shadetol) # Only 9 of the 16 species found in both datasets were assigned to the same shade tolerance guild in both datasets.



# Merge John's data (14 Dec) ----------------------------------------------


fp <- 'C:/Users/Q/Google Drive/ForestLight/data'
growthdf <- read.csv(file.path(fp, 'Growth2.csv'), stringsAsFactors = FALSE)
lightdf <- read.csv(file.path(fp, 'Light.csv'), stringsAsFactors = FALSE)

# Fix mismatches in data frames:
# Replace spaces with underscores in the growth Site name
growthdf$Site <- gsub(' ', '_', growthdf$Site)

mergeddf <- merge(growthdf, lightdf, by=c("Site", "Site_code", "Line", "Main_Stem", "LTER_Tag", "ENQ_Tag"), all.x=T, all.y=F)


# Second attempt to merge data (16 Dec) -----------------------------------

setwd('C:/Users/Q/Google Drive/ForestLight/data/Data To Merge/')
bcidf <- read.csv('BCILight.csv', stringsAsFactors = FALSE)
growthdf <- read.csv('BiomassAndGrowthSTM2014Grady!.csv', stringsAsFactors = FALSE)
harvdf <- read.csv('HarvLight.csv', stringsAsFactors = FALSE)

# Combine bci and harvard
lightdf <- rbind(bcidf, harvdf[, names(bcidf)])

# Add NA's to enqtag for lightdf
lightdf$ENQ_Tag[lightdf$ENQ_Tag == ''] <- NA
lightdf$Main_Stem <- as.numeric(lightdf$Main_Stem)
lightdf$Site <- gsub(' ', '_', lightdf$Site)
lightdfsub <- subset(lightdf, !is.na(CII))

# Correct Site name on growthdf
growthdf$Site <- gsub(' ', '_', growthdf$Site)

growthdf$Main_Stem <- as.numeric(growthdf$Main_Stem)

# Rename columns of growthdf to match those of lightdf
names(growthdf)[names(growthdf) == 'LTER_Tag_.'] <- 'LTER_Tag'
names(growthdf)[names(growthdf) == 'ENQ_Tag_.'] <- 'ENQ_Tag'
names(growthdf)[names(growthdf) == 'Site_code'] <- 'Site_Code'

growthdf$origidx <- 1:nrow(growthdf)

mergeddf <- merge(growthdf, lightdfsub, by=c("Site", "Site_Code", "Line", "Main_Stem", "LTER_Tag", "ENQ_Tag"), all.x=T, all.y=F, sort = F)

t1 <- table(mergeddf$origidx)
problemrows <- as.numeric(names(t1[t1>1]))
growthdf[problemrows,]
mergeddf[mergeddf$origidx %in% problemrows, ]

lightdf[which(lightdf$ENQ_Tag == '5186'),]

# Remove problematic row from lightdf
lightdf <- lightdf[-6846,]

mergeddf <- merge(growthdf[,-88], lightdfsub, by=c("Site", "Site_Code", "Line", "Main_Stem", "LTER_Tag", "ENQ_Tag"), all.x=T, all.y=F, sort = F)

# Sort merged df
mergeddf <- mergeddf[with(mergeddf, order(Site, Site_Code, Line, LTER_Tag, ENQ_Tag)), ]

write.csv(mergeddf, file = 'GrowthLightMerged.csv', row.names=FALSE)


library(dplyr)

merge2 <- left_join(growthdf[,-88], lightdfsub)

tagnames <- c('Main_Stem','LTER_Tag')
mergeddf_nolter <- merge(growthdf[, !names(growthdf) %in% c('origidx',tagnames)], lightdfsub[, !names(lightdfsub) %in% tagnames], all.x=T, all.y=F, sort=F)


# 20 Dec: merge the lter and the enq tags separately

mergelter <- merge(subset(growthdf[,-88], !is.na(LTER_Tag)), lightdfsub, by = c('Site', 'Site_Code', 'Line', 'Main_Stem', 'LTER_Tag'), all.x=T, all.y=F, sort=F)
mergeenq <- merge(subset(growthdf[,-88], !is.na(ENQ_Tag)), lightdfsub, by = c('Site', 'Site_Code', 'Line', 'Main_Stem', 'ENQ_Tag'), all.x=T, all.y=F, sort=F)

mergeenq <- subset(mergeenq, is.na(LTER_Tag.x))

mergeenq <- mergeenq[,!names(mergeenq) %in% c('LTER_Tag.x', 'LTER_Tag.y')]
mergelter <- mergelter[,!names(mergelter) %in% c('ENQ_Tag.x', 'ENQ_Tag.y')]

mergeboth <- merge(mergeenq, mergelter, all=T)
mergeboth <- mergeboth[with(mergeboth, order(Site, Site_Code, Line, LTER_Tag, ENQ_Tag)), ]

# Reorder columns
mergeboth <- mergeboth[,c(1:3, 90, 89, 4:88)]
write.csv(mergeboth, file = 'GrowthLightMerged.csv', row.names=FALSE)
