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
