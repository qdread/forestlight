# Load FIA data to calculate growth and mortality.
# Uses Massachusetts data

fp <- 'C:/Users/Q/Dropbox/projects/forestlight/'

ma <- read.csv(file.path(fp, 'MA/MA_TREE.csv'), stringsAsFactors = FALSE)
nh <- read.csv(file.path(fp, 'NH/NH_TREE.csv'), stringsAsFactors = FALSE)
vt <- read.csv(file.path(fp, 'VT/VT_TREE.csv'), stringsAsFactors = FALSE)
me <- read.csv(file.path(fp, 'ME/ME_TREE.csv'), stringsAsFactors = FALSE)
ct <- read.csv(file.path(fp, 'CT/CT_TREE.csv'), stringsAsFactors = FALSE)
ri <- read.csv(file.path(fp, 'RI/RI_TREE.csv'), stringsAsFactors = FALSE)
ny <- read.csv(file.path(fp, 'NY/NY_TREE.csv'), stringsAsFactors = FALSE)

fiaall <- rbind(ma, nh, vt, me, ct, ri, ny)

# Summarize number of measurements
library(dplyr)

nmeas <- fiaall %>% group_by(STATECD, PLOT, SUBP, TREE) %>% summarize(ndia = sum(!is.na(DIA)))

calcrgr <- function(x) {
  if (sum(!is.na(x$DIA)) > 1) {
    xsub <- subset(x, !is.na(DIA))
    gr <- diff(xsub$DIA)/diff(xsub$INVYR)
    rgr <- gr/(xsub$DIA[-1])
    return(rgr)
  }
  else {
    return(NA)
  }
  
}

rgrs <- fiaall %>% group_by(STATECD, PLOT, SUBP, TREE) %>%
  do(rgr = calcrgr(.))