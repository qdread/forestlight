# Load FIA data to calculate growth and mortality.
# Uses Massachusetts data

fp <- 'C:/Users/Q/Dropbox/projects/forestlight/'

ma <- read.csv(file.path(fp, 'MA/MA_TREE.csv'), stringsAsFactors = FALSE)

# Summarize number of measurements
library(dplyr)

nmeas <- ma %>% group_by(PLOT, SUBP, TREE) %>% summarize(ndia = sum(!is.na(DIA)))

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

rgrs <- ma %>% group_by(PLOT, SUBP, TREE) %>%
  do(rgr = calcrgr(.))