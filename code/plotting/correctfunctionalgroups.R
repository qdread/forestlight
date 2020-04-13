# Correct the functional groups data so that there is no confusing correction code in the final version of the workflow that goes to publication.

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users/jgradym/Google_Drive/ForestLight'))

fgbci <- read.table(file.path(gdrive_path, 'data/Ruger/fgroups_dynamics_new.txt'), stringsAsFactors = FALSE)
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5)) # correct in advance?
fgbci$PC_slow_to_fast <- -fgbci$X1new
fgbci$PC_breeder_to_pioneer <- fgbci$X2new

readr::write_csv(fgbci, file.path(gdrive_path, 'data/functionalgroups.csv'))