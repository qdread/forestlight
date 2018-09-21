# Combine predicted bin csvs into one

pbins <- lapply(1:140, function(i) read.csv(paste0('~/forestlight/stanoutput/fitinfo/predbin_',i,'.csv'), stringsAsFactors = FALSE))
pbins <- do.call(rbind, pbins)
write.csv(pbins[,-1], '~/forestlight/stanoutput/fitinfo/predbin_all.csv', row.names = FALSE)