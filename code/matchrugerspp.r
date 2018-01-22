# Match taxa in Nadja's dataset with Condit's dataset
all_condit_spp <- unique(bci.stem7$sp)
fgbci <- mutate(fgbci, sp = tolower(sp))

taxmatch <- fgbci$sp %in% all_condit_spp 
fgbci[!taxmatch, ] # Two don't match.
table(all_condit_spp %in% fgbci$sp) # 280 are in.