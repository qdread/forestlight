# "deprecated" script where normalization totals are calculate
# no longer needed with better correction factor

# Production totals

alltreeprod <- map2_dfr(alltreedat, years, ~ data.frame(year = .y, fg = 'alltree', production = sum(.x$production)))
fgprod <- map2_dfr(alltreedat, years, function(x, y) {
  x <- x %>%
    mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg',fg))) %>%
    group_by(fg) %>%
    summarize(production = sum(production))
  data.frame(year = y, x)	
})

prodtotals <- rbind(alltreeprod, fgprod) %>%
  arrange(year, fg)

write_csv(prodtotals, file.path(fpdump, 'production_total.csv'))

# Total production for all trees that have light values
alltreeprod <- sum(alltree_light_95$production)
fgprod <- alltree_light_95 %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg',fg))) %>%
  group_by(fg) %>%
  summarize(production = sum(production))

prodtotals <- data.frame(year = 1995, 
                         rbind(c(fg = 'alltree', production = alltreeprod), fgprod))

write_csv(prodtotals, file.path(fpdump, 'production_total_lighttrees.csv'))

# Light received totals
alltreelightrec <- sum(dat95$light_received)
fglightrec <- dat95 %>%
  group_by(fg) %>%
  summarize(light_received = sum(light_received))

lightrectotals <- data.frame(year = 1995, 
                             rbind(c(fg = 'alltree', light_received = alltreelightrec), fglightrec))

write_csv(lightrectotals, file.path(fpdump, 'lightrec_total.csv'))

# Crown volume totals
alltreecrownvol <- sum(dat95$crownvolume)
fgcrownvol <- dat95 %>%
  group_by(fg) %>%
  summarize(crownvolume = sum(crownvolume))

crownvolumetotals <- data.frame(year = 1995, 
                                rbind(c(fg = 'alltree', crownvolume = alltreecrownvol), fgcrownvol))

write_csv(crownvolumetotals, file.path(fpdump, 'crownvol_total.csv'))
