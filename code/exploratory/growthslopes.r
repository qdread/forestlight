# Test log log regression

library(tidyverse)

ga_vs_la <- alltree_light_95 %>% 
  mutate(GA = production/crownarea, LA = light_received/crownarea) %>%
  group_by(fg) %>%
  do(lm = lm(log(GA)~log(LA), data = .)) %>%
  mutate(slope = coef(lm)[2])

gv_vs_lv <- alltree_light_95 %>% 
  mutate(GV = production/crownvolume, LV = light_received/crownvolume) %>%
  group_by(fg) %>%
  do(lm = lm(log(GV)~log(LV), data = .)) %>%
  mutate(slope = coef(lm)[2])

dat <- alltree_light_95 %>%
  transmute(G = production, A = crownarea, V = crownvolume, L = light_received, GA = G/A, LA = L/A, GV = G/V, LV = L/V)

set.seed(333)
dat2 <- dat[sample(nrow(dat), 1000),]

ggplot(dat2, aes(x = LA, y = GA)) + geom_point() + scale_x_log10() + scale_y_log10()
ggplot(dat2, aes(x = LV, y = GV)) + geom_point() + scale_x_log10() + scale_y_log10()

apply(dat2, 2, range)

ggplot(dat) +
  geom_point(aes(x = LA, y = GA)) +
  geom_point(aes(x = LV, y = GV), color = 'red') +
  scale_x_log10() + scale_y_log10()
