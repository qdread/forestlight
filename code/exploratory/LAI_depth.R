#from Kitijama, 2005 Annals of Botany, using datathief 1.7 to digitize data
library(tidyverse)
library(broom)
library(lme4)

lai <- read_csv('/Users/jgradym/Google_Drive/ForestLight/data/data_forplotting/LAI_Depth.csv')
ggplot(data = lai,aes( x = Depth, y = LAI, fill = Species, color = Species )) +
  geom_point(size = 6, shape = 21, color = "black") + theme_plant  +
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Species, color = Species))

# with species
lm1 <- lm(log(LAI) ~ log(Depth)+ Species, data = lai)
summary(lm1) # 0.351, r2 = 0.99

# without species
lm2 <- lm(log(LAI) ~ log(Depth)+ Species, data = lai)
summary(lm2) # 0.54, r2 = 0.38

# calculated separately
lm3 <- lm(log(LAI) ~ log(Depth), data = lai2)
summary(lm3)
lm_each <-  lai %>%
  nest(-Species) %>% #group variable
  mutate(
    fit = map(data, ~ lm(log(LAI) ~ log(Depth), data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)%>%
  filter(term != '(Intercept)')
lm_each
mean(lm_each$estimate) #0.331


# Mixed Model
lai2 <- lai %>% filter(Species != "Cecropia") #Cecropia is weird

# all 5 species
mlm <- lmer(log(LAI) ~ log(Depth) + (log(Depth) |Species), lai) 
mlm # 0.33 sloep

#no Cecropia 
mlm2 <- lmer(log(LAI) ~ log(Depth) + (log(Depth) |Species), lai2) 
mlm2 #0.37 slope

