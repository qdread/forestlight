#from Kitijama, 2005 Annals of Botany, using datathief 1.7
library(tidyverse)
library(broom)

lai <- read_csv('/Users/jgradym/Google_Drive/ForestLight/data/data_forplotting/LAI_Depth.csv')
ggplot(data = lai,aes( x = Depth, y = LAI, fill = Species, color = Species )) +
  geom_point(size = 6, shape = 21, color = "black") + theme_plant  +
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Species, color = Species))
lm1 <- lm(log(LAI) ~ log(Depth)+ Species, data = lai)
summary(lm1) # 0.35
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#  log(Depth)         0.353899   0.025583  13.834  < 2e-16 ***

#Residual standard error: 0.2065 on 55 degrees of freedom
#Multiple R-squared:  0.9562,	Adjusted R-squared:  0.9522 
#F-statistic: 240.2 on 5 and 55 DF,  p-value: < 2.2e-16


lm_each <-  lai %>%
  nest(-Species) %>% #group variable
  mutate(
    fit = map(data, ~ lm(log(LAI) ~ log(Depth), data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)%>%
  filter(term != '(Intercept)')
lm_each
#> lm_each
# A tibble: 5 x 8
#Species    data              fit    term       estimate std.error statistic  p.value
#<chr>      <list>            <list> <chr>         <dbl>     <dbl>     <dbl>    <dbl>
# 1 Castilla   <tibble [9 × 2]>  <lm>   log(Depth)    0.420    0.0275    15.2   1.26e- 6
#2 Luehea     <tibble [11 × 2]> <lm>   log(Depth)    0.323    0.0218    14.8   1.26e- 7
#3 Anacardium <tibble [19 × 2]> <lm>   log(Depth)    0.389    0.0248    15.7   1.54e-11
#4 Antirrhoea <tibble [13 × 2]> <lm>   log(Depth)    0.357    0.0171    20.9   3.33e-10
#5 Cecropia   <tibble [9 × 2]>  <lm>   log(Depth)    0.202    0.226      0.891 4.02e- 1

growth_diam_slopes_under_10
  lm(log(LAI) ~ log(Depth)+ Species, data = lai)
summary(lm1) # 0.35

trans <- data.frame("LAI" = c(0.44, 3.8, 8), "Trans" = c(60, 20.5, 5.4))
lm2 <- lm(log(Trans) ~ log(LAI), data = trans)
summary(lm2) #-0.75
