library(tidyverse)
D <- seq(1,100, 1)
A_fast <- 1000*D^-1.8
A_slow <- 10000*D^-2.2 
ratio <- A_fast/A_slow 
data <- as_tibble(cbind(D, A_fast, A_slow, ratio))
data
ggplot(data = data, aes(x = D, y = ratio)) +
  geom_point(shape =20)+
  scale_y_log10() +
  scale_x_log10() + theme_classic()
#slope should be 0.4
lm1 <- lm(log(ratio) ~ log(D), data = data)
summary(lm1) #bingo
#log(D)       4.000e-01