library(ggplot2)
ggplot(data.frame(x = c(0,100), y = c(0, 10)), aes(x=x)) + theme_bw() +
  scale_x_log10() + scale_y_log10() +
  stat_function(fun = dweibull, n = 1001, args = list(shape = 0.3, scale = 0.1), color = 'blue') +
  stat_function(fun = dweibull, n = 1001, args = list(shape = 0.4, scale = 0.1)) +
  stat_function(fun = dweibull, n = 1001, args = list(shape = 0.5, scale = 0.1), color = 'red')
