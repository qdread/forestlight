curve(dlnorm(x, meanlog = 3.783, sdlog = 1.071)/dlnorm(x, meanlog = 3.002, sdlog = 0.955), from = 0.1, to = 100000, log = 'xy', lwd = 2)
curve(0.02*x^1.05, from = 0.1, to = 100000, log = 'xy', add = TRUE, lwd = 2, col = 'red')
text(c(100,3000), c(100,10), c('Ratio of two lognormals','Power law\nor ratio of two power laws'), col = c('black','red'))
