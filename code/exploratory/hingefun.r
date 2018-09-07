f2piece <- function(x, a, b, c, b1, tau) exp ( log(-a*x^(-b) + c) + b1 * (log(x) - log(tau)) * (x >= tau) )

curve(f2piece(x, 3, 0.5, 5, 2, 3), from = 1, to = 100, log = 'xy')

# See
# https://andrewgelman.com/2017/05/19/continuous-hinge-function-bayesian-modeling/
logistic_hinge <- function(x, x0, a, b0, b1, delta) {
  return(a + b0 * (x - x0) + (b1 - b0) * delta * log(1 + exp((x - x0) / delta)))
}

logistic_hinge_loglog <- function(x, x0, a, b0, b1, delta) {
  x <- log(x)
  return(exp( a + b0 * (x - x0) + (b1 - b0) * delta * log(1 + exp((x - x0) / delta))) )
}


curve(logistic_hinge(x, 5, 5, 3, 1, 0.5), from = 0.1, to = 20)

x_pred <- exp(seq(log(1.1), log(300), length.out=50))

logy_pred <- logistic_hinge(log(x_pred), 1, 1, 3, 1, 0.1)
exp(logy_pred)
plot(x_pred, exp(logy_pred), type = 'l', log = 'xy')

curve(logistic_hinge_loglog(x, hinge_coefs[1], hinge_coefs[2], hinge_coefs[3], hinge_coefs[4], hinge_coefs[5]), from = 1.1, to = 300, log = 'xy', lwd = 1.5)
points(x, y, col = adjustcolor('blue', alpha.f = 0.2))
