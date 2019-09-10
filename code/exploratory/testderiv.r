# can we get derivatives of this function since it is piecewise

pdf_3part <- function(x, xmin, alpha_low, alpha_mid, alpha_high, tau_low, tau_high) {
  C_con_low <- tau_low ^ -(alpha_mid + alpha_low)
  C_con_high <- tau_high ^ (alpha_high - alpha_mid)
  C_norm <- ( (C_con_low / alpha_low) * (tau_low ^ alpha_low - xmin ^ alpha_low) + (1 / alpha_mid) * (tau_low ^ -alpha_mid - tau_high ^ -alpha_mid) + (C_con_high / alpha_high) * (tau_high ^ -alpha_high) ) ^ -1
  
  prob <- case_when(
    x < tau_low ~ C_con_low * C_norm * ( x ^ (alpha_low - 1) ),
    x >= tau_low & x <= tau_high ~ C_norm * ( x ^ - (alpha_mid + 1) ),
    x > tau_high ~ C_con_high * C_norm * ( x ^ - (alpha_high + 1) )
  )
  return(prob)
}

xs <- exp(seq(log(1.2),log(315),length.out = 101))

ys <- pdf_3part(x=xs, xmin=1.1, alpha_low=0.85, alpha_mid=1.26, alpha_high=3.98, tau_low=2.4, tau_high=39.2)

ggplot(data.frame(x=xs,y=ys), aes(x=x,y=y)) + geom_line() + scale_x_log10() + scale_y_log10()

library(pracma)

dx1 <- fderiv(f = pdf_3part, x = xs, xmin=1.1, alpha_low=0.85, alpha_mid=1.26, alpha_high=3.98, tau_low=2.4, tau_high=39.2)
dx2 <- fderiv(f = pdf_3part, x = xs, n = 2, xmin=1.1, alpha_low=0.85, alpha_mid=1.26, alpha_high=3.98, tau_low=2.4, tau_high=39.2)

ggplot(data.frame(x=xs,y=ys), aes(x=x,y=y)) + geom_line() + scale_x_log10() 
ggplot(data.frame(x=xs,y=dx1), aes(x=x,y=y)) + geom_line() + scale_x_log10()
ggplot(data.frame(x=xs,y=dx2), aes(x=x,y=y)) + geom_line() + scale_x_log10()

