# This script gives an example of the procedure to construct a prediction interval
# for a linear regression model using a bootstrap method.  The method is the one
# described in Section 6.3.3 of Davidson and Hinckley (1997),
# _Bootstrap Methods and Their Application_.


#rm(list=ls())
set.seed(12344321)
library(MASS)
library(Hmisc)

# Generate bivariate regression data
x <- runif(n=100,min=0,max=100)
y <- 1 + x + (rexp(n=100,rate=0.25)-4)

my.reg <- lm(y~x)
summary(my.reg)

# Predict y for x=78:
y.p <- coef(my.reg)["(Intercept)"] + coef(my.reg)["x"]*78
y.p

# Create adjusted residuals
leverage <- influence(my.reg)$hat
my.s.resid <- residuals(my.reg)/sqrt(1-leverage)
my.s.resid <- my.s.resid - mean(my.s.resid)


reg <- my.reg
s <- my.s.resid

the.replication <- function(reg,s,x_Np1=0){
  # Make bootstrap residuals
  ep.star <- sample(s,size=length(reg$residuals),replace=TRUE)
  
  # Make bootstrap Y
  y.star <- fitted(reg)+ep.star
  
  # Do bootstrap regression
  x <- model.frame(reg)[,2]
  bs.reg <- lm(y.star~x)
  
  # Create bootstrapped adjusted residuals
  bs.lev <- influence(bs.reg)$hat
  bs.s   <- residuals(bs.reg)/sqrt(1-bs.lev)
  bs.s   <- bs.s - mean(bs.s)
  
  # Calculate draw on prediction error
  xb.xb <- coef(my.reg)["(Intercept)"] - coef(bs.reg)["(Intercept)"] 
  xb.xb <- xb.xb + (coef(my.reg)["x"] - coef(bs.reg)["x"])*x_Np1
  return(unname(xb.xb + sample(bs.s,size=1)))
}

# Do bootstrap with 10,000 replications
ep.draws <- replicate(n=10000,the.replication(reg=my.reg,s=my.s.resid,x_Np1=78))

# Create prediction interval
y.p+quantile(ep.draws,probs=c(0.05,0.95))

# prediction interval using normal assumption
predict(my.reg,newdata=data.frame(x=78),interval="prediction",level=0.90)


# Quick and dirty Monte Carlo to see which prediction interval is better
# That is, what are the 5th and 95th percentiles of Y_{N+1}
# 
# To do it properly, I guess we would want to do the whole procedure above
# 10,000 times and then see what percentage of the time each prediction 
# interval covered Y_{N+1}

y.np1 <- 1 + 78 + (rexp(n=10000,rate=0.25)-4)
quantile(y.np1,probs=c(0.05,0.95))