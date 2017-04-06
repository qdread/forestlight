# Energy equivalence plots using simulated data

# Data are simulated to achieve perfect energy equivalence by diameter size class

# abundance equation: n = n0 * r^-2
# production equation: b = b0 * r^2
# standing biomass equation: m = m0 * r^(8/3)

# production by diameter equation: pdiam = n0 * b0 * r^0 = n0 * b0
# production by mass equation: pmass = 

# Analytical solution for finite bin width

pbyradius <- function(n0, b0, rl, ru) (1/3) * n0 * b0 * (rl^(-1) - ru^(-1)) * (ru^3 - rl^3)
pbymass <- function(n0, b0, m0, ml, mu) (1/3) * n0 * b0 * m0^2 * (ml^(-3/8) - mu^(-3/8)) * (mu^(9/8) - ml^(9/8))

pbyradius_fixedbin <- function(n0, b0, r, deltar) (1/3) * n0 * b0 * (r^(-1) - (r+deltar)^(-1)) * ((r+deltar)^3 - r^3)
pbymass_fixedbin <- function(n0, b0, m0, r, deltar) (1/3) * n0 * b0 * (mass(m0, r)^(-1) - mass(m0, r+deltar)^(-1)) * (mass(m0, r+deltar)^3 - mass(m0, r)^3)

mass <- function(m0, r) m0^(3/8) * r^(3/8)  
  
n0 <- 1
b0 <- 1
m0 <- 1

binwidth <- 1

radii <- seq(1, 100, by = binwidth)
masses <- m0 * radii^(3/8)

radius_bin_lower <- radii[-length(radii)]
radius_bin_upper <- radii[-1]
radius_bin_midpoint <- (radius_bin_lower + radius_bin_upper)/2
mass_bin_lower <- masses[-length(masses)]
mass_bin_upper <- masses[-1]
mass_bin_midpoint <- (mass_bin_lower + mass_bin_upper)/2

p_r <- pbyradius(n0, b0, radius_bin_lower, radius_bin_upper)
p_m <- pbymass(n0, b0, m0, mass_bin_lower, mass_bin_upper)

plot(radius_bin_midpoint, p_r, log = 'xy', ylim = c(0.1, max(p_r)))
plot(mass_bin_midpoint, p_m, log = 'xy', ylim = c(1e-7, max(p_m)))

pdf('C:/Users/Q/Google Drive/ForestLight/figs/ee_solution_finitebin.pdf', height=5, width=5)
curve(pbyradius_fixedbin(n0=1,b0=1,r=x,deltar=1), from=1, to=100, log='xy', ylim=c(0.1, 1.1), xlab = 'Radius', ylab = 'Radius-binned production')
curve(pbymass_fixedbin(n0=1,b0=1,m0=1,r=x,deltar=1), from=1, to=100, log='xy', ylim=c(1e-5, 0.1), xlab = 'Radius', ylab = 'Mass-binned production')
dev.off()

# Better plot with ggplot, showing the effect of declining deltar.
library(ggplot2)
pdf('C:/Users/Q/Google Drive/ForestLight/figs/ee_solution_finitebin.pdf', height=5, width=5)
ggplot(data.frame(r = c(1,100)), aes(x=r)) +
  stat_function(fun = pbyradius_fixedbin, args = list(n0=1,b0=1,deltar=1), aes(colour='1'), xlim=c(0.1,10)) +
  stat_function(fun = pbyradius_fixedbin, args = list(n0=1,b0=1,deltar=0.5), aes(colour='0.5'), xlim=c(0.1,10)) +
  stat_function(fun = pbyradius_fixedbin, args = list(n0=1,b0=1,deltar=0.1), aes(colour='0.1'), xlim=c(0.1,10)) +
  scale_y_log10(name = 'Size-class production', limits = c(0.001, 1.3)) + scale_x_log10(name = 'Radius') +
  scale_color_manual(name=expression(Delta*r), values = c('red','blue','green')) +
  theme_bw() +
  ggtitle('Binned by radius')

ggplot(data.frame(r = c(1,100)), aes(x=r)) +
  stat_function(fun = pbymass_fixedbin, args = list(n0=1,b0=1,m0=1,deltar=1), aes(colour='1'), xlim=c(0.1,10)) +
  stat_function(fun = pbymass_fixedbin, args = list(n0=1,b0=1,m0=1,deltar=0.5), aes(colour='0.5'), xlim=c(0.1,10)) +
  stat_function(fun = pbymass_fixedbin, args = list(n0=1,b0=1,m0=1,deltar=0.1), aes(colour='0.1'), xlim=c(0.1,10)) +
  scale_y_log10(name = 'Size-class production') + scale_x_log10(name = 'Mass') +
  scale_color_manual(name=expression(Delta*r), values = c('red','blue','green')) +
  theme_bw() +
  ggtitle('Binned by mass')
dev.off()


#######################
# Corrected

pbyradius_fixedbin <- function(n0, b0, r, deltar) (1/3) * n0 * b0 * (r^(-1) - (r+deltar)^(-1)) * ((r+deltar)^3 - r^3)
pbymass_fixedbin <- function(n0, b0, m, deltam) (16/7) * n0 * b0 * ((m+deltam)^(1/4) - m^(1/4)) * ((m+deltam)^(7/4) - m^(7/4))

ggplot(data.frame(m = c(1,100)), aes(x=m)) +
  stat_function(fun = pbymass_fixedbin, args = list(n0=1,b0=1,deltam=1), aes(colour='1'), xlim=c(0.1,10)) +
  stat_function(fun = pbymass_fixedbin, args = list(n0=1,b0=1,deltam=0.5), aes(colour='0.5'), xlim=c(0.1,10)) +
  stat_function(fun = pbymass_fixedbin, args = list(n0=1,b0=1,deltam=0.1), aes(colour='0.1'), xlim=c(0.1,10)) +
  scale_y_log10(name = 'Size-class production') + scale_x_log10(name = 'Mass') +
  scale_color_manual(name=expression(Delta*m), values = c('red','blue','green')) +
  theme_bw() +
  ggtitle('Binned by mass')

########################
# Bin widths
bw <- function(k, r, deltar) k^(-8/3) * ((r + deltar)^(8/3) - r^(8/3))
bw(1, 1:10, 1)


########################
# Added 06 April: Plot with different bin widths
# Assume energy equivalence by radius

pm <- function(n0, b0, k, m, deltam) n0 * b0 * k * ( (m + deltam)^(3/8) - m^(3/8) )

library(ggplot2)
p <- ggplot(data.frame(m = c(1, 100000)), aes(x=m)) +
  scale_x_log10(name = 'Mass') + 
  theme_bw() +
  stat_function(geom = 'line', fun = pm, args = list(n0 = 100000, b0 = 0.1, k = 1, deltam = 1), aes(colour = '1')) +
  stat_function(geom = 'line', fun = pm, args = list(n0 = 100000, b0 = 0.1, k = 1, deltam = 0.1), aes(colour = '0.1')) +
  stat_function(geom = 'line', fun = pm, args = list(n0 = 100000, b0 = 0.1, k = 1, deltam = 0.01), aes(colour = '0.01')) +
  stat_function(geom = 'line', fun = pm, args = list(n0 = 100000, b0 = 0.1, k = 1, deltam = 0.001), aes(colour = '0.001')) +
  stat_function(geom = 'line', fun = pm, args = list(n0 = 100000, b0 = 0.1, k = 1, deltam = 0.0001), aes(colour = '0.0001')) +
  scale_color_manual(name=expression(Delta*m), values = rainbow(5)) 

p + scale_y_log10(name = 'Total production')
p + scale_y_continuous(name = 'Total production')


############################
# Further added on 06 April. Plot density scaling relationship with finite bin width.

nm <- function(n0, k, m, deltam) 4 * n0 * k * ( (m+deltam)^(1/4) - m^(1/4) )

p <- ggplot(data.frame(m = c(1, 100000)), aes(x=m)) +
  scale_x_log10(name = 'Mass') + 
  theme_bw() +
  stat_function(geom = 'line', fun = nm, args = list(n0 = 100000, k = 1, deltam = 1), aes(colour = '1')) +
  stat_function(geom = 'line', fun = nm, args = list(n0 = 100000, k = 1, deltam = 0.1), aes(colour = '0.1')) +
  stat_function(geom = 'line', fun = nm, args = list(n0 = 100000, k = 1, deltam = 0.01), aes(colour = '0.01')) +
  stat_function(geom = 'line', fun = nm, args = list(n0 = 100000, k = 1, deltam = 0.001), aes(colour = '0.001')) +
  stat_function(geom = 'line', fun = nm, args = list(n0 = 100000, k = 1, deltam = 0.0001), aes(colour = '0.0001')) +
  scale_color_manual(name=expression(Delta*m), values = rainbow(5)) 

p + scale_y_log10(name = 'Density')
p + scale_y_continuous(name = 'Density')
