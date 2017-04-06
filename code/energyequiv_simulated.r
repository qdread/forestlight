# Mass & diameter energy equivalence with some simulated data

# Details of simulation
rmin <- 1
rmax <- 100
binwidth <- 0.1

# Constants
n0 <- 100000
b0 <- 0.1
k <- 1

r <- seq(from=rmin, to=rmax, by=binwidth)
n <- round(n0 * r^-2)
b <- b0 * r^2
m <- k * r^(8/3)

radii <- rep(r, n)
masses <- rep(m, n)
productions <- rep(b, n)

# Energy equivalence by radius
plot(r, n*b, ylim=c(0.1, max(n*b)), log='xy')

# Bin linearly by mass into 100 bins
masscuts <- seq(0, max(m), length.out=length(r) + 1)
masslower <- masscuts[-(length(r) + 1)]
massupper <- masscuts[-1]
massupper[length(r)] <- massupper[length(r)] + 1

mb <- rep(0, length(r))

for (i in 1:length(r)) {
  mb[i] <- which(m[i] >= masslower & m[i] < massupper)
}

mbfact <- as.numeric(factor(mb, labels=1:length(unique(mb))))

massbins <- rep(mbfact, n)
prodbinnedbymass <- tapply(productions, massbins, sum)

plot(m[unique(mb)], prodbinnedbymass, log='xy')
