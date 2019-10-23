ITD <- function(P, M, A) {
  # 1. calculate critical canopy height z*
  
  # add the crown areas of each cohort starting with the tallest
  M <- M[order(M$z, decreasing = TRUE), ]
  totalcrownA <- cumsum(M$crownA)
  # find the height of the last cohort added before plot area was reached
  # if plot area is not reached, all cohorts are in the canopy
  zstar <- if (sum(M$crownA) >= A) M$z[which(totalcrownA >= A)[1]]
  else min(M$z)
  # create logical for canopy status
  M$overstory <- M$z >= zstar
  
  # 2. implement mortality
  
  # assign canopy or understory mortality parameter to each cohort in each species
  mu <- apply(cbind(M$overstory, P$muL[M$species], P$muD[M$species]), 1, function(x) if (x[1]) x[2] else x[3])
  # get new abundance for each cohort
  M$w <- M$w * (1 - mu)
  
  # 3. implement growth
  
  # assign canopy or understory mortality parameter to each cohort in each species
  G <- apply(cbind(M$overstory, P$GL[M$species], P$GD[M$species]), 1, function(x) if (x[1]) x[2] else x[3])
  # get new diameter for each cohort
  M$dbh <- M$dbh + G
  # get new height for each cohort with height-dbh allometry
  M$z <- P$alpha[M$species] * M$dbh ^ (P$beta[M$species])
  # get new crown area for each cohort with crown area-dbh allometry
  M$crownA <- M$w * pi * (1/10000) * (P$phi[M$species] * M$dbh) ^ 2
  # sort the list of cohorts by height
  M <- M[order(M$z, decreasing = TRUE), ]
  # add the crown areas of each cohort starting with the tallest
  totalcrownA <- cumsum(M$crownA)
  # find the height of the last cohort added before plot area was reached
  # if plot area is not reached, all cohorts are in the canopy
  zstar <- if (sum(M$crownA) >= A) M$z[which(totalcrownA >= A)[1]]
  else min(M$z)
  # create logical for canopy status
  M$overstory <- M$z >= zstar
  
  # 4. implement reproduction
  
  # calculate proportion exposed crown area for each species
  ACanopy <- with(subset(M, overstory == TRUE), sum(crownA))
  AjCanopy <- rep(0, nrow(P))
  for (j in 1:nrow(P)) {
    Aj <- with(subset(M, species == j), sum(crownA))
    AjCanopy[j] <- Aj - with(subset(M, species == j & overstory == FALSE), sum(crownA))
  }
  AjCanopy[is.na(AjCanopy)] <- 0
  # multiply proportion by fecundity to get number of individuals in next cohort
  wnew <- P$F * (AjCanopy/ACanopy)
  seedlings <- data.frame(species = 1:nrow(P),
                          cohort = max(M$cohort) + 1,
                          dbh = 0,
                          w = wnew,
                          z = 0,
                          crownA = 0,
                          overstory = FALSE)
  # add new cohorts to output matrix
  M <- rbind(M, seedlings)
}



# test:
P <- data.frame(GD = c(0.1, 0.05),
                GL = c(1, 2),
                muD = c(0.06, 0.4),
                muL = c(0.01, 0.035),
                phi = c(1, 1),
                alpha = c(5, 5),
                beta = c(0.3, 0.3),
                F = c(1000,1000))
M1 <- data.frame(species=c(1,2),
                 cohort=c(1,1),
                 dbh=P$GL,
                 w=c(1000,1000),
                 z=with(P, alpha*GL^beta),
                 crownA=10 * pi * (1/10000) * (P$phi * P$GL) ^ 2,
                 overstory=c(TRUE,TRUE))

A <-1e6

M <- M1
for (i in 1:1000) M <- ITD(P, M, A)

table(M$species)
sum(M$crownA)

# Plot the final diameter ratios.

# Find log bins for dbh
bin_edges <- 10^(seq(log10(1), log10(max(M$dbh)+0.1), length.out = 21))
bin_midpoints <- (bin_edges[1:20] + bin_edges[2:21])/2
shade_abund <- numeric(20)
gap_abund <- numeric(20)

for (i in 1:(length(bin_edges)-1)) {
  shade_abund[i] <- with(subset(M, species == 1 & dbh >= bin_edges[i] & dbh < bin_edges[i+1]), sum(w))
  gap_abund[i] <- with(subset(M, species == 2 & dbh >= bin_edges[i] & dbh < bin_edges[i+1]), sum(w))
  
}

library(cowplot)
ggplot(data = data.frame(diameter = bin_midpoints, ratio = shade_abund/gap_abund),
       aes(x = diameter, y = ratio)) +
  geom_point() + scale_x_log10() + scale_y_log10()

ggsave('C:/Users/Q/google_drive/ForestLight/figs/PPA_prediction_1.png', height=5, width=5)
