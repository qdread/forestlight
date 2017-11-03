# Parameters and state variables of ITD (ideal tree distribution) model
# The ITD is a special case of the PPA (perfect plasticity approximation)

# Parameters
# ----------

# For each parameter, the first value represents shade-tolerant tree and the second value represents gap specialists.

### Parameters that should differ between shade and gap specialists

# GD: Growth rate in dbh units (cm) per timestep (y) in the dark (below the critical height z*, all trees below z* receive less light than the trees above z*)
# GL: Growth rate in dbh units (cm) per timestep (y) in the light (above the critical height z*)
# muD: Mortality rate, proportion per timestep (y) in the dark (below z*)
# muL: Mortality rate, proportion per timestep (y) in the light (above z*)

### Parameters that should (for now) be the same between shade and gap specialists

# phi: Allometric constant relating dbh to crown area. Crown area is proportional to phi * dbh^2
# alpha: Allometric constant relating dbh to height. Height is proportional to alpha * dbh^beta
# beta: Allometric constant relating dbh to height. Height is proportional to alpha * dbh^beta
# F: Fecundity rate, number of new individuals with dbh=0 per reproducing individual generated per timestep (y)

P <- data.frame(GD = c(0.1, 0.05),
                GL = c(1, 2),
                muD = c(0.06, 0.4),
                muL = c(0.01, 0.035),
                phi = c(1, 1),
                alpha = c(5, 5),
                beta = c(0.3, 0.3),
                F = c(1000,1000))
				
# State variables
# ---------------

# The state variable at each timestep is a matrix. Each guild (shade=1, gap=2) is split up into cohorts, each cohort has an abundance w and a dbh.
# The dbh is used to calculate the height and total crown area of the cohort.
# The height at which the total crown area equals the ground area is the critical height. 
# In other words, start with the tallest tree and work down, summing the crown area as you go, until you meet or exceed the ground area of the forest stand.
# That means trees of that critical height and above are filling the canopy. Everything shorter is shaded.
# Trees equal or greater than the critical height have overstory=TRUE.

### Starting values of the state variables

M1 <- data.frame(species=c(1,2),
                 cohort=c(1,1),
                 dbh=P$GL,
                 w=c(1000,1000),
                 z=with(P, alpha*GL^beta),
                 crownA=10 * pi * (1/10000) * (P$phi * P$GL) ^ 2,
                 overstory=c(TRUE,TRUE))
				 
A <- 1e6 # Plot area (a constant)
n_timesteps <- 1000 # number of years to run

# Run simulation
# --------------

# Source function
source('~/GitHub/forestlight/itd/itdfunction.r')

M <- M1 										# Initialize state variables
for (i in 1:n_timesteps) M <- ITD(P, M, A)		# Run simulation with given parameters and initial condition
# (this just takes a few seconds to run 1000 timesteps)

# Plot ratio of shade to gap across size classes
# ----------------------------------------------

n_bins <- 20

# Generate log bins for dbh
bin_edges <- 10^(seq(log10(1), log10(max(M$dbh)+0.1), length.out = n_bins + 1))
bin_midpoints <- (bin_edges[1:n_bins] + bin_edges[2:(n_bins + 1)])/2
shade_abund <- numeric(n_bins)
gap_abund <- numeric(n_bins)

for (i in 1:(length(bin_edges)-1)) {
  shade_abund[i] <- with(subset(M, species == 1 & dbh >= bin_edges[i] & dbh < bin_edges[i+1]), sum(w))
  gap_abund[i] <- with(subset(M, species == 2 & dbh >= bin_edges[i] & dbh < bin_edges[i+1]), sum(w))
}

# Data frame with midpoint of each bin and ratio of shade to gap abundance
plot_data <- data.frame(diameter = bin_midpoints, ratio = shade_abund/gap_abund)

# Draw plot
require(cowplot)

ggplot(data = plot_data,
       aes(x = diameter, y = ratio)) +
  geom_line() + geom_point() + scale_x_log10() + scale_y_log10()