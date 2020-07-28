# Diag. plot

# fitted line fn
twoseg_log <- function(x, alpha, beta, tau) {
  x2 <- ifelse(log10(x) < tau, 0, 1)
  10 ^ (alpha + beta[1] * log10(x) + beta[2] * (log10(x) - tau) * x2)
}

threeseg_log <- function(x, alpha, beta, tau_low, tau_high) {
  x2a <- ifelse(log10(x) < tau_low, 0, 1)
  x2b <- ifelse(log10(x) < tau_high, 0, 1)
  10 ^ (alpha + beta[1] * log10(x) + beta[2] * (log10(x) - tau_low) * x2a + beta[3] * (log10(x) - tau_high) * x2b)
}

xpred <- logseq(1, 350, 50)

ypred <- twoseg_log(xpred, alpha = 2.95, beta = c(-1.2, -2.19), tau = 1.73)
ypred3seg <- threeseg_log(xpred, alpha = 2.85, beta = c(-0.96, -0.58, -2.13), tau_low = 0.94, tau_high = 1.88)

plot(dat_all$x, dat_all$y, log = 'xy')
points(xpred, ypred, type = 'l', col = 'blue')
points(xpred, ypred3seg, type = 'l', col = 'red')

# Light

xpredlight <- logseq(1, 400, 50)

ypred2seglight <- twoseg_log(xpredlight, alpha = 1.6, beta = c(1.03, -2.12), tau = 0.63)
ypred3seglight <- threeseg_log(xpredlight, alpha = 1.47, beta = c(1.64, -2.42, -0.49), tau_low = 0.47, tau_high = 1.49)

plot(dat_all_light$x, dat_all_light$y, log = 'xy', ylim = c(1, 150))
points(xpredlight, ypred2seglight, type = 'l', col = 'blue')
points(xpredlight, ypred3seglight, type = 'l', col = 'red')

# Mixed model: diam. Extract coefficients
coef_2segmixed <- summary(fit_2seg_mixed_all)[[1]][,'mean']
# Reshape coefs.
coef_2segmixed_df <- data.frame(fg = 1:5, 
                                alpha = coef_2segmixed[grep('coef_alpha', names(coef_2segmixed))],
                                beta_low = coef_2segmixed[grep('coef_beta_low', names(coef_2segmixed))],
                                beta_high = coef_2segmixed[grep('coef_beta_high', names(coef_2segmixed))],
                                tau = coef_2segmixed[grep('coef_tau', names(coef_2segmixed))]
                                )

ypred2segmixed <- coef_2segmixed_df %>%
  group_by(fg) %>%
  group_modify(~ data.frame(x = xpred, ypred = twoseg_log(x = xpred, alpha = .$alpha, beta = c(.$beta_low, .$beta_high), tau = .$tau)))

ggplot() +
  geom_point(aes(x = bin_midpoint, y = richness_by_bin_width, color = fg), data = bin_x_fg_use) +
  geom_line(aes(x = x, y = ypred, color = fg), data = ypred2segmixed %>% ungroup %>% mutate(fg = paste0('fg',fg))) +
  theme_bw() +
  scale_x_log10(name = 'Diameter') + scale_y_log10('Richness per bin')

# Mixed model 3 segment: diam. Extract coefficients
coef_3segmixed <- summary(fit_3seg_mixed_all)[[1]][,'mean']
# Reshape coefs.
coef_3segmixed_df <- data.frame(fg = 1:5, 
                                alpha = coef_3segmixed[grep('coef_alpha', names(coef_3segmixed))],
                                beta_low = coef_3segmixed[grep('coef_beta_low', names(coef_3segmixed))],
                                beta_mid = coef_3segmixed[grep('coef_beta_mid', names(coef_3segmixed))],
                                beta_high = coef_3segmixed[grep('coef_beta_high', names(coef_3segmixed))],
                                tau_low = coef_3segmixed[grep('coef_tau_low', names(coef_3segmixed))],
                                tau_high = coef_3segmixed[grep('coef_tau_high', names(coef_3segmixed))]
)

ypred3segmixed <- coef_3segmixed_df %>%
  group_by(fg) %>%
  group_modify(~ data.frame(x = xpred, ypred = threeseg_log(x = xpred, alpha = .$alpha, beta = c(.$beta_low, .$beta_mid, .$beta_high), tau_low = .$tau_low, tau_high = .$tau_high)))

ggplot() +
  geom_point(aes(x = bin_midpoint, y = richness_by_bin_width, color = fg), data = bin_x_fg_use) +
  geom_line(aes(x = x, y = ypred, color = fg), data = ypred3segmixed %>% ungroup %>% mutate(fg = paste0('fg',fg)) %>% filter(ypred > 1e-2)) +
  theme_bw() +
  scale_x_log10(name = 'Diameter') + scale_y_log10('Richness per bin')
