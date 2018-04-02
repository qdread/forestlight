# Fit function to individual trees' production as a function of light received per crown area

library(cowplot)

ggplot(alltree_light_95, aes(x = light_received/crownarea, y = production/crownarea)) +
  geom_hex() +
  scale_x_log10() + scale_y_log10() +
  scale_fill_gradient(low = 'gray90', high = 'gray10')


ggplot(alltree_light_95, aes(x = light_received, y = production)) +
  geom_hex() +
  scale_x_log10() + scale_y_log10() +
  scale_fill_gradient(low = 'gray90', high = 'gray10')

# Fit linear model and fit nonlinear least squares model in each case.

dat <- alltree_light_95 %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea)

lm1 <- lm(log10(production_area) ~ log10(light_area), data = dat)
nls1 <- nls(production_area ~ L/(1 + exp(-k*(light_area-x0))), 
            data = dat, 
            start = list(L = 2, x0 = 10, k = 1))
nls2 <- nls(production_area ~ beta0 * light_area ^ beta1, data = dat, start = list(beta0 = -1, beta1 = 1))
nls3 <- nls(production_area ~ (beta0 * light_area ^ beta1) * (a + light_area ^ b + c), data = dat, start = list(beta0 = -1, beta1 = 1, a = -0.1, b = 0.1, c = -0.1))

lmcoef <- coef(lm1)
nlscoef <- coef(nls1)
nls2coef <- coef(nls2)

lm_fn <- function(x, b0, b1) 10^(b0) * x^b1
nls_fn <- function(x, L, x0, k) L/(1 + exp(-k*(x - x0)))
nls2_fn <- function(x, beta0, beta1) beta0 * x^beta1

with(alltree_light_95, plot(light_received/crownarea, production/crownarea, log='xy'))
curve(lm_fn(x, b0=lmcoef[1], b1=lmcoef[2]), col = 'red', add = TRUE)
curve(nls_fn(x, L = nlscoef[1], x0 = nlscoef[2], k = nlscoef[3]), col = 'blue', add = TRUE)
curve(nls2_fn(x, beta0 = nls2coef[1], beta1 = nls2coef[2]), col = 'green', add = TRUE)


ggplot(alltree_light_95, aes(x = light_received/crownarea, y = production/crownarea)) +
  geom_hex() +
  geom_abline(slope = lmcoef[2], intercept=lmcoef[1], color = 'red') +
  scale_x_log10() + scale_y_log10() +
  scale_fill_gradient(low = 'gray90', high = 'gray10') 
  
