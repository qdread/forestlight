# Fit function to individual trees' production as a function of light received per crown area

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'
load(file.path(fpdata, 'rawdataobj_22jan.r'))

library(dplyr)
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
nls3 <- nls(production_area ~ (b1 * light_area ^ a) / (1 + b2 * exp(k * light_area)),
            data = dat,
            start = list(a = -1, b1 = 0.1, b2 = 0.1, k = -1))

lmcoef <- coef(lm1)
nlscoef <- coef(nls1)
nls2coef <- coef(nls2)

lm_fn <- function(x, b0, b1) 10^(b0) * x^b1
nls_fn <- function(x, L, x0, k) L/(1 + exp(-k*(x - x0)))
nls2_fn <- function(x, beta0, beta1) beta0 * x^beta1

source('code/allfunctions27july.r')

fakebin_across_years <- function(dat_values, dat_classes, edges, mean_type = 'geometric', n_census = 5) {
  qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  # add some padding just in case
  mins <- edges$bin_min
  mins[1] <- 0
  maxes <- edges$bin_max
  maxes[length(maxes)] <- Inf
  
  binstats <- t(sapply(1:length(mins), function(x) {
    indivs <- dat_values[dat_classes >= mins[x] & dat_classes < maxes[x]]
    if (mean_type == 'geometric') {
      mean_n <- exp(mean(log(indivs)))
      sd_n <- exp(sd(log(indivs)))
      ci_width <- qnorm(0.975) * sd(log(indivs)) / sqrt(length(indivs))
      ci_min <- exp(mean(log(indivs)) - ci_width)
      ci_max <- exp(mean(log(indivs)) + ci_width)
      
    } else {
      mean_n <- mean(indivs)
      sd_n <- sd(indivs)
      ci_width <- qnorm(0.975) * sd(indivs) / sqrt(length(indivs))
      ci_min <- mean_n - ci_width
      ci_max <- mean_n + ci_width
    }
    c(mean_n_individuals = length(indivs) / n_census,
      mean = mean_n, 
      sd = sd_n,
      quantile(indivs, probs = qprobs),
      ci_min = ci_min,
      ci_max = ci_max)
  }))
  dimnames(binstats)[[2]] <- c('mean_n_individuals','mean', 'sd', 'q025', 'q25', 'median', 'q75', 'q975', 'ci_min', 'ci_max')
  data.frame(bin_midpoint = edges$bin_midpoint,
             bin_min = edges$bin_min,
             bin_max = edges$bin_max,
             binstats)
}

numbins <- 20
light_bins <- logbin(x = dat$light_area, n = numbins)

prod_light_bin <- fakebin_across_years(dat_values = dat$production_area,
                                       dat_classes = dat$light_area,
                                       edges = light_bins,
                                       n_census = 1)

with(alltree_light_95, plot(light_received/crownarea, production/crownarea, log='xy', type = 'n'))
for (i in 1:numbins) {
  lines(x = rep(prod_light_bin$bin_midpoint[i], 2), y = prod_light_bin[i, c('q025','q975')])
  lines(x = rep(prod_light_bin$bin_midpoint[i], 2), y = prod_light_bin[i, c('q25','q75')], lwd = 3)
}
points(x = prod_light_bin$bin_midpoint, y = prod_light_bin$median, cex = 1.5)
curve(lm_fn(x, b0=lmcoef[1], b1=lmcoef[2]), col = 'red', add = TRUE)
curve(nls_fn(x, L = nlscoef[1], x0 = nlscoef[2], k = nlscoef[3]), col = 'blue', add = TRUE)
curve(nls2_fn(x, beta0 = nls2coef[1], beta1 = nls2coef[2]), col = 'green', add = TRUE)


ggplot(alltree_light_95, aes(x = light_received/crownarea, y = production/crownarea)) +
  geom_hex() +
  geom_abline(slope = lmcoef[2], intercept=lmcoef[1], color = 'red') +
  scale_x_log10() + scale_y_log10() +
  scale_fill_gradient(low = 'gray90', high = 'gray10') 
  


# Bayesian way ------------------------------------------------------------

library(purrr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3) 


# Fit logistic curve with Stan.
mod_logistic <- stan_model(file='stan/logistic_curve.stan')
mod_logistic2 <- stan_model(file='stan/logistic2.stan')
mod_logistic3 <- stan_model(file='stan/logistic3.stan')

n_sub <- 5000
set.seed(919)

sample_rows <- sample(nrow(alltree_light_95), n_sub)

dat_logistic <- with(alltree_light_95[sample_rows, ], 
                     list(N = n_sub, x = light_received/crownarea, y = production/crownarea))

dat_logistic_full <- with(alltree_light_95, 
                     list(N = nrow(alltree_light_95), x = light_received/crownarea, y = production/crownarea))

NC <- 3
NI <- 6000
NW <- 5000

fit_logistic <- sampling(mod_logistic, data = dat_logistic, chains = NC, iter = NI, warmup = NW, seed = 20202)
fit_logistic2 <- sampling(mod_logistic2, data = dat_logistic, chains = NC, iter = NI, warmup = NW, seed = 112)
fit_logistic3 <- sampling(mod_logistic3, data = dat_logistic, chains = NC, iter = NI, warmup = NW, seed = 113)

library(bayesplot)
mcmc_trace(as.array(fit_logistic), pars = c('L','x0','k'))
mcmc_trace(as.array(fit_logistic2), pars = c('a', 'b1', 'b2', 'k'))
mcmc_trace(as.array(fit_logistic3), pars = c('G', 'b1', 'k'))

sum_logistic <- summary(fit_logistic)
sum_logistic[[1]][c('L', 'k', 'x0'), ]

sum_logistic2 <- summary(fit_logistic2)
sum_logistic2[[1]][c('a', 'b1', 'b2', 'k'), ]

sum_logistic3 <- summary(fit_logistic3)
sum_logistic3[[1]][c('G', 'b1', 'k'), ]

# Create plot with data, fit, and credible interval of fit.

getpar <- function(fit) as.data.frame(do.call(cbind, extract(fit)))

par_logistic <- getpar(fit_logistic)
par_logistic2 <- getpar(fit_logistic2)
par_logistic3 <- getpar(fit_logistic3)

fn_logistic <- function(x, L, x0, k, ...) L/(1 + exp(-k*(x - x0)))
fn_logistic2 <- function(x, a, b1, b2, k, ...) (b1 * x ^ -a) / (1 + b2 * exp(-k * x))
fn_logistic3 <- function(x, G, b1, k, ...) G * (1 - b1 * exp(-k * x)) ^ 3


light_pred <- exp(seq(log(1.1), log(412), length.out = 101))

pred_logistic <- do.call(rbind, pmap(par_logistic, fn_logistic, x = light_pred))
pred_logistic2 <- do.call(rbind, pmap(par_logistic2, fn_logistic2, x = light_pred))
pred_logistic3 <- do.call(rbind, pmap(par_logistic3, fn_logistic3, x = light_pred))

getpred <- function(x, p, n, modelname) {
  out <- x %>%
    as.data.frame %>%
    map(~ quantile(., probs = p)) %>%
    do.call(rbind, .) %>%
    as.data.frame %>%
    setNames(nm = n)
  data.frame(model = modelname, cbind(light_area = light_pred, out))
}
  
pred_logistic_quant <- getpred(pred_logistic, c(0.025, 0.25, 0.5, 0.75, 0.975), c('q025', 'q25', 'q50', 'q75', 'q975'), 'logistic')
pred_logistic_quant2 <- getpred(pred_logistic2, c(0.025, 0.25, 0.5, 0.75, 0.975), c('q025', 'q25', 'q50', 'q75', 'q975'), 'modified logistic')
pred_logistic_quant3 <- getpred(pred_logistic3, c(0.025, 0.25, 0.5, 0.75, 0.975), c('q025', 'q25', 'q50', 'q75', 'q975'), 'bertalanffy')

pred_logistic_all <- rbind(pred_logistic_quant, pred_logistic_quant2, pred_logistic_quant3)

p_median <- ggplot(prod_light_bin) +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q25, yend = q75), size = 0.75) +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q025, yend = q975)) +
  geom_point(aes(x = bin_midpoint, y = median)) +
  geom_line(data = pred_logistic_all, aes(x = light_area, y = q50, color = model, group = model)) +
  geom_line(data = pred_logistic_all, aes(x = light_area, y = q025, color = model, group = model), linetype = 'dashed') +
  geom_line(data = pred_logistic_all, aes(x = light_area, y = q975, color = model, group = model), linetype = 'dashed') +
  scale_x_log10(name = 'Incoming light per area (W m-2)') + scale_y_log10('Production per area (kg y-1 m-2)') +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        legend.position = c(0.2, 0.8))

p_mean <- ggplot(prod_light_bin) +
  geom_errorbar(aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max), width = 0.1) +
  geom_point(aes(x = bin_midpoint, y = mean)) +
  geom_line(data = pred_logistic_all, aes(x = light_area, y = q50, color = model, group = model)) +
  geom_line(data = pred_logistic_all, aes(x = light_area, y = q025, color = model, group = model), linetype = 'dashed') +
  geom_line(data = pred_logistic_all, aes(x = light_area, y = q975, color = model, group = model), linetype = 'dashed') +
  scale_x_log10(name = 'Incoming light per area (W m-2)') + scale_y_log10('Production per area (kg y-1 m-2)') +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        legend.position = c(0.2, 0.8))

p_hist <- ggplot(prod_light_bin, aes(x = bin_midpoint, y = mean_n_individuals)) + geom_col() +
  scale_x_log10() +
  theme_void()

library(gridExtra)
png('C:/Users/Q/google_drive/ForestLight/figs/credible_interval_plots/production_vs_light_1995_median.png', height = 5, width = 6, res = 400, units = 'in')
grid.arrange(
  grobs = list(p_hist, p_median),
  widths = c(1, 24),
  heights = c(1, 5),
  layout_matrix = rbind(c(NA, 1),
                        c(2, 2))
)
dev.off()

png('C:/Users/Q/google_drive/ForestLight/figs/credible_interval_plots/production_vs_light_1995_mean.png', height = 5, width = 6, res = 400, units = 'in')
grid.arrange(
  grobs = list(p_hist, p_mean),
  widths = c(1, 24),
  heights = c(1, 5),
  layout_matrix = rbind(c(NA, 1),
                        c(2, 2))
)
dev.off()


# Fit by FG ---------------------------------------------------------------

# Create data

make_logistic_data <- function(x, n_sub) {
  if (n_sub < nrow(x)) {
    sample_rows <- sample(nrow(x), n_sub)
  } else {
    sample_rows <- 1:nrow(x)
    n_sub <- nrow(x)
  }
  with(x[sample_rows, ],
       list(N = n_sub, x = light_received/crownarea, y = production/crownarea))
}

n_sub <- 5000
set.seed(919)

dat <- mutate(dat, fg = paste0('fg', fg))
dat$fg[dat$fg == 'fgNA'] <- 'unclassified'

dat_logistic_fg <- dat %>%
  group_by(fg) %>%
  do(data = make_logistic_data(., n_sub))

NC <- 3
NI <- 6000
NW <- 5000

fit_bert_fg <- list()

for (i in 1:6) {
print(i)
fit_bert_fg[[i]] <- sampling(mod_logistic3, data = dat_logistic_fg$data[[i]], chains = NC, iter = NI, warmup = NW, seed = 117+i)

}

sum_bert_fg <- lapply(fit_bert_fg, function(x) summary(x)[[1]][c('G', 'b1', 'k'), ])

getpar <- function(fit) as.data.frame(do.call(cbind, extract(fit)))

par_bert_fg <- lapply(fit_bert_fg, getpar)

fn_logistic3 <- function(x, G, b1, k, ...) G * (1 - b1 * exp(-k * x)) ^ 3

light_pred <- exp(seq(log(1.1), log(412), length.out = 101))

pred_bert_fg <- lapply(par_bert_fg, function(z) do.call(rbind, pmap(z, fn_logistic3, x = light_pred)))

getpred <- function(x, p, n, modelname) {
  out <- x %>%
    as.data.frame %>%
    map(~ quantile(., probs = p)) %>%
    do.call(rbind, .) %>%
    as.data.frame %>%
    setNames(nm = n)
  data.frame(model = modelname, cbind(light_area = light_pred, out))
}

pred_bert_quant_fg <- lapply(pred_bert_fg, function(z) getpred(z, c(0.025, 0.25, 0.5, 0.75, 0.975), c('q025', 'q25', 'q50', 'q75', 'q975'), 'bertalanffy'))

pred_bert_quant_fg <- cbind(fg = rep(c('fg1','fg2','fg3','fg4','fg5','unclassified'), each = length(light_pred)),
                            do.call(rbind, pred_bert_quant_fg))

# Bin data by functional group to create the plot.
numbins <- 20
light_bins <- logbin(x = dat$light_area, n = numbins)

prod_light_bin_fg <- dat %>%
  group_by(fg) %>%
  do(bin = fakebin_across_years(dat_values = .$production_area,
                                       dat_classes = .$light_area,
                                       edges = light_bins,
                                       n_census = 1))

prod_light_bin_fg <- cbind(fg = rep(c('fg1','fg2','fg3','fg4','fg5','unclassified'), each = numbins),
                           do.call(rbind, prod_light_bin_fg$bin)) %>%
  filter(complete.cases(.))

x_loc <- 'bin_midpoint'

p_median <- ggplot(prod_light_bin_fg %>% filter(!fg %in% "unclassified")) +
  geom_segment(aes_string(x = x_loc, xend = x_loc, y = 'q25', yend = 'q75', group = 'fg'), size = 0.75) +
  geom_segment(aes_string(x = x_loc, xend = x_loc, y = 'q025', yend = 'q975', group = 'fg')) +
  geom_point(aes_string(x = x_loc, y = 'median', color = 'fg', group = 'fg')) +
  geom_line(data = pred_bert_quant_fg %>% filter(!fg %in% "unclassified"), aes(x = light_area, y = q50, color = fg, group = fg)) +
  geom_line(data = pred_bert_quant_fg %>% filter(!fg %in% "unclassified"), aes(x = light_area, y = q025, color = fg, group = fg), linetype = 'dashed') +
  geom_line(data = pred_bert_quant_fg %>% filter(!fg %in% "unclassified"), aes(x = light_area, y = q975, color = fg, group = fg), linetype = 'dashed') +
  scale_x_log10(name = 'Incoming light per area (W m-2)') + scale_y_log10('Production per area (kg y-1 m-2)') +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        legend.position = c(0.2, 0.8))

p_mean <- ggplot(prod_light_bin_fg %>% filter(!fg %in% "unclassified")) +
  geom_errorbar(aes_string(x = x_loc, ymin = 'ci_min', ymax = 'ci_max', group = 'fg'), width = 0.1) +
  geom_point(aes_string(x = x_loc, y = 'mean', color = 'fg', group = 'fg')) +
  geom_line(data = pred_bert_quant_fg %>% filter(!fg %in% "unclassified"), aes(x = light_area, y = q50, color = fg, group = fg)) +
  geom_line(data = pred_bert_quant_fg %>% filter(!fg %in% "unclassified"), aes(x = light_area, y = q025, color = fg, group = fg), linetype = 'dashed') +
  geom_line(data = pred_bert_quant_fg %>% filter(!fg %in% "unclassified"), aes(x = light_area, y = q975, color = fg, group = fg), linetype = 'dashed') +
  scale_x_log10(name = 'Incoming light per area (W m-2)') + scale_y_log10('Production per area (kg y-1 m-2)') +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        legend.position = c(0.2, 0.8))

slope_fg <- lapply(sum_bert_fg, function(x) c(q50 = x[3, 6], q025 = x[3,4], q975 = x[3,8]))
slope_fg <- data.frame(fg = c('fg1','fg2','fg3','fg4','fg5','unclassified'), do.call(rbind, slope_fg))

p_slope <- ggplot(slope_fg[1:5,], aes(x = fg, y = q50, ymin = q025, ymax = q975)) +
  geom_pointrange() +
  theme_classic() + 
  labs(x = 'Functional group', y = 'Bertalanffy slope')

fpfig <- 'C:/Users/Q/google_drive/ForestLight/figs/credible_interval_plots'
ggsave(file.path(fpfig, 'production_vs_light_1995_median_allFG.png'), p_median, height=5, width=6, dpi=400)
ggsave(file.path(fpfig, 'production_vs_light_1995_mean_allFG.png'), p_mean, height=5, width=6, dpi=400)
ggsave(file.path(fpfig, 'production_vs_light_1995_slopes.png'), p_slope, height=3, width=4, dpi=400)

# Histogram of number of individuals by fg
ggplot(prod_light_bin_fg %>% filter(!fg %in% "unclassified"), aes(x = bin_midpoint, y = mean_n_individuals)) + 
  theme_bw() + geom_col() + facet_grid(fg ~ ., scales = 'free_y') + scale_x_log10(name = 'Incoming light per area') + labs(y = 'Number of individuals in bin')
ggsave(file.path(fpfig, 'histograms_bylight.png'), height = 8, width = 5, dpi = 400)



# Plot individual functional groups with raw data -------------------------

prod_light_bin_fg <- filter(prod_light_bin_fg, !fg %in% 'unclassified')
pred_bert_quant_fg <- filter(pred_bert_quant_fg, !fg %in% 'unclassified')

p_median_byfg <- ggplot(prod_light_bin_fg) +
  facet_wrap(~ fg, scales = 'free') +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q25, yend = q75), size = 0.75) +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q025, yend = q975)) +
  geom_point(aes(x = bin_midpoint, y = median)) +
  geom_point(aes(x = bin_midpoint, y = mean), color = 'dodgerblue') +
  geom_line(data = pred_bert_quant_fg, aes(x = light_area, y = q50), color = 'red') +
  geom_line(data = pred_bert_quant_fg, aes(x = light_area, y = q025), linetype = 'dashed', color = 'red') +
  geom_line(data = pred_bert_quant_fg, aes(x = light_area, y = q975), linetype = 'dashed', color = 'red') +
  scale_x_log10(name = expression(paste('Incoming light per area (W m'^-2,')',sep=''))) + scale_y_log10(name = expression(paste('Growth per area (kg y'^-1, ' m'^-2,')', sep=''))) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA))

# Raw data
# Set color so that the data points look roughly the same amount of transparent.
p_raw_byfg <- ggplot(dat %>% filter(!fg %in% 'unclassified')) +
  facet_wrap(~ fg, scales = 'free') +
  geom_point(aes(x = light_area, y = production_area, alpha = fg)) +
  geom_line(data = pred_bert_quant_fg, aes(x = light_area, y = q50), color = 'red') +
  geom_line(data = pred_bert_quant_fg, aes(x = light_area, y = q025), linetype = 'dashed', color = 'red') +
  geom_line(data = pred_bert_quant_fg, aes(x = light_area, y = q975), linetype = 'dashed', color = 'red') +
  scale_alpha_manual(values = c(0.1, 0.08, 0.02, 0.02, 0.004)) +
  scale_x_log10(name = expression(paste('Incoming light per area (W m'^-2,')',sep=''))) + scale_y_log10(name = expression(paste('Growth per area (kg y'^-1, ' m'^-2,')', sep=''))) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA),
        legend.position = 'none')

fpfig <- 'C:/Users/Q/google_drive/ForestLight/figs/credible_interval_plots'
ggsave(file.path(fpfig, 'production_vs_light_1995_separatefg_median.png'), p_median_byfg, height=6, width=9, dpi=400)
ggsave(file.path(fpfig, 'production_vs_light_1995_separatefg_rawdata.png'), p_raw_byfg, height=6, width=9, dpi=400)
