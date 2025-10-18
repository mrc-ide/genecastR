# deploy.R
#
# Author: Bob Verity
# Date: 2025-06-22
#
# Purpose:
# Sandbox for genecastR package.
#
# ------------------------------------------------------------------

library(ggplot2)
library(tidyverse)

set.seed(1)

# ------------------------------------------------------------------

# define parameters
s <- 0.5
p0 <- 0.1
N <- 1e3
generations <- 52*5
t_step <- 7
n_samp <- 100
n_pop <- 9
sigma <- 1/sqrt(N)
dx <- 0.01
nx <- 1/dx + 1
sigma_true <- sqrt(365 / (N*dt))

# simulate prevalence
df_sim_WF <- sim_WF(s = s, p0 = p0, N = N, generations = generations, n_pop = n_pop, t_step = t_step)

df_sim_WF |>
  ggplot() + theme_bw() +
  geom_line(aes(x = t, y = p, col = factor(pop))) +
  ggtitle("Wright-Fisher simulation")

# or simulate from diffusion approximation
df_sim_diff <- sim_diffusion(s = s, p0 = p0, sigma = sigma_true, generations = generations,
                             n_pop = n_pop, t_step = t_step)

df_sim_diff |>
  ggplot() + theme_bw() +
  geom_line(aes(x = t, y = p, col = factor(pop))) +
  ggtitle("Diffusion simulation")

# sample data
df_sample <- expand_grid(pop = 1:n_pop, gen = seq(0, 52*5, by = 16), n_samp = n_samp)
df_data <- sample_prev(df_sim = df_sim_WF, df_sample = df_sample)

# plot raw data
plot_prev(df_data = df_data)

# estimate params via ML
df_ML <- estimate_ML(df_data = df_data, s_range = c(-2, 2))
df_ML

# check ML fit
df_sim_WF |>
  left_join(df_ML) |>
  mutate(t_numeric = as.numeric(t - min(t)) / 365,
         p_logistic = get_logistic(s = s, p0 = p0, t = t_numeric)) |>
  ggplot() + theme_bw() +
  geom_line(aes(x = t, y = p), col = "firebrick2") +
  geom_line(aes(x = t, y = p_logistic), col = "#2171B5") +
  geom_pointrange(aes(x = t, y = p_est, ymin = CI_lower, ymax = CI_upper), size = 0.1,
                  data = get_CIs(df_data)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  facet_wrap(~pop)

# estimate parameters via MCMC
mcmc <- run_mcmc(df_data = df_data,
                 burnin = 1e2,
                 samples = 1e3,
                 dt = 7)

#drjacoby::plot_trace(mcmc)
#drjacoby::plot_pairs(mcmc)

plot_mcmc_scatter(mcmc)
plot_mcmc_density(mcmc)

estimate_Bayesian(mcmc)

# solve HMM
df_mcmc_draws <- sample_mcmc(mcmc, n_draws = 100)
df_post <- get_posterior_prev(df_data = df_data,
                              s = df_mcmc_draws$s,
                              #s = 0.00,
                              sigma = df_mcmc_draws$sigma,
                              #sigma <- 0.5,
                              dt = 7,
                              dx = dx)

# plot
#plot_prev(df_post = df_post, df_data = df_data)
plot_prev(df_post = df_post, df_trajectory = df_sim_WF, df_data = df_data)


# sample size based on MOE
sample_date <- as.Date("2020-05-20")
get_MOE(pop = 1, t = sample_date, n_samp = 100, df_post = df_post)
get_sample_size_MOE(pop = 1, t = sample_date, MOE = 0.1, n_max = 1e3, df_post = df_post)

# sample size based on z-test
get_power_ztest(pop = 1, t = sample_date, n_samp = 100, prev_thresh = 0.05, df_post = df_post)
get_sample_size_ztest(pop = 1, t = sample_date, prev_thresh = 0.05, df_post = df_post)

