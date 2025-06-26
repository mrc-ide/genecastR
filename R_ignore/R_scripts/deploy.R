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

#set.seed(9)

# ------------------------------------------------------------------

# define parameters
s <- 0.01
p0 <- 0.2
N <- 1e3
t_min <- 1
t_max <- 365*2
n_samp <- 100
n_pop <- 9
sigma <- 1/sqrt(N)
dx <- 0.01
nx <- 1/dx + 1

# QUESTION
#20 pm for 15 sites?
#40 pm for 8 sites?
# prior on selection coefficient for k13
# starting frequency of 7.5%

# REMAINING TODO
# - Assertions/checks
# - MCMC posterior plots
# - Wald CI
# - Example data
# - Continuous integration
# - Branch and version control
# - Tests
# - Code coverage
# - Package documentation
#     - Estimating selection
#     - Predicting allele prevalence
#     - Power and sample size calculation, with & without data

# simulate prevalence
df_sim_WF <- sim_WF(s = s, p0 = p0, N = N, t_max = t_max, n_pop = n_pop)

# sample data
df_data <- expand_grid(pop = 1:n_pop, t = 50*(1:4), n_samp = n_samp) |>
  sample_prev(df_sim = df_sim_WF)

# estimate params via ML
df_ML <- estimate_ML(df_data = df_data)

df_sim_WF |>
  left_join(df_ML) |>
  mutate(p_logistic = get_logistic(s = s, p0 = p0, t = t)) |>
  ggplot() + theme_bw() +
  geom_line(aes(x = t, y = p), col = grey(0.7)) +
  geom_line(aes(x = t, y = p_logistic), col = "red") +
  geom_pointrange(aes(x = t, y = p_est, ymin = CI_lower, ymax = CI_upper), data = df_data) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  facet_wrap(~pop)

# estimate parameters via MCMC
mcmc <- run_mcmc(df_data = df_data,
                 burnin = 1e2,
                 samples = 1e3)

#drjacoby::plot_trace(mcmc)
#drjacoby::plot_pairs(mcmc)

mcmc$output |>
  ggplot() + theme_bw() +
  geom_density_2d_filled(aes(x = s, y = sigma), contour_var = "density", bins = 9) +
  scale_fill_brewer(palette = "Reds") +
  scale_x_continuous(limits = c(0, 0.015), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.08), expand = c(0, 0)) +
  guides(fill = "none")

estimate_Bayesian(mcmc)

# solve HMM
df_mcmc_draws <- sample_mcmc(mcmc, n_draws = 100)
df_post <- get_posterior_prev(df_data = df_data,
                              s = df_mcmc_draws$s,
                              #s = 0.00,
                              sigma = df_mcmc_draws$sigma,
                              #sigma <- 0.03,
                              dx = dx,
                              t_min = t_min,
                              t_max = t_max)

# plot
#plot_prev(df_post = df_post)
#plot_prev(df_post = df_post, df_data = df_data)
plot_prev(df_post = df_post, df_trajectory = df_sim_WF, df_data = df_data)


# sample size based on MOE
get_MOE(pop = 5, t_samp = 200, n_samp = 20, df_post = df_post)
get_sample_size_MOE(pop = 1, t_samp = 300, MOE = 0.1, n_max = 1e3, df_post = df_post)

# sample size based on z-test
get_power_ztest(pop = 2, t_samp = 300, n_samp = 18, prev_thresh = 0.05, df_post = df_post)
get_sample_size_ztest(pop = 2, t_samp = 300, prev_thresh = 0.05, df_post = df_post)

