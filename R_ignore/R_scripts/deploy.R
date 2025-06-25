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

#set.seed(2)

# ------------------------------------------------------------------

# define parameters
s <- 0.005
p0 <- 0.1
N <- 1e3
t_max <- 365*2
n_samp <- 100
n_pop <- 3
sigma <- 1/sqrt(N)
dx <- 0.01
nx <- 1/dx + 1

# QUESTION
#20 pm for 15 sites?
#40 pm for 8 sites?
# prior on selection coefficient for k13
# starting frequency of 7.5%

# simulate prevalence
df_sim_WF <- sim_WF(s = s, p0 = p0, N = N, t_max = t_max, n_pop = n_pop)

# sample data
df_data <- expand_grid(pop = 1:n_pop, t = 100*(1:4), n_samp = n_samp) |>
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

data_list <- df_data |>
  select(pop, t, n_samp, n_pos) |>
  mutate(t = as.integer(t),
         n_samp = as.integer(n_samp),
         n_pos = as.integer(n_pos)) |>
  filter(pop == 1) |>
  #filter(t == 100) |>
  group_split(pop, .keep = FALSE)

x_init <- rep(1/nx, nx)

t0 <- Sys.time()
for (i in 1:100) {
  loglike_cpp11(params = c(s = 0.01, sigma = 0.02),
                data = data_list,
                misc = list(dx = 0.01,
                            x_init = x_init))
}
Sys.time() - t0

t0 <- Sys.time()
for (i in 1:100) {
  loglike_R(params = c(s = 0.01, sigma = 0.02),
            data = list(df = data_list),
            misc = list(dx = 0.01,
                        x_init = x_init,
                        t_min = 1))
}
Sys.time() - t0




# estimate parameters via MCMC
mcmc <- run_mcmc(df_data = df_data,
                 burnin = 1e2,
                 samples = 1e2)

#drjacoby::plot_trace(mcmc)
#drjacoby::plot_pairs(mcmc)

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
plot_prev(df_post = df_post, df_trajectory = df_sim_WF, df_data = df_data)

# sample size based on MOE
get_MOE(pop = 1, t_samp = 500, n_samp = 100, df_post = df_post)
get_sample_size_MOE(pop = 1, t_samp = 500, MOE = 0.1, n_max = 1e3, df_post = df_post)

# sample size based on z-test
get_power_ztest(pop = 3, t_samp = 600, n_samp = 50, prev_thresh = 0.05, df_post = df_post)
get_sample_size_ztest(pop = 1, t_samp = 600, prev_thresh = 0.05, df_post = df_post)

