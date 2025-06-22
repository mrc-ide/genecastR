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

# define params
s <- 0.01
p0 <- 0.2
N <- 1e3
t_max <- 365
n_deme <- 3
sigma <- 1/sqrt(N)
dx <- 0.01

df_sim_WF <- sim_WF(s = s, p0 = p0, N = N, t_max = t_max, n_deme = n_deme)

df_sim_diff <- sim_diffusion(s = s, p0 = p0, sigma = sigma, dx = dx, t_max = t_max, n_deme = n_deme)

df_sample <- expand_grid(pop = 1:2, t = c(50, 100, 150), n_samp = 1e2)
sample_prev(df_sim = df_sim_WF, df_sample = df_sample)
