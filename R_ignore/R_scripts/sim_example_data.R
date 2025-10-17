# sim_example_data.R
#
# Author: Bob Verity
# Date: 2025-10-17
#
# Inputs: (none)
#
# Outputs: data/example_data.rda
#
# Purpose:
# Creates example dataset to be included with package by simulating from WF model.
#
# ------------------------------------------------------------------

library(tidyverse)
library(here)

set.seed(13)

# define simulation parameters
s <- 0.5
p0 <- 0.1
N <- 1e3
generations <- 52*5 + 1
samp_generations <- 1 + 52*(0:5)
n_samp <- length(samp_generations)
n_pop <- 3
t_start <- as.Date("2020-01-01")
t_step <- 7
n_samp <- sample(seq(50, 150, by = 10), n_samp, replace = TRUE)

# simulate prevalence
df_sim_WF <- sim_WF(s = s, p0 = p0, N = N, generations = generations,
                    n_pop = n_pop, t_start = t_start, t_step = t_step)

df_sim_WF |>
  ggplot() + theme_bw() +
  geom_line(aes(x = t, y = p, col = factor(pop))) +
  ggtitle("Wright-Fisher simulation")

# sample data
df_sample <- expand_grid(pop = 1:n_pop, gen = samp_generations) |>
  mutate(n_samp = rep(n_samp, times = n_pop))
df_data <- sample_prev(df_sim = df_sim_WF, df_sample = df_sample)

# plot
plot_prev(df_data = df_data)

# format and save
example_data <- df_data |>
  select(pop, t, n_samp, n_pos)

save(example_data, file = here("data", "example_data.rda"))
