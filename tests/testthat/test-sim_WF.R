test_that("sim_WF returns expected structure and types", {
  set.seed(1)
  result <- sim_WF(s = 0.1, p0 = 0.5, N = 100, generations = 52, n_pop = 3,
                   t_start = as.Date("2020-01-01", t_step = 7))

  expect_s3_class(result, "data.frame")
  expect_named(result, c("pop", "gen", "t", "p"))
  expect_type(result$pop, "integer")
  expect_type(result$gen, "integer")
  expect_type(result$t, "double")
  expect_type(result$p, "double")
})

test_that("sim_WF output has correct number of rows", {
  set.seed(1)
  generations <- 20
  n_pop <- 5
  result <- sim_WF(s = 0, p0 = 0.5, N = 100, generations = generations, n_pop = n_pop,
                   t_start = as.Date("2020-01-01", t_step = 7))
  expect_equal(nrow(result), generations * n_pop)
})

test_that("allele frequencies are within [0, 1]", {
  set.seed(1)
  result <- sim_WF(s = 0.05, p0 = 0.3, N = 100, generations = 15, n_pop = 2,
                   t_start = as.Date("2020-01-01", t_step = 7))
  expect_true(all(result$p >= 0))
  expect_true(all(result$p <= 1))
})

test_that("sim_WF throws error for invalid inputs", {
  set.seed(1)
  expect_error(sim_WF(s = "strong", p0 = 0.5, N = 100, generations = 10, n_pop = 1))
  expect_error(sim_WF(s = 0.1, p0 = -0.1, N = 100, generations = 10, n_pop = 1))
  expect_error(sim_WF(s = 0.1, p0 = 0.5, N = 0, generations = 10, n_pop = 1))
  expect_error(sim_WF(s = 0.1, p0 = 0.5, N = 100, generations = 0, n_pop = 1))
  expect_error(sim_WF(s = 0.1, p0 = 0.5, N = 100, generations = 10, n_pop = -2))
})
