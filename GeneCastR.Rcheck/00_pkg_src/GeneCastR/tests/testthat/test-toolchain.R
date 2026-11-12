test_that("the compiled toolchain is present", {
  tc <- gc_toolchain()
  expect_type(tc$openmp, "logical")
  expect_gte(tc$max_threads, 1L)
  expect_gte(tc$cxx_standard, 201703)
  if (!tc$openmp) {
    skip("built without OpenMP; the per-deme loop will run serially")
  }
})
