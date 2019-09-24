context("Jacobian matrix creation")

source("model_dead_externals.R")
test_that("correct JM created with both dead and external arguments", {
  expect_equal(getJacobian(model),
               JM)
})
