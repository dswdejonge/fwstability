context("Jacobian matrix creation")

source("models/model_dead.R")
test_that("correct JM created with dead compartments", {
  expect_equal(getJacobian(model),
               JM)
})
test_that("getJacobian incorporates vector diagonal correctly", {
  expect_equal(getJacobian(model2),
               JM2)
})
test_that("getJacobian calculates the diagonal from flux values", {
  expect_equal(getJacobian(model3),
               JM3)
})

source("models/model_dead_externals.R")
test_that("correct JM created with both dead and external arguments", {
  expect_equal(getJacobian(model),
               JM)
})

source("models/model_multiple_dead.R")
test_that("correct JM created with multiple dead compartments but one defecation", {
  expect_equal(getJacobian(model),
               JM)
})

source("models/model_mort_def_same_comp.R")
test_that("correct JM created with defecation and mortality into same compartment", {
  expect_equal(getJacobian(model),
               JM)
})

source("models/model_multiple_def.R")
test_that("correct JM created with multiple dead and defecation compartments", {
  expect_equal(getJacobian(model),
               JM)
})
