context("Jacobian matrix creation")

source("models/model_regular.R")
test_that("the function works without optional arguments", {
  expect_equal(getJacobian(model, diagonal = 0), JM)
  # netto FM
  expect_equal(getJacobian(model1, diagonal = 0), JM1)
})

source("models/model_dead.R")
test_that("correct JM created with dead compartments", {
  expect_equal(getJacobian(model, diagonal = 0), JM)
  # netto FM
  expect_equal(getJacobian(model4, diagonal = 0), JM4)
})
test_that("getJacobian incorporates vector diagonal correctly", {
  expect_equal(getJacobian(model2, diagonal = model2$diagonal), JM2)
})
test_that("getJacobian calculates the diagonal from flux values", {
  expect_equal(getJacobian(model3, diagonal = model3$diagonal), JM3)
})
test_that("a matrix with dead compartments gets normalized correctly", {
  expect_error(normalizeJacobian(JM1),
               "No zeroes may be present on the diagonal when normalizing the Jacobian matrix.")
  expect_equal(normalizeJacobian(JM2, allzero = F), JMnorm)
  expect_equal(normalizeJacobian(JM2, dead_names = "DETRITUS"), JMnorm2)
})


source("models/model_externals.R")
test_that("correct JM created with external compartments", {
  expect_equal(getJacobian(model, diagonal = 0), JM)
})

source("models/model_dead_externals.R")
test_that("correct JM created with both dead and external arguments", {
  expect_equal(getJacobian(model, diagonal = 0), JM)
})

source("models/model_multiple_dead.R")
test_that("correct JM created with multiple dead compartments but one defecation", {
  expect_equal(getJacobian(model, diagonal = model$diagonal), JM)
})

source("models/model_mort_def_same_comp.R")
test_that("correct JM created with defecation and mortality into same compartment", {
  expect_equal(getJacobian(model, diagonal = 0),
               JM)
})

source("models/model_multiple_def.R")
test_that("correct JM created with multiple dead and defecation compartments", {
  expect_equal(getJacobian(model, diagonal = 0), JM)
})

source("models/model_Omni_Pred.R")
test_that("correct JM created with reciprocal predation and detritus", {
  expect_equal(getJacobian(model, diagonal = 0), JM)
})

source("models/model_use_diff_dead.R")
test_that("correct JM if species deposit in different compartments", {
  expect_equal(getJacobian(model, diagonal = 0), JM)
})

source("models/model_use_diff_dead2.R")
test_that("correct JM if species deposit in different compartments", {
  expect_equal(getJacobian(model, diagonal = 0), JM)
})
