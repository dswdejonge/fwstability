context("Jacobian matrix creation")

source("models/model_regular.R")
test_that("the function works without optional arguments", {
  expect_equal(getJacobian(model), JM)
  # netto FM
  expect_equal(getJacobian(model1), JM1)
})

source("models/model_dead.R")
test_that("correct JM created with dead compartments", {
  expect_equal(getJacobian(model), JM)
  # netto FM
  expect_equal(getJacobian(model4), JM4)
})
test_that("getJacobian incorporates vector diagonal correctly", {
  expect_equal(getJacobian(model2), JM2)
})
test_that("getJacobian calculates the diagonal from flux values", {
  expect_equal(getJacobian(model3), JM3)
})

source("models/model_externals.R")
test_that("correct JM created with external compartments", {
  expect_equal(getJacobian(model), JM)
})

source("models/model_dead_externals.R")
test_that("correct JM created with both dead and external arguments", {
  expect_equal(getJacobian(model), JM)
})

source("models/model_multiple_dead.R")
test_that("correct JM created with multiple dead compartments but one defecation", {
  expect_equal(getJacobian(model), JM)
})

source("models/model_mort_def_same_comp.R")
test_that("correct JM created with defecation and mortality into same compartment", {
  expect_equal(getJacobian(model),
               JM)
})

source("models/model_multiple_def.R")
test_that("correct JM created with multiple dead and defecation compartments", {
  expect_equal(getJacobian(model), JM)
})

source("models/model_Omni_Pred.R")
test_that("correct JM created with reciprocal predation and detritus", {
  expect_equal(getJacobian(model), JM)
})

source("models/model_deRuiter1995small.R")
test_that("correct JM created with data from De Ruiter et al. 1995", {
  expect_equal(getJacobian(model), JM)
})

# Take value visually from paper
source("models/model_deRuiter1995.R")
JM <- getJacobian(model)
a1 <- JM["Predatory_collembola", "Predatory_mites"]
a2 <- JM["Predatory_mites", "Predatory_collembola"]

test_that("correct JM created with data from De Ruiter et al. 1995", {
  expect_equal(a1 > 0.1 & a1 < 0.15, TRUE)
  expect_equal(a2 < 0 & a2 > -1, TRUE)
})

source("models/model_Lotka_Volterra.R")
test_that("The package works with a dynamic LV function", {
  expect_equal(getStability(getJacobian(model)), answer)
})

source("models/model_use_diff_dead.R")
test_that("correct JM if species deposit in different compartments", {
  expect_equal(getJacobian(model), JM)
})
