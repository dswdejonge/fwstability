context("Extracting data from a LIM")
source("models/model_LIM.R")

tagTest <- c(1,2,3,4,5)
names(tagTest) <- c("myAtag", "myBtag", "myC", "myD", "myE")
exp_result <- tagTest[1:2]; names(exp_result) <- c("MYA", "MYB")
test_that("the getTag function extracts the correct information", {
  expect_equal(getTag(tagTest, tag = "tag"), exp_result)
})

test_that("the Flowmatrix function works with parallel flows", {
  #expect_equal(Flowmatrix(lim, web = lim_solved$X), FM)
  expect_equal(getFlowMatrix(readLIM), FM)
  expect_equal(getFlowMatrix(readLIM, web = lim_solved$X), FM)
})

test_that("the getNettoFM provides correct answer", {
  expect_equal(getNettoFM(FM, c("DEADLABILE", "DEADREFRAC")), FM1)
})

test_that("the getVariables function gives right answer", {
  expect_equal(getVariables(readLIM, web = lim_solved$X), variables)
  expect_equal(getVariables(readLIM), variables)
})

test_that("the right conversion efficiencies are extracted", {
  expect_equal(getCE(FM = FM, vars = variables, lim = lim, aTag = "ass", gTag = "growth")$AE,
               AE)
  expect_equal(getCE(FM = FM, vars = variables, lim = lim, aTag = "ass", gTag = "growth")$GE,
               GE)
  expect_equal(getCE(FM = FM, vars = variables, lim = lim)$AE,
               AE)
  expect_equal(getCE(FM = FM, vars = variables, lim = lim)$GE,
               GE)
})

test_that("the right Jacobian matrix is produced from the LIM", {
  expect_equal(getJacobian(model)$JM, JM)
  # Use netto FM
  expect_equal(getJacobian(model1)$JM, JM1)
})

