context("Jacobian matrix creation")
library(fwstability)

# Test data
FM <- matrix(1:9, nrow = 3, ncol = 3)
rownames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
BM <- c(10, 20, 30)
AE <- c(0.1, 0.2, 0.3)
GE <- c(0.1, 0.2, 0.3)



test_that("the function works with and without optional arguments", {
  expect_equal(dim(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE)),
               c(3,3))
  expect_equal(dim(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE, dead = "DETRITUS")),
               c(3,3))
  expect_equal(dim(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE, externals = "DETRITUS")),
               c(2,2))
})


test_that("the function throws an error if not all required data is given", {
  expect_error(getJacobian(BM = BM, AE = AE, GE = GE, dead = "DETRITUS"),
               "argument \"FM\" is missing, with no default")
  expect_error(getJacobian(FM = FM, AE = AE, GE = GE, dead = "DETRITUS"),
               "argument \"BM\" is missing, with no default")
  expect_error(getJacobian(FM = FM, BM = BM, GE = GE, dead = "DETRITUS"),
               "argument \"AE\" is missing, with no default")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, dead = "DETRITUS"),
               "argument \"GE\" is missing, with no default")
})


test_that("the function only executes with a square matrix", {
  expect_error(getJacobian(FM[,1:2], BM, AE, GE, dead = "DETRITUS"),
               "flow matrix is not square")
})


test_that("the function only executes if all biomasses are greater than zero", {
  expect_error(getJacobian(FM, BM = c(10,0,30), AE, GE, dead = "DETRITUS"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobian(FM, BM = c(10,-10,30), AE, GE, dead = "DETRITUS"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobian(FM, BM = c(10,NA,30), AE, GE, dead = "DETRITUS"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobian(FM, BM = c("10","NA","30"), AE, GE, dead = "DETRITUS"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")

})
