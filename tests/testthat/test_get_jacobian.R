context("Jacobian matrix creation")
library(fwstability)

# Test data
FM <- matrix(c(0, 3, 8, 7, 0, 0, 4, 4, 0), nrow = 3, ncol = 3)
BM <- c(30, 20, 10)
AE <- c(0.1, 0.2, 0.3)
GE <- c(0.1, 0.2, 0.3)

rownames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
names(BM) <- c("DETRITUS", "PLANT", "ANIMAL")
names(AE) <- c("DETRITUS", "PLANT", "ANIMAL")
names(GE) <- c("DETRITUS", "PLANT", "ANIMAL")

results <- c(0,
             (FM[2,1] - FM[1,2] + FM[2,3]*(1-AE[3])) / BM[2],
             (FM[3,1] - FM[1,3]) / BM[3],
             AE[2] * GE[2] * FM[1,2] / BM[1],
             0,
             -FM[2,3] / BM[3],
             AE[3] * GE[3] * FM[1,3] / BM[1],
             AE[3] * GE[3] * FM[2,3] / BM[2],
             0
             )
JM <- matrix(results, nrow = 3, ncol = 3)
rownames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")

test_that("the function works with and without optional arguments", {
  expect_equal(dim(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE)),
               c(3,3))
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE, dead = "DETRITUS"),
               JM)
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


#test_that("the function checks the indices of dead compartments against the presence of NA in the AE and GE vector", {

#})
