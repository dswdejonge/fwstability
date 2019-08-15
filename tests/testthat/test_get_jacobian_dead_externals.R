context("Jacobian matrix creation with dead and external compartments")
### Both dead compartments and externals.

# Proper data format
fwnames <- c("DETRITUS", "PLANT", "ANIMAL", "CO2")
FM <- matrix(c(0, 3, 8, 1, 7, 0, 0, 0, 4, 4, 0, 0, 0, 0, 1, 0), nrow = 4, ncol = 4)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
BM <- c(30, 20, 10) ; names(BM) <- fwnames[1:3]
AE <- c(NA, 0.2, 0.3) ; names(AE) <- fwnames[1:3]
GE <- c(NA, 0.2, 0.3) ; names(GE) <- fwnames[1:3]
JM <- matrix(c(0,
               (FM[2,1] - FM[1,2] + FM[2,3]*(1-AE[3])) / BM[2],
               (FM[3,1] - FM[1,3]) / BM[3],
               AE[2] * GE[2] * FM[1,2] / BM[1],
               0,
               -FM[2,3] / BM[3],
               AE[3] * GE[3] * FM[1,3] / BM[1],
               AE[3] * GE[3] * FM[2,3] / BM[2],
               0
), nrow = 3, ncol = 3)
rownames(JM) <- fwnames[1:3]
colnames(JM) <- fwnames[1:3]

test_that("the function works correctly with both dead and external arguments", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = "DETRITUS", externals = "CO2"),
               JM)
})

test_that("the function only executes when the dead and external compartments have an existing name", {
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = "CARCASS", externals = "CO2"),
               "the names of the dead compartments are unknown")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = "DETRITUS", externals = "CARCASS"),
               "the names of the external compartments are unknown")
})
