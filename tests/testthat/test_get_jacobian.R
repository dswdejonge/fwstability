context("Jacobian matrix creation with no optional arguments")
# No optional arguments

#### Proper data format ####
fwnames <- c("PLANT", "WORM", "ANT")
FM <- matrix(c(0, 3, 2, 5, 0, 3, 0, 5, 0), nrow = 3, ncol = 3)
rownames(FM) <- fwnames ; colnames(FM) <- fwnames
BM <- c(30, 20, 10) ; names(BM) <- fwnames
AE <- c(0.1, 0.2, 0.3) ; names(AE) <- fwnames
GE <- c(0.1, 0.2, 0.3) ; names(GE) <- fwnames
JM <- matrix(c(0,
               AE[1] * GE[1] * FM[2,1] / BM[2] + -FM[1,2] / BM[2],
               AE[1] * GE[1] * FM[3,1] / BM[3],

               AE[2] * GE[2] * FM[1,2] / BM[1] + -FM[2,1] / BM[1],
               0,
               AE[2] * GE[2] * FM[3,2] / BM[3] + -FM[2,3] / BM[3],

               -FM[3,1] / BM[1],
               AE[3] * GE[3] * FM[2,3] / BM[2] + -FM[3,2] / BM[2],
               0
), nrow = 3, ncol = 3)
rownames(JM) <- fwnames ; colnames(JM) <- fwnames

test_that("the function works without optional arguments", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE),
               JM)
})
