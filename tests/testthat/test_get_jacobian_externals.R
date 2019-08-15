context("Jacobian matrix creation with externals")
### Externals

# Proper data format
fwnames <- c("PLANT", "WORM", "ANT", "CO2")
FM <- matrix(c(0, 3, 2, 1, 5, 0, 3, 0, 0, 5, 0, 0, 0, 0, 1, 0), nrow = 4, ncol = 4)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
BM <- c(30, 20, 10) ; names(BM) <- fwnames[1:3]
AE <- c(0.1, 0.2, 0.3) ; names(AE) <- fwnames[1:3]
GE <- c(0.1, 0.2, 0.3) ; names(GE) <- fwnames[1:3]
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
rownames(JM) <- fwnames[1:3]
colnames(JM) <- fwnames[1:3]

test_that("the function works with external compartments", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE, externals = "CO2"),
               JM)
})
