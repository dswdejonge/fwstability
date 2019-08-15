context("Jacobian matrix creation with dead and external compartments")
### Both dead compartments and externals.
# Detritus takes up CO2.
# Plants and animals take up detritus.
# Plants are eaten by animals.
# Animals and plants both deposit material into detritus compartment.
# Animals respire into the external CO2 compartment.
#
#                 V---3---\
# CO2 -1-> DETRITUS -7-> PLANT
#  ^            ^ \ 4     | 4
#  \ 1        8 \ v       v
#   -----------  #  ANIMAL  #

fwnames <- c("DETRITUS", "PLANT", "ANIMAL", "CO2")
FM <- matrix(c(0, 3, 8, 1, 7, 0, 0, 0, 4, 4, 0, 0, 0, 0, 1, 0), nrow = 4, ncol = 4)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
BM <- c(30, 20, 10) ; names(BM) <- fwnames[1:3]
AE <- c(NA, 0.2, 0.3) ; names(AE) <- fwnames[1:3]
GE <- c(NA, 0.2, 0.3) ; names(GE) <- fwnames[1:3]
dead <- list("DETRITUS", "Def")
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
                           dead = dead, externals = "CO2"),
               JM)
})
