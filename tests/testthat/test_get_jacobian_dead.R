context("Jacobian matrix creation with dead compartments")
# Dead compartment, no externals
# Plants and animals both take up detritus and deposit into detritus.
# Animals also eat plants.
#
#       /---------4----------v
# DETRITUS -7-> PLANT -4-> ANIMAL
#  ^______3______/ <____8___/
#
## Proper data format
fwnames <- c("DETRITUS", "PLANT", "ANIMAL")
FM <- matrix(c(0, 3, 8, 7, 0, 0, 4, 4, 0), nrow = 3, ncol = 3)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
BM <- c(30, 20, 10) ; names(BM) <- fwnames
AE <- c(NA, 0.2, 0.3) ; names(AE) <- fwnames
GE <- c(NA, 0.2, 0.3) ; names(GE) <- fwnames
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
rownames(JM) <- fwnames ; colnames(JM) <- fwnames

test_that("the function works with dead compartments", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE, dead = dead),
               JM)
})

# Proper data format, include vector diagonal
DIAG <- c(-1,-2,-3)
JM2 <- JM; JM2[c(1,5,9)] <- DIAG

test_that("the function incorporates vector diagonal correctly", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = dead, diagonal = DIAG),
               JM2)
})

# Proper data format, include vector calculated from model
MR <- c(NA, 10, 5) ; names(MR) <- fwnames
JM3 <- JM; JM3[c(1,5,9)] <- c(-1/BM[1] * (AE[2] * FM[1,2] + AE[3] * FM[1,3]),
                              -MR[2] / BM[2],
                              -MR[3] / BM[3]
)
test_that("the function calculates the diagonal from flux values", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = dead, diagonal = "model", MR = MR),
               JM3)
})
