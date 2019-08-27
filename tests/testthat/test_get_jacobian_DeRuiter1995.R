context("Jacobian matrix creation with data from De Ruiter et al. 1995")

# Proper data format, from De Ruiter et al. 1995
# Bacteria, Fungi and Enchytraeids take up detritus.
# Enchytraeids also feed on Fungi and Bacteria.
#
#     FUN
#    /   \
# DET -- ENCH
#    \   /
#     BAC

fwnames <- c("DET", "BAC", "FUN", "ENCH")
FM <- matrix(c(
  0, 1606, 40.8, 11.1,
  0, 0,    0,    10.1,
  0, 0,    0,    0.09,
  0, 0,    0,     0),
  nrow = 4, ncol = 4, byrow = T)
rownames(FM) <- fwnames ; colnames(FM) <- fwnames
BM <- c(2500, 227.5, 2.13, 0.43) ; names(BM) <- fwnames
AE <- c(NA, 1.0, 1.0, 0.25) ; names(AE) <- fwnames
GE <- c(NA, 0.30, 0.30, 0.40) ; names(GE) <- fwnames
dead <- list(names = "DET", def = "noDef")
model <- list(
  type = "EF", FM = FM, BM = BM, AE = AE, GE = GE, dead = dead
)
# Expected answer
JM <- matrix(c(0,
               - FM[1,2] / BM[2],
               - FM[1,3] / BM[3],
               - FM[1,4] / BM[4],

               AE[2] * GE[2] * FM[1,2] / BM[1],
               0,
               0,
               - FM[2,4] / BM[4],

               AE[3] * GE[3] * FM[1,3] / BM[1],
               0,
               0,
               - FM[3,4] / BM[4],

               AE[4] * GE[4] * FM[1,4] / BM[1],
               AE[4] * GE[4] * FM[2,4] / BM[2],
               AE[4] * GE[4] * FM[3,4] / BM[3],
               0
), nrow = 4, ncol = 4)
rownames(JM) <- fwnames ; colnames(JM) <- fwnames

test_that("the function works with data from literature", {
  expect_equal(getJacobian(model),
               JM)
})
