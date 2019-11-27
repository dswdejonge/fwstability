context("Diagonal value calculation")

#### getDetritusDiagonal

# Proper data
FM <- matrix(c(0, 3, 8, 7, 0, 0, 4, 4, 0), nrow = 3, ncol = 3)
rownames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
BM <- c(30, 20, 10) ; names(BM) <- c("DETRITUS", "PLANT", "ANIMAL")
AE <- c(NA, 0.2, 0.3) ; names(AE) <- c("DETRITUS", "PLANT", "ANIMAL")
GE <- c(NA, 0.2, 0.3) ; names(GE) <- c("DETRITUS", "PLANT", "ANIMAL")
dead <- list(names = "DETRITUS")
result <- -1/BM[1] * (AE[2] * FM[1,2] + AE[3] * FM[1,3])

test_that("getDiagonalDetritus works with correct data input", {
  expect_equal(getDiagonalDetritus(FM = FM, BM = BM, AE = AE, dead = dead),
               result)
})



#### getDiagonal
# Proper data
FM <- matrix(c(0, 3, 8, 7, 0, 0, 4, 4, 0), nrow = 3, ncol = 3)
rownames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
BM <- c(30, 20, 10) ; names(BM) <- c("DETRITUS", "PLANT", "ANIMAL")
AE <- c(NA, 0.2, 0.3) ; names(AE) <- c("DETRITUS", "PLANT", "ANIMAL")
GE <- c(NA, 0.2, 0.3) ; names(GE) <- c("DETRITUS", "PLANT", "ANIMAL")
MR <- c(NA, 10/BM[2], 5/BM[3]) ; names(MR) <- c("DETRITUS", "PLANT", "ANIMAL")
result <- c(-1/BM[1] * (AE[2] * FM[1,2] + AE[3] * FM[1,3]),
            -MR[2],
            -MR[3]
            )

test_that("getDiagonal works with correct data input", {
  expect_equal(getDiagonal(MR = MR, FM = FM, BM = BM, AE = AE, dead = dead),
               result)
})

# Error: optional dead compartment is present, but not all other required data.
test_that("getDiagonal asks for all required data for dead compartments", {
  expect_error(getDiagonal(MR = MR, BM = BM, dead = dead),
               "please provide all required data to calculate dead diagonal values")
})
