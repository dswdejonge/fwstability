context("Diagonal value calculation")

#### getSpeciesDiagonal

# Proper data
MR <- c(15,10,5) ; names(MR) <- c("PLANT", "WORM", "ANT")
BM <- c(30, 20, 10) ; names(BM) <- c("PLANT", "WORM", "ANT")
result <- -MR/BM ; names(result) <- c("PLANT", "WORM", "ANT")

test_that("getDiagonalSpecies works with correct input data", {
  expect_equal(getDiagonalSpecies(MR = MR, BM = BM),
               result)
})

# Error: the vector lengths are unequal
test_that("getDiagonalSpecies needs equal input vector lengths", {
  expect_error(getDiagonalSpecies(MR = MR, BM = BM[1:2]),
               "input vectors have unequal lengths")
  expect_error(getDiagonalSpecies(MR = MR[1:2], BM = BM),
               "input vectors have unequal lengths")
})

# Error: the vectors are not numeric
test_that("getDiagonalSpecies needs numeric input vectors", {
  expect_error(getDiagonalSpecies(MR = MR, BM = c("a","b","c")),
               "input vectors must be numeric")
  expect_error(getDiagonalSpecies(MR = c("a", "b", "c"), BM = BM),
               "input vectors must be numeric")
})

# Error: the vectors are unnamed
MR1 <- MR ; names(MR1) <- NULL
BM1 <- BM ; names(BM1) <- NULL
test_that("getDiagonalSpecies only executes with named vectors", {
  expect_error(getDiagonalSpecies(MR = MR1, BM = BM),
               "input vectors must be named")
  expect_error(getDiagonalSpecies(MR = MR, BM = BM1),
               "input vectors must be named")
})

# Error: the names do not match or in the wrong order
MR2 <- MR ; names(MR2) <- names(MR)[c(2,3,1)]
BM2 <- BM; names(BM2) <- c("a", "b", "c")
test_that("getDiagonalSpecies only executes with matching names of vectors", {
  expect_error(getDiagonalSpecies(MR = MR2, BM = BM),
               "names of vectors do not match")
  expect_error(getDiagonalSpecies(MR = MR, BM = BM2),
               "names of vectors do not match")
})


#### getDetritusDiagonal

# Proper data
FM <- matrix(c(0, 3, 8, 7, 0, 0, 4, 4, 0), nrow = 3, ncol = 3)
rownames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
BM <- c(30, 20, 10) ; names(BM) <- c("DETRITUS", "PLANT", "ANIMAL")
AE <- c(NA, 0.2, 0.3) ; names(AE) <- c("DETRITUS", "PLANT", "ANIMAL")
GE <- c(NA, 0.2, 0.3) ; names(GE) <- c("DETRITUS", "PLANT", "ANIMAL")
result <- -1/BM[1] * (AE[2] * FM[1,2] + AE[3] * FM[1,3])

test_that("getDiagonalDetritus works with correct data input", {
  expect_equal(getDiagonalDetritus(FM = FM, BM = BM, AE = AE, dead = "DETRITUS"),
               result)
})



#### getDiagonal
# Proper data
MR <- c(NA, 10, 5) ; names(MR) <- c("DETRITUS", "PLANT", "ANIMAL")
FM <- matrix(c(0, 3, 8, 7, 0, 0, 4, 4, 0), nrow = 3, ncol = 3)
rownames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
BM <- c(30, 20, 10) ; names(BM) <- c("DETRITUS", "PLANT", "ANIMAL")
AE <- c(NA, 0.2, 0.3) ; names(AE) <- c("DETRITUS", "PLANT", "ANIMAL")
GE <- c(NA, 0.2, 0.3) ; names(GE) <- c("DETRITUS", "PLANT", "ANIMAL")
result <- c(-1/BM[1] * (AE[2] * FM[1,2] + AE[3] * FM[1,3]),
            -MR[2] / BM[2],
            -MR[3] / BM[3]
            )

test_that("getDiagonal works with correct data input", {
  expect_equal(getDiagonal(MR = MR, FM = FM, BM = BM, AE = AE, dead = "DETRITUS"),
               result)
})
