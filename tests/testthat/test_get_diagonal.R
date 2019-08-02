context("Diagonal value calculation")

#### getSpeciesDiagonal

# Proper data
MR <- c(10,20,30)
BM <- c(3,2,1)

test_that("getDiagonalSpecies works with correct input data", {
  expect_equal(getDiagonalSpecies(MR = MR, BM = BM),
               -MR/BM)
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
