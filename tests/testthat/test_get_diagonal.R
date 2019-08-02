context("Diagonal value calculation")

# Proper data
MR <- c(10,20,30)
BM <- c(3,2,1)

test_that("getDiagonalSpecies works with correct input data", {
  expect_equal(getDiagonalSpecies(MR = MR, BM = BM),
               -MR/BM)
})
