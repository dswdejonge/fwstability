context("Evaluating JM stability")

# source: https://www.scss.tcd.ie/~dahyotr/CS1BA1/SolutionEigen.pdf
JM <- matrix(c(1,3,6,-3,-5,-6,3,3,4),
             nrow = 3, ncol = 3)
rownames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")
#### Method 1 - eigenvalue
stability <- 4

test_that("the getStability function works by finding max real part eigenvalues", {
  expect_equal(getStability(JM, method = "eigenvalue"),
               stability)
})

#### Method 2 - scalar

# no dead compartments
mortalities <- c(1, 1, 1)
names(mortalities) <- c("DETRITUS", "PLANT", "ANIMAL")
stability <- 2.068195

test_that("the getStability function works by finding mortality scalar", {
  expect_equal(getStability(JM, method = "scalar", mortalities = mortalities) > 2.06 &
               getStability(JM, method = "scalar", mortalities = mortalities) < 2.07,
               TRUE)
})

# including dead compartment
mortalities <- c(NA, 1, 1)
names(mortalities) <- c("DETRITUS", "PLANT", "ANIMAL")
stability <- 2.9

test_that("the getStability function works by finding mortality scalar including detritus", {
  expect_equal(getStability(JM, method = "scalar",
                            mortalities = mortalities, dead = "DETRITUS"),
               stability)
})
