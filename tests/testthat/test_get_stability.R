context("Evaluating JM stability")

# Example matrix
# source: https://www.scss.tcd.ie/~dahyotr/CS1BA1/SolutionEigen.pdf
JM <- matrix(c(1,3,6,-3,-5,-6,3,3,4),
             nrow = 3, ncol = 3)
rownames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")


#### Method 1 - eigenvalue
# Correct answer
stability <- 4
test_that("the default getStability provides the correct answer", {
  expect_equal(getStability(JM, method = "eigenvalue"),
               stability)
  expect_equal(getStability(JM),
               stability)
  expect_equal(getStability(t(JM)),
               stability)
})


#### Method 2 - scalar
MR <- c(1, 1, 1)
names(MR) <- c("DETRITUS", "PLANT", "ANIMAL")
stability <- 2.068195

test_that("the scalar getStability provides the correct answer", {
  expect_equal(getStability(JM, method = "scalar", MR = MR) > 2.06 &
               getStability(JM, method = "scalar", MR = MR) < 2.07,
               TRUE)
  expect_equal(getStability(t(JM), method = "scalar", MR = MR) > 2.06 &
                 getStability(t(JM), method = "scalar", MR = MR) < 2.07,
               TRUE)
})

# Correct answer - including dead compartment
MRNA <- c(NA, 1, 1)
names(MRNA) <- c("DETRITUS", "PLANT", "ANIMAL")
stability <- 2.9

test_that("the scalar getStability provides the correct answer including detritus", {
  expect_equal(getStability(JM, method = "scalar",
                            MR = MRNA, dead = "DETRITUS"),
               stability)
  expect_equal(getStability(t(JM), method = "scalar",
                            MR = MRNA, dead = "DETRITUS"),
               stability)
})

# Method 3: Initial
stability <- 0.5*max(Re(eigen(JM + t(JM))$values))
test_that("the initial getStability provides the correct answer", {
  expect_equal(getStability(JM, method = "initial"),
               stability)
  expect_equal(getStability(t(JM), method = "initial"),
               stability)
})
