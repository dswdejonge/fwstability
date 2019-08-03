context("Evaluating JM stability")

JM <- matrix(c(0, -0.25, 0.002, -0.09, 0, -0.49, -0.07, -0.13, 0),
             nrow = 3, ncol = 3)

test_that("the getStability function works", {
  expect_equal(is.numeric(getStability(JM)),
               TRUE)
})
