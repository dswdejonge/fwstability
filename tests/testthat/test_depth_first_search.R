context("Test depth-first-search algorithms and wrapper.")

set.seed(1994)
AM <- matrix(c(
  0, 0, 1, 1, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 1, 0,
  1, 0, 0, 0, 1, 1, 0, 0,
  1, 1, 0, 0, 0, 1, 1, 0,
  1, 0, 1, 0, 0, 0, 0, 1,
  1, 0, 1, 1, 0, 0, 0, 1,
  0, 1, 0, 1, 0, 0, 0, 1,
  0, 1, 0, 0, 1, 1, 1, 0
), byrow = T, nrow = 8, ncol = 8)
AM1 <- AM ; AM1[3] <- 5

test_that("the function only executes with right data input", {
  expect_error(dfs(AM[,1:4]), "Adjacency matrix must be square.")
  expect_error(dfs(AM1), "Adjancency matrix can only contain 0 and 1.")
  expect_error(dfs(AM, k = "foo"), "k must be an integer.")
  expect_error(dfs(AM, output = 3), "output should contain string with filename for output.")
  expect_error(dfs())
})

# Manual test.
# dfs(AM, k = 3, output = "tests/testthat/allLoops")
#file <- "tests/testthat/allLoops_k=3.txt"
#test_that("All loops of certain length are found", {
#  expect_equal(
#    sum(read.delim(file,sep = " ")[,1:3]),
#    sum(1,1,2,2,2,3,4,4,7,7,6,6,7,4,8))
#})
