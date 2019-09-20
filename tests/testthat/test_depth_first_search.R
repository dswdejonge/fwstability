context("Test depth-first-search algorithms and wrapper.")

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

#dfs(AM)
#dfs(AM, k = 3)
