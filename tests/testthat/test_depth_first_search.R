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

FM <- AM
FM[which(FM == 1)] <- runif(length(which(FM == 1)), min = 0.01, max = 10)
model <- list(
  type = "EF",
  FM = FM,
  BM = runif(dim(FM)[1], min = 5, max = 25),
  AE = runif(dim(FM)[1], min = 0, max = 1),
  GE = runif(dim(FM)[1], min = 0, max = 1),
  dead = list(names = c(1,2)),
  i = T
)



test_that("the function only executes with right data input", {
  expect_error(dfs(AM[,1:4]), "Adjacency matrix must be square.")
  expect_error(dfs(AM1), "Adjancency matrix can only contain 0 and 1.")
  expect_error(dfs(AM, k = "foo"), "k must be an integer.")
  expect_error(dfs(AM, output = 3), "output should contain string with filename for output.")
  expect_error(dfs())
})

#dfs(AM)
#dfs(AM, k = 3)
