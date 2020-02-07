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

#dfs(AM)
#dfs(AM, k = 3)

fwnames <- LETTERS[1:dim(AM)[1]]
FM <- AM
FM[which(FM == 1)] <- runif(length(which(FM == 1)), min = 0.01, max = 10)
BM <- runif(dim(FM)[1], min = 5, max = 25)
AE <- runif(dim(FM)[1], min = 0, max = 1)
GE <- runif(dim(FM)[1], min = 0, max = 1)
MR <- runif(dim(FM)[1], min = 0, max = 2)
rownames(FM) <- fwnames ; colnames(FM) <- fwnames
names(BM) <- fwnames ; names(AE) <- fwnames ; names(GE) <- fwnames
model <- list(
  type = "EF",
  FM = FM,
  BM = BM,
  AE = AE,
  GE = GE,
  dead = list(names = "A", frac = FM)
)
JM <- getJacobian(model)

test_that("the function only executes with right data input", {
  expect_error(dfs(AM[,1:4]), "Adjacency matrix must be square.")
  expect_error(dfs(AM1), "Adjancency matrix can only contain 0 and 1.")
  expect_error(dfs(AM, k = "foo"), "k must be an integer.")
  expect_error(dfs(AM, output = 3), "output should contain string with filename for output.")
  expect_error(dfs())
})
k <- 3
IS <- JM[which(JM != 0)][1:k]

test_that("getFeedback works", {
  expect_equal(getFeedback(IS), prod(IS))
})

test_that("getLoopWeight works", {
  expect_equal(getLoopWeight(IS), abs(prod(IS)) ^ (1 / k))
  expect_equal(getLoopWeight(IS, MR), abs(prod(IS) / prod(MR)) ^ (1 / k))
})

