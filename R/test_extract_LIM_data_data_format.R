context("How the extractLIMdata function handles different data formats.")

# getTag
tag <- "tag"
mynumbers <- c(1, 2, 3)
mytext <- c("a", "b", "c")
mylist <- list(1, 2, 3)
mynamesSingle <- c("tagA", "B", "C")
mynamesDouble <- c("tagA", "B", "tagC")

vars1 <- mynumbers ; names(vars1) <- mynamesSingle ; res1 <- vars1[1] ; names(res1) <- "A"
vars2 <- mynumbers ; names(vars2) <- mynamesDouble ; res2 <- vars2[c(1,3)] ; names(res2) <- c("A", "C")
vars3 <- mytext ; names(vars3) <- mynamesSingle ; res3 <- vars3[1] ; names(res3) <- "A"
vars4 <- mytext ; names(vars4) <- mynamesDouble ; res4 <- vars4[c(1,3)] ; names(res4) <- c("A", "C")

test_that("the getTag function extract the right information", {
  expect_equal(getTag(vars = vars1, tag = tag), res1)
  expect_equal(getTag(vars = vars2, tag = tag), res2)
  expect_equal(getTag(vars = vars3, tag = tag), res3)
  expect_equal(getTag(vars = vars4, tag = tag), res4)
  expect_error(getTag(vars = vars1),
               "argument \"tag\" is missing, with no default")
  expect_error(getTag(vars = mylist, tag = tag),
               "getTag only accepts named vectors for argument \"vars\"")
  expect_error(getTag(vars = mynumbers, tag = tag),
               "getTag only accepts named vectors for argument \"vars\"")
  expect_error(getTag(vars = vars1, tag = 1),
               "getTag only accepts a string for argument \"tag\"")
  expect_error(getTag(vars = vars1, tag = c("tag", "B")),
               "getTag only accepts a string for argument \"tag\"")
})
