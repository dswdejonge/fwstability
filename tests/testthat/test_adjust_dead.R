context("Adjust the list \"dead\" input")


answer1 <- list(names = "DET", def = NULL)
answer2 <- list(names = c("DET","NUT"), def = NULL)
answer3 <- list(names = c("DET", "NUT"), def = c("Def", "noDef"))

test_that("the list \"dead\" gets adjusted correctly", {
  expect_equal(adjustDeadInput("DET"), answer1)
  expect_equal(adjustDeadInput(list("DET")), answer1)
  expect_equal(adjustDeadInput(c("DET","NUT")), answer2)
  expect_equal(adjustDeadInput(list(c("DET", "NUT"))), answer2)
  expect_equal(adjustDeadInput(list(c("DET", "NUT"), c("Def", "noDef"))), answer3)
})

