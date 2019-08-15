context("Adjust the list \"dead\" input")


answer1 <- list(names = "DET", def = NULL, frac = NULL)
answer2 <- list(names = c("DET","NUT"), def = NULL, frac = NULL)
answer3 <- list(names = c("DET", "NUT"), def = c("Def", "noDef"), frac = NULL)
answer4 <- list(names = c("DET", "NUT"), frac = c(0.5, 0.5), def = NULL)

test_that("the list \"dead\" gets adjusted correctly", {
  expect_equal(adjustDeadInput("DET"), answer1)
  expect_equal(adjustDeadInput(list("DET")), answer1)
  expect_equal(adjustDeadInput(c("DET","NUT")), answer2)
  expect_equal(adjustDeadInput(list(c("DET", "NUT"))), answer2)
  expect_equal(adjustDeadInput(list(c("DET", "NUT"), c("Def", "noDef"))), answer3)
  expect_equal(adjustDeadInput(list(c("DET", "NUT"), c(0.5, 0.5))), answer4)
})

