context("Adjust the list \"dead\" input")


answer1 <- list(names = "DET", def = NULL, frac = NULL)
answer2 <- list(names = c("DET","NUT"), def = NULL, frac = NULL)
answer3 <- list(names = c("DET", "NUT"), def = c("Def", "noDef"), frac = NULL)

test_that("the list \"dead\" gets adjusted correctly", {
  expect_equal(adjustDeadInput(list(names = "DET")), answer1)
  expect_equal(adjustDeadInput(list(names = c("DET","NUT"))), answer2)
  expect_equal(adjustDeadInput(list(names = c("DET", "NUT"), def = c("Def", "noDef"))), answer3)
  expect_error(adjustDeadInput("DET"),
               "argument \"dead\" must be a named list")
  expect_error(adjustDeadInput(list("DET")),
               "argument \"dead\" must be a named list")
  expect_error(adjustDeadInput(list(def = "noDef")),
               "\"names\" element is required in the \"dead\" list")
  expect_error(adjustDeadInput(list(names = "DET", def = "Def", frac = 0, extra = 0)),
               "the list \"dead\" should have 3 elements at most")
  expect_error(adjustDeadInput(list(names = c("DET", "NUT"), def = c("foo", "noDef"))),
               "the second element of the list \"dead\" may only contain the strings \"Def\" and \"noDef\"")
})

