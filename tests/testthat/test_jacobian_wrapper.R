context("getJacobian wrapper function")

library(LIM)
readLIM <- Read(system.file("extdata", "foodweb2.lim", package = "fwstability"))
model <- list(
  type = "LIM",
  LIM = readLIM
)

test_that("getJacobian wrapper gives error with unknown model", {
  expect_error(getJacobian(list(type = "foo")),
               "Unknown model input")
  expect_message(getJacobian(model),
                 "No model solutions given, LIM resolved by minimizing sum of squares.")
})
