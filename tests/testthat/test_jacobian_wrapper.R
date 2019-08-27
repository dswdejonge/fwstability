context("getJacobian wrapper function")

test_that("getJacobian wrapper gives error with unknown model", {
  expect_error(getJacobian(list(type = "foo")),
               "Unknown model input")
})
