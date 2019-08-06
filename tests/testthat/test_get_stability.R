context("Evaluating JM stability")

# Example matrix
# source: https://www.scss.tcd.ie/~dahyotr/CS1BA1/SolutionEigen.pdf
JM <- matrix(c(1,3,6,-3,-5,-6,3,3,4),
             nrow = 3, ncol = 3)
rownames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")



#### Method 1 - eigenvalue
# Correct answer
stability <- 4
test_that("the default getStability provides the correct answer", {
  expect_equal(getStability(JM, method = "eigenvalue"),
               stability)
  expect_equal(getStability(JM),
               stability)
})

# Error: the JM is not square, numeric, or named
JM1 <- JM[1:2,]
JM2 <- matrix(c("a","b","c","d"), nrow = 2, ncol = 2)
JM3 <- JM ; rownames(JM3) <- NULL
JM4 <- JM ; colnames(JM4) <- NULL
JM5 <- JM ; rownames(JM5) <- NULL ; colnames(JM5) <- NULL
JM6 <- JM ; rownames(JM6) <- rownames(JM6)[c(2,1,3)]
JM8 <- JM ; rownames(JM8)[1] <- "nonexistent"
test_that("the function only executes with a square, numeric, and correctly named JM", {
  expect_error(getStability(JM1),
               "Jacobian matrix is not square")
  expect_error(getStability(JM2),
               "Jacobian matrix must be numeric")
  expect_error(getStability(JM3),
               "all required vectors and matrices must be named")
  expect_error(getStability(JM4),
               "all required vectors and matrices must be named")
  expect_error(getStability(JM5),
               "all required vectors and matrices must be named")
  expect_error(getStability(JM6),
               "row names and column names of Jacobian matrix do not match")
  expect_error(getStability(JM8),
               "row names and column names of Jacobian matrix do not match")
})

# Error: the chosen method is unknown
test_that("the function only executes with available methods", {
  expect_error(getStability(JM, method = "nonexistent"),
               "unknown method chosen")
})

# Error: the diagonal contains NAs
JM7 <- JM ; diag(JM7)[1] <- NA
test_that("the eigenvalue method doesn't execute with diagonal NAs", {
  expect_error(getStability(JM7),
               "for the eigenvalue method the diagonal cannot contain NAs")
})

# Warning: mortalities and/or dead compartments are given, which is irrelevant for the
# eigenvalue method
test_that("a warning is produced when data irrelevant to the method is given", {
  expect_warning(getStability(JM, method = "eigenvalue", mortalities),
                "given mortality values or dead compartments are irrelevant for the eigenvalue method")
  expect_warning(getStability(JM, method = "eigenvalue", dead = "DETRITUS"),
                 "given mortality values or dead compartments are irrelevant for the eigenvalue method")
  expect_warning(getStability(JM, method = "eigenvalue", mortalities, dead = "DETRITUS"),
                 "given mortality values or dead compartments are irrelevant for the eigenvalue method")
  expect_equal(getStability(JM, method = "eigenvalue", mortalities, dead = "DETRITUS"),
                 stability)
})




#### Method 2 - scalar

# Correct answer - no dead compartments
mortalities <- c(1, 1, 1)
names(mortalities) <- c("DETRITUS", "PLANT", "ANIMAL")
stability <- 2.068195

test_that("the scalar getStability provides the correct answer", {
  expect_equal(getStability(JM, method = "scalar", mortalities = mortalities) > 2.068 &
               getStability(JM, method = "scalar", mortalities = mortalities) < 2.069,
               TRUE)
})

# Correct answer - including dead compartment
mortalitiesNA <- c(NA, 1, 1)
names(mortalitiesNA) <- c("DETRITUS", "PLANT", "ANIMAL")
stability <- 2.9

test_that("the scalar getStability provides the correct answer including detritus", {
  expect_equal(getStability(JM, method = "scalar",
                            mortalities = mortalitiesNA, dead = "DETRITUS"),
               stability)
})

# Error: the mortalities vector is not numeric, named correctly, or has <=0 values
mort <- c("a", "b", "c") ; names(mort) <- names(mortalities)
mort2 <- mortalities[1:2]
mort3 <- mortalities ; names(mort3) <- NULL
mort4 <- mortalities ; names(mort4) <- names(mortalities)[c(2,3,1)]
mort5 <- mortalities ; names(mort5)[1] <- "nonexistent"
mort7 <- c(1,1,0) ; names(mort7) <- names(mortalities)
mort8 <- c(-1,-1,-1) ; names(mort8) <- names(mortalities)

test_that("the scalar method only works with a correct mortality vector", {
  expect_error(getStability(JM, method = "scalar", mortalities = mort),
               "the mortalities vector contains values equal or smaller than zero, or is non-numeric")
  expect_error(getStability(JM, method = "scalar", mortalities = mort7),
               "the mortalities vector contains values equal or smaller than zero, or is non-numeric")
  expect_error(getStability(JM, method = "scalar", mortalities = mort8),
               "the mortalities vector contains values equal or smaller than zero, or is non-numeric")
  expect_error(getStability(JM, method = "scalar", mortalities = mort2),
               "the names and their order must be equal in all named vectors and matrices")
  expect_error(getStability(JM, method = "scalar", mortalities = mort3),
               "all required vectors and matrices must be named")
  expect_error(getStability(JM, method = "scalar", mortalities = mort4),
               "the names and their order must be equal in all named vectors and matrices")
  expect_error(getStability(JM, method = "scalar", mortalities = mort5),
               "the names and their order must be equal in all named vectors and matrices")
})

# Error: the name of dead compartment is non-existent
test_that("the scalar method only works with an existing dead compartment", {
  expect_error(getStability(JM, method = "scalar", mortalities, dead = "foo"),
               "the names of the dead compartments are unknown")
})

# Error: the NA values in the mortality vactor are of non-dead compartment
mort6 <- mortalities ; mort6[2] <- NA
test_that("non-dead compartments cannot have NA mortalities with scalar method", {
  expect_error(getStability(JM, method = "scalar",
                            mortalities = mort6, dead = "DETRITUS"),
               "mortalities values are set to NA for non-dead compartments")
})

# Error: the method is scalar, but there is no input
test_that("the scalar method asks for the correct input", {
  expect_error(getStability(JM, method = "scalar"),
               "mortalities vector required for the scalar method")
  expect_error(getStability(JM, method = "scalar", dead = "DETRITUS"),
               "mortalities vector required for the scalar method")
})

