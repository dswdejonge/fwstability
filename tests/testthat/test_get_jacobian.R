context("Jacobian matrix creation")

### No optional arguments

# Proper data format
FM <- matrix(c(0, 3, 2, 5, 0, 3, 0, 5, 0), nrow = 3, ncol = 3)
rownames(FM) <- c("PLANT", "WORM", "ANT")
colnames(FM) <- c("PLANT", "WORM", "ANT")
BM <- c(30, 20, 10) ; names(BM) <- c("PLANT", "WORM", "ANT")
AE <- c(0.1, 0.2, 0.3) ; names(AE) <- c("PLANT", "WORM", "ANT")
GE <- c(0.1, 0.2, 0.3) ; names(GE) <- c("PLANT", "WORM", "ANT")
JM <- matrix(c(0,
               AE[1] * GE[1] * FM[2,1] / BM[2] + -FM[1,2] / BM[2],
               AE[1] * GE[1] * FM[3,1] / BM[3],

               AE[2] * GE[2] * FM[1,2] / BM[1] + -FM[2,1] / BM[1],
               0,
               AE[2] * GE[2] * FM[3,2] / BM[3] + -FM[2,3] / BM[3],

               -FM[3,1] / BM[1],
               AE[3] * GE[3] * FM[2,3] / BM[2] + -FM[3,2] / BM[2],
               0
), nrow = 3, ncol = 3)
rownames(JM) <- c("PLANT", "WORM", "ANT")
colnames(JM) <- c("PLANT", "WORM", "ANT")

test_that("the function works without optional arguments", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE),
               JM)
})




### Dead compartment, no externals

# Proper data format
FM <- matrix(c(0, 3, 8, 7, 0, 0, 4, 4, 0), nrow = 3, ncol = 3)
rownames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(FM) <- c("DETRITUS", "PLANT", "ANIMAL")
BM <- c(30, 20, 10) ; names(BM) <- c("DETRITUS", "PLANT", "ANIMAL")
AE <- c(NA, 0.2, 0.3) ; names(AE) <- c("DETRITUS", "PLANT", "ANIMAL")
GE <- c(NA, 0.2, 0.3) ; names(GE) <- c("DETRITUS", "PLANT", "ANIMAL")
JM <- matrix(c(0,
             (FM[2,1] - FM[1,2] + FM[2,3]*(1-AE[3])) / BM[2],
             (FM[3,1] - FM[1,3]) / BM[3],
             AE[2] * GE[2] * FM[1,2] / BM[1],
             0,
             -FM[2,3] / BM[3],
             AE[3] * GE[3] * FM[1,3] / BM[1],
             AE[3] * GE[3] * FM[2,3] / BM[2],
             0
             ), nrow = 3, ncol = 3)
rownames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")

test_that("the function works with dead compartments", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE, dead = "DETRITUS"),
               JM)
})

# Proper data format, include vector diagonal
DIAG <- c(-1,-2,-3)
JM2 <- JM; JM2[c(1,5,9)] <- DIAG

test_that("the function incorporates vector diagonal correctly", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = "DETRITUS", diagonal = DIAG),
               JM2)
})

# Error: Diagonal vector has the wrong length
DIAG2 <- c(-1,-2,-3,-4)

test_that("the function only executes with correct diagonal length", {
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = "DETRITUS", diagonal = DIAG2),
               "given diagonal has incorrect length")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = "DETRITUS", diagonal = DIAG2[1:2]),
               "given diagonal has incorrect length")
})

# Error: Diagonal vector is not numeric
DIAG3 <- c("a", "b", "c")
test_that("the function only executes with a numeric diagonal", {
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = "DETRITUS", diagonal = DIAG3),
               "given diagonal not numeric")
})

# Error: Physiological values exist for dead compartments
AEe <- c(0.1, 0.2, 0.3) ; names(AEe) <- c("DETRITUS", "PLANT", "ANIMAL")
GEe <- c(0.1, 0.2, 0.3) ; names(GEe) <- c("DETRITUS", "PLANT", "ANIMAL")

test_that("dead compartments cannot have physiological values", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AEe, GE = GEe, dead = "DETRITUS"),
               JM)
  expect_warning(getJacobian(FM = FM, BM = BM, AE = AEe, GE = GEe, dead = "DETRITUS"),
                 "physiological values set to NA for dead compartments")
})

# Error: not all required data is provided
test_that("the function throws an error if not all required data is given", {
  expect_error(getJacobian(BM = BM, AE = AE, GE = GE, dead = "DETRITUS"),
               "argument \"FM\" is missing, with no default")
  expect_error(getJacobian(FM = FM, AE = AE, GE = GE, dead = "DETRITUS"),
               "argument \"BM\" is missing, with no default")
  expect_error(getJacobian(FM = FM, BM = BM, GE = GE, dead = "DETRITUS"),
               "argument \"AE\" is missing, with no default")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, dead = "DETRITUS"),
               "argument \"GE\" is missing, with no default")
})

# Error: the flowmatrix is not squared
test_that("the function only executes with a square matrix", {
  expect_error(getJacobian(FM[,1:2], BM, AE, GE, dead = "DETRITUS"),
               "flow matrix is not square")
})

# Error: the flow matrix and vectors are not named
FM2 <- FM; rownames(FM2) <- NULL; colnames(FM2) <- NULL
FM3 <- FM; rownames(FM3) <- NULL
FM4 <- FM; colnames(FM4) <- NULL
BM2 <- BM; names(BM2) <- NULL
AE2 <- AE; names(AE2) <- NULL
GE2 <- GE; names(GE2) <- NULL

test_that("the function only executes when all vectors and matrices are named", {
  expect_error(getJacobian(FM = FM2, BM = BM, AE = AE, GE = GE, dead = "DETRITUS"),
               "all required vectors and matrices must be named")
  expect_error(getJacobian(FM = FM3, BM = BM, AE = AE, GE = GE, dead = "DETRITUS"),
               "all required vectors and matrices must be named")
  expect_error(getJacobian(FM = FM4, BM = BM, AE = AE, GE = GE, dead = "DETRITUS"),
               "all required vectors and matrices must be named")
  expect_error(getJacobian(FM = FM, BM = BM2, AE = AE, GE = GE, dead = "DETRITUS"),
               "all required vectors and matrices must be named")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE2, GE = GE, dead = "DETRITUS"),
               "all required vectors and matrices must be named")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE2, dead = "DETRITUS"),
               "all required vectors and matrices must be named")
})

# Error: Names of some vectors are different or in wrong order
FM5 <- FM; rownames(FM5) <- c("A", "B", "C"); colnames(FM5) <- c("A", "B", "C")
FM6 <- FM; rownames(FM6) <- c("A", "B", "C")
FM7 <- FM; colnames(FM7) <- c("A", "B", "C")
BM4 <- BM; names(BM4) <- c("A", "B", "C")
BM5 <- BM; names(BM5) <- names(BM)[c(2,3,1)]
AE3 <- AE; names(AE3) <- c("A", "B", "C")
GE3 <- GE; names(GE3) <- c("A", "B", "C")

test_that("all vectors and matrices have the same names", {
  expect_error(getJacobian(FM = FM5, BM = BM, AE = AE, GE = GE, dead = "DETRITUS"),
               "the names and their order must be equal in all named vectors and matrices")
  expect_error(getJacobian(FM = FM6, BM = BM, AE = AE, GE = GE, dead = "DETRITUS"),
               "row names and column names of flow matrix do not match")
  expect_error(getJacobian(FM = FM7, BM = BM, AE = AE, GE = GE, dead = "DETRITUS"),
               "row names and column names of flow matrix do not match")
  expect_error(getJacobian(FM = FM, BM = BM4, AE = AE, GE = GE, dead = "DETRITUS"),
               "the names and their order must be equal in all named vectors and matrices")
  expect_error(getJacobian(FM = FM, BM = BM5, AE = AE, GE = GE, dead = "DETRITUS"),
               "the names and their order must be equal in all named vectors and matrices")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE3, GE = GE, dead = "DETRITUS"),
               "the names and their order must be equal in all named vectors and matrices")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE3, dead = "DETRITUS"),
               "the names and their order must be equal in all named vectors and matrices")
})

# Error: biomasses have incorrect values
BM6 <- BM; BM6[2] <- 0
BM7 <- BM; BM7[2] <- -10
BM8 <- BM; BM8[2] <- NA
BM9 <- BM; BM9[2] <- "foo"
test_that("the function only executes if all biomasses are greater than zero", {
  expect_error(getJacobian(FM, BM = BM6, AE, GE, dead = "DETRITUS"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobian(FM, BM = BM7, AE, GE, dead = "DETRITUS"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobian(FM, BM = BM8, AE, GE, dead = "DETRITUS"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobian(FM, BM = BM9, AE, GE, dead = "DETRITUS"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
})


### Externals
#test_that("the function works with external compartments", {
#  expect_equal(dim(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE, externals = "DETRITUS")),
#               c(2,2))
#})


### Both dead compartments and externals.
test_that("the function only executes when the dead and external compartments have a existing name", {
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE, dead = "CARCASS"),
               "the names of the dead compartments are unknown")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE, externals = "CARCASS"),
               "the names of the external compartments are unknown")
})


