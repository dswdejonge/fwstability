context("Jacobian matrix creation with wrong data format")

# Proper data input
fwnames <- c("DETRITUS", "PLANT", "ANIMAL", "CO2")
FM <- matrix(c(0, 3, 8, 1, 7, 0, 0, 0, 4, 4, 0, 0, 0, 0, 1, 0), nrow = 4, ncol = 4)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
BM <- c(30, 20, 10) ; names(BM) <- fwnames[1:3]
AE <- c(NA, 0.2, 0.3) ; names(AE) <- fwnames[1:3]
GE <- c(NA, 0.2, 0.3) ; names(GE) <- fwnames[1:3]
dead <- list(names = "DETRITUS", def = "Def")
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
rownames(JM) <- fwnames[1:3]
colnames(JM) <- fwnames[1:3]


# Incomplete data
test_that("the function only executes when the dead and external compartments have an existing name", {
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = list(names = "CARCASS"), externals = "CO2"),
               "the names of the dead compartments are unknown")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = dead, externals = "CARCASS"),
               "the names of the external compartments are unknown")
})

# Error: Diagonal vector has the wrong length
DIAG2 <- c(-1,-2,-3,-4)

test_that("the function only executes with correct diagonal length", {
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = dead, diagonal = DIAG2, externals = "CO2"),
               "given diagonal has incorrect length")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = dead, diagonal = DIAG2[1:2], externals = "CO2"),
               "given diagonal has incorrect length")
})

# Error: Diagonal vector is not numeric
DIAG3 <- c("a", "b", "c")
test_that("the function only executes with a numeric diagonal", {
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = dead, diagonal = DIAG3, externals = "CO2"),
               "given diagonal not numeric")
})

# Error: Physiological values exist for dead compartments
AEe <- c(0.1, 0.2, 0.3) ; names(AEe) <- fwnames[1:3]
GEe <- c(0.1, 0.2, 0.3) ; names(GEe) <- fwnames[1:3]

test_that("dead compartments cannot have physiological values", {
  expect_equal(getJacobian(FM = FM, BM = BM, AE = AEe, GE = GEe,
                           dead = dead, externals = "CO2"),
               JM)
  expect_warning(getJacobian(FM = FM, BM = BM, AE = AEe, GE = GEe,
                             dead = dead, externals = "CO2"),
                 "physiological values set to NA for dead compartments")
})

# Error: not all required data is provided
test_that("the function throws an error if not all required data is given", {
  expect_error(getJacobian(BM = BM, AE = AE, GE = GE, dead = dead, externals = "CO2"),
               "argument \"FM\" is missing, with no default")
  expect_error(getJacobian(FM = FM, AE = AE, GE = GE, dead = dead, externals = "CO2"),
               "argument \"BM\" is missing, with no default")
  expect_error(getJacobian(FM = FM, BM = BM, GE = GE, dead = dead, externals = "CO2"),
               "argument \"AE\" is missing, with no default")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, dead = dead, externals = "CO2"),
               "argument \"GE\" is missing, with no default")
})

# Error: the flowmatrix is not squared
test_that("the function only executes with a square matrix", {
  expect_error(getJacobian(FM[,1:2], BM, AE, GE, dead = dead),
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
  expect_error(getJacobian(FM = FM2, BM = BM, AE = AE, GE = GE, dead = dead),
               "all required vectors and matrices must be named")
  expect_error(getJacobian(FM = FM3, BM = BM, AE = AE, GE = GE, dead = dead),
               "all required vectors and matrices must be named")
  expect_error(getJacobian(FM = FM4, BM = BM, AE = AE, GE = GE, dead = dead),
               "all required vectors and matrices must be named")
  expect_error(getJacobian(FM = FM, BM = BM2, AE = AE, GE = GE, dead = dead, externals = "CO2"),
               "all required vectors and matrices must be named")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE2, GE = GE, dead = dead, externals = "CO2"),
               "all required vectors and matrices must be named")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE2, dead = dead, externals = "CO2"),
               "all required vectors and matrices must be named")
})

# Error: Names of some vectors are different or in wrong order
FM5 <- FM; rownames(FM5) <- c("A", "B", "C", "D"); colnames(FM5) <- c("A", "B", "C", "D")
FM6 <- FM; rownames(FM6) <- c("A", "B", "C", "D")
FM7 <- FM; colnames(FM7) <- c("A", "B", "C", "D")
BM4 <- BM; names(BM4) <- c("A", "B", "C")
BM5 <- BM; names(BM5) <- names(BM)[c(2,3,1)]
AE3 <- AE; names(AE3) <- c("A", "B", "C")
GE3 <- GE; names(GE3) <- c("A", "B", "C")

test_that("all vectors and matrices have the same names", {
  expect_error(getJacobian(FM = FM5, BM = BM, AE = AE, GE = GE, dead = dead),
               "the names and their order must be equal in all named vectors and matrices")
  expect_error(getJacobian(FM = FM6, BM = BM, AE = AE, GE = GE, dead = dead),
               "row names and column names of flow matrix do not match")
  expect_error(getJacobian(FM = FM7, BM = BM, AE = AE, GE = GE, dead = dead),
               "row names and column names of flow matrix do not match")
  expect_error(getJacobian(FM = FM, BM = BM4, AE = AE, GE = GE, dead = dead, externals = "CO2"),
               "the names and their order must be equal in all named vectors and matrices")
  expect_error(getJacobian(FM = FM, BM = BM5, AE = AE, GE = GE, dead = dead, externals = "CO2"),
               "the names and their order must be equal in all named vectors and matrices")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE3, GE = GE, dead = dead, externals = "CO2"),
               "the names and their order must be equal in all named vectors and matrices")
  expect_error(getJacobian(FM = FM, BM = BM, AE = AE, GE = GE3, dead = dead, externals = "CO2"),
               "the names and their order must be equal in all named vectors and matrices")
})

# Error: biomasses have incorrect values
BM6 <- BM; BM6[2] <- 0
BM7 <- BM; BM7[2] <- -10
BM8 <- BM; BM8[2] <- NA
BM9 <- BM; BM9[2] <- "foo"
test_that("the function only executes if all biomasses are greater than zero", {
  expect_error(getJacobian(FM, BM = BM6, AE, GE, dead = dead, externals = "CO2"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobian(FM, BM = BM7, AE, GE, dead = dead, externals = "CO2"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobian(FM, BM = BM8, AE, GE, dead = dead, externals = "CO2"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobian(FM, BM = BM9, AE, GE, dead = dead, externals = "CO2"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
})
