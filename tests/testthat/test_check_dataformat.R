context("Jacobian matrix creation with wrong data format")

# Proper data input
fwnames <- c("DETRITUS", "PLANT", "ANIMAL", "CO2")
FM <- matrix(c(0, 3, 8, 1, 7, 0, 0, 0, 4, 4, 0, 0, 0, 0, 1, 0), nrow = 4, ncol = 4)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
BM <- c(30, 20, 10) ; names(BM) <- fwnames[1:3]
AE <- c(NA, 0.2, 0.3) ; names(AE) <- fwnames[1:3]
GE <- c(NA, 0.2, 0.3) ; names(GE) <- fwnames[1:3]
dead <- list(names = "DETRITUS", frac = FM[1:3, 1:3])
# Expected answer
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

# Test input
test_that("example is proper and order doesn't mater", {
  expect_equal(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE, GE = GE,
                                     dead = dead, externals = "CO2", diagonal = 0),
               JM)
  expect_equal(getJacobianEnergyFlux(FM = FM, BM = BM[c(3,2,1)], AE = AE, GE = GE,
                                     dead = dead, externals = "CO2", diagonal = 0),
               JM)
  expect_equal(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE[c(3,2,1)], GE = GE,
                                     dead = dead, externals = "CO2", diagonal = 0),
               JM)
  expect_equal(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE, GE = GE[c(3,2,1)],
                                     dead = dead, externals = "CO2", diagonal = 0),
               JM)
})

# Incomplete data
test_that("the function only executes when the dead and external compartments have an existing name", {
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = list(names = "CARCASS", frac = FM[1:3,1:3]), externals = "CO2"),
               "the names of the dead compartments are unknown")
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = dead, externals = "CARCASS"),
               "the names of the external compartments are unknown")
})

# Error: Diagonal vector has the wrong length
DIAG2 <- c(-1,-2,-3,-4)

test_that("the function only executes with correct diagonal length", {
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = dead, diagonal = DIAG2, externals = "CO2"),
               "given diagonal has incorrect length")
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = dead, diagonal = DIAG2[1:2], externals = "CO2"),
               "given diagonal has incorrect length")
})

# Error: Diagonal vector is not numeric
DIAG3 <- c("a", "b", "c")
test_that("the function only executes with a numeric diagonal", {
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE, GE = GE,
                           dead = dead, diagonal = DIAG3, externals = "CO2"),
               "given diagonal not numeric")
})

# Error: Physiological values exist for dead compartments
AEe <- c(0.1, 0.2, 0.3) ; names(AEe) <- fwnames[1:3]
GEe <- c(0.1, 0.2, 0.3) ; names(GEe) <- fwnames[1:3]

test_that("dead compartments cannot have physiological values", {
  expect_warning(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AEe, GE = GEe,
                             dead = dead, externals = "CO2"),
                 "physiological values set to NA for dead compartments")
})

# Error: not all required data is provided
test_that("the function throws an error if not all required data is given", {
  expect_error(getJacobianEnergyFlux(BM = BM, AE = AE, GE = GE, dead = dead, externals = "CO2"),
               "argument \"FM\" is missing, with no default")
  expect_error(getJacobianEnergyFlux(FM = FM, AE = AE, GE = GE, dead = dead, externals = "CO2"),
               "argument \"BM\" is missing, with no default")
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, GE = GE, dead = dead, externals = "CO2"),
               "argument \"AE\" is missing, with no default")
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE, dead = dead, externals = "CO2"),
               "argument \"GE\" is missing, with no default")
})

# Error: the flowmatrix is not squared
test_that("the function only executes with a square matrix", {
  expect_error(getJacobianEnergyFlux(FM[,1:2], BM, AE, GE, dead = dead),
               "Input matrix is not square")
})

# Error: the flow matrix and vectors are not named
FM2 <- FM; rownames(FM2) <- NULL; colnames(FM2) <- NULL
FM3 <- FM; rownames(FM3) <- NULL
FM4 <- FM; colnames(FM4) <- NULL
BM2 <- BM; names(BM2) <- NULL
AE2 <- AE; names(AE2) <- NULL
GE2 <- GE; names(GE2) <- NULL

test_that("the function only executes when all vectors and matrices are named", {
  expect_error(getJacobianEnergyFlux(FM = FM2, BM = BM, AE = AE, GE = GE, dead = dead),
               "Input matrix must have named rows and columns.")
  expect_error(getJacobianEnergyFlux(FM = FM3, BM = BM, AE = AE, GE = GE, dead = dead),
               "Input matrix must have named rows and columns.")
  expect_error(getJacobianEnergyFlux(FM = FM4, BM = BM, AE = AE, GE = GE, dead = dead),
               "Input matrix must have named rows and columns.")
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM2, AE = AE, GE = GE, dead = dead, externals = "CO2"),
               "All required vectors must be named.")
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE2, GE = GE, dead = dead, externals = "CO2"),
               "All required vectors must be named.")
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE, GE = GE2, dead = dead, externals = "CO2"),
               "All required vectors must be named.")
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
  # This test also gives an R base warning message
  expect_error(getJacobianEnergyFlux(FM = FM5, BM = BM, AE = AE, GE = GE, dead = dead),
               "The names must be equal in all named vectors and matrices.")
  expect_error(getJacobianEnergyFlux(FM = FM6, BM = BM, AE = AE, GE = GE, dead = dead),
               "Input matrix must have same names in rows and columns.")
  expect_error(getJacobianEnergyFlux(FM = FM7, BM = BM, AE = AE, GE = GE, dead = dead),
               "Input matrix must have same names in rows and columns.")
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM4, AE = AE, GE = GE, dead = dead, externals = "CO2"),
               "The names must be equal in all named vectors and matrices.")
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE3, GE = GE, dead = dead, externals = "CO2"),
               "The names must be equal in all named vectors and matrices.")
  expect_error(getJacobianEnergyFlux(FM = FM, BM = BM, AE = AE, GE = GE3, dead = dead, externals = "CO2"),
               "The names must be equal in all named vectors and matrices.")
})

# Error: biomasses have incorrect values
BM6 <- BM; BM6[2] <- 0
BM7 <- BM; BM7[2] <- -10
BM8 <- BM; BM8[2] <- NA
BM9 <- BM; BM9[2] <- "foo"
test_that("the function only executes if all biomasses are greater than zero", {
  expect_error(getJacobianEnergyFlux(FM, BM = BM6, AE, GE, dead = dead, externals = "CO2"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobianEnergyFlux(FM, BM = BM7, AE, GE, dead = dead, externals = "CO2"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobianEnergyFlux(FM, BM = BM8, AE, GE, dead = dead, externals = "CO2"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  expect_error(getJacobianEnergyFlux(FM, BM = BM9, AE, GE, dead = dead, externals = "CO2"),
               "biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
})

# Error: AE and GE have incorrect values
AE1 <- AE; AE1[3] <- -1
AE2 <- AE; AE2[3] <- 2
GE1 <- GE; GE1[3] <- -1
GE2 <- GE; GE2[3] <- 2
test_that("the function only executes if all AE and GE are 0-1.", {
  expect_error(getJacobianEnergyFlux(FM, BM = BM, AE1, GE, dead = dead, externals = "CO2"),
               "assimilation and growth efficiencies must lie between 0 and 1")
  expect_error(getJacobianEnergyFlux(FM, BM = BM, AE2, GE, dead = dead, externals = "CO2"),
               "assimilation and growth efficiencies must lie between 0 and 1")
  expect_error(getJacobianEnergyFlux(FM, BM = BM, AE, GE1, dead = dead, externals = "CO2"),
               "assimilation and growth efficiencies must lie between 0 and 1")
  expect_error(getJacobianEnergyFlux(FM, BM = BM, AE, GE2, dead = dead, externals = "CO2"),
               "assimilation and growth efficiencies must lie between 0 and 1")
})

# Error: dead is wrong data format
test_that("the list \"dead\" format gets checked correctly", {
  expect_error(checkDeadFormat("DET"),
               "argument \"dead\" must be a named list")
  expect_error(checkDeadFormat(list("DET")),
               "argument \"dead\" must be a named list")
  expect_error(checkDeadFormat(list(def = "noDef")),
               "\"names\" element is required in the \"dead\" list")
  expect_error(checkDeadFormat(list(names = "DET", def = "Def", frac = 0, extra = 0)),
               "the list \"dead\" should have 2 elements at most")
  expect_error(checkDeadFormat(list(names = c("DET", "NUT"))),
               "the \"frac\" element is required in the \"dead\" list")
})

### Stability
# ***********
# Example matrix
JM <- matrix(c(1,3,6,-3,-5,-6,3,3,4),
             nrow = 3, ncol = 3)
rownames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")
colnames(JM) <- c("DETRITUS", "PLANT", "ANIMAL")
MR <- c(1, 1, 1)
names(MR) <- c("DETRITUS", "PLANT", "ANIMAL")

# Error: the JM is not square, numeric, or named
JM1 <- JM[1:2,]
JM2 <- matrix(c("a","b","c","d"), nrow = 2, ncol = 2)
rownames(JM2) <- c("DETRITUS", "PLANT")
colnames(JM2) <- c("DETRITUS", "PLANT")
JM3 <- JM ; rownames(JM3) <- NULL
JM4 <- JM ; colnames(JM4) <- NULL
JM5 <- JM ; rownames(JM5) <- NULL ; colnames(JM5) <- NULL
JM6 <- JM ; rownames(JM6) <- rownames(JM6)[c(2,1,3)]
JM8 <- JM ; rownames(JM8)[1] <- "nonexistent"
test_that("the function only executes with a square, numeric, and correctly named JM", {
  expect_error(getStability(JM1),
               "Input matrix is not square")
  expect_error(getStability(JM2),
               "Input matrix must be numeric")
  expect_error(getStability(JM3),
               "Input matrix must have named rows and columns.")
  expect_error(getStability(JM4),
               "Input matrix must have named rows and columns.")
  expect_error(getStability(JM5),
               "Input matrix must have named rows and columns.")
  expect_error(getStability(JM6),
               "Input matrix must have same names in rows and columns.")
  expect_error(getStability(JM8),
               "Input matrix must have same names in rows and columns.")
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
               "NAs and NaNs not allowed in matrix")
})

# Warning: mortalities and/or dead compartments are given, which is irrelevant for the
# eigenvalue method
mortalities <- c(1, 1, 1) ; names(mortalities) <- c("DETRITUS", "PLANT", "ANIMAL")
test_that("a warning is produced when data irrelevant to the method is given", {
  expect_warning(getStability(JM, method = "eigenvalue", mortalities),
                 "given mortality values or dead compartments are irrelevant for the eigenvalue method")
  expect_warning(getStability(JM, method = "eigenvalue", dead = "DETRITUS"),
                 "given mortality values or dead compartments are irrelevant for the eigenvalue method")
  expect_warning(getStability(JM, method = "eigenvalue", mortalities, dead = "DETRITUS"),
                 "given mortality values or dead compartments are irrelevant for the eigenvalue method")
})

# Error: the MR vector is not numeric, named correctly, or has <=0 values
mort <- c("a", "b", "c") ; names(mort) <- names(MR)
mort2 <- MR[1:2]
mort3 <- MR ; names(mort3) <- NULL
mort4 <- MR ; names(mort4) <- names(MR)[c(2,3,1)]
mort5 <- MR ; names(mort5)[1] <- "nonexistent"
mort7 <- c(1,1,0) ; names(mort7) <- names(MR)
mort8 <- c(-1,-1,-1) ; names(mort8) <- names(MR)

test_that("the scalar method only works with a correct mortality vector", {
  expect_error(getStability(JM, method = "scalar", MR = mort),
               "the MR vector contains values equal or smaller than zero, or is non-numeric")
  expect_error(getStability(JM, method = "scalar", MR = mort7),
               "the MR vector contains values equal or smaller than zero, or is non-numeric")
  expect_error(getStability(JM, method = "scalar", MR = mort8),
               "the MR vector contains values equal or smaller than zero, or is non-numeric")
  # This test also gives an R base warning message
  expect_error(getStability(JM, method = "scalar", MR = mort2),
               "The names must be equal in all named vectors and matrices.")
  expect_error(getStability(JM, method = "scalar", MR = mort3),
               "All required vectors must be named.")
  expect_error(getStability(JM, method = "scalar", MR = mort5),
               "The names must be equal in all named vectors and matrices.")
})

# Error: the name of dead compartment is non-existent
test_that("the scalar method only works with an existing dead compartment", {
  expect_error(getStability(JM, method = "scalar", MR, dead = "foo"),
               "the names of the dead compartments are unknown")
})

# Error: the NA values in the mortality vactor are of non-dead compartment
mort6 <- MR ; mort6[2] <- NA
test_that("non-dead compartments cannot have NA MR with scalar method", {
  expect_error(getStability(JM, method = "scalar",
                            MR = mort6, dead = "DETRITUS"),
               "Mortality rates of non-dead compartments cannot be NA.")
})

# Error: the method is scalar, but there is no input
test_that("the scalar method asks for the correct input", {
  expect_error(getStability(JM, method = "scalar"),
               "MR vector required for the scalar method")
  expect_error(getStability(JM, method = "scalar", dead = "DETRITUS"),
               "MR vector required for the scalar method")
})

