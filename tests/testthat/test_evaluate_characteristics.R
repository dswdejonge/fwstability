context("Assessing stabilizing characteristics")
require(ineq)

# Example matrix
fwnames <- c("DET", "BAC", "FUN", "ENCH")
FM <- matrix(c(
  0, 1606, 40.8, 11.1,
  0, 0,    0,    10.1,
  0, 0,    0,    0.09,
  0, 0,    0,     0),
  nrow = 4, ncol = 4, byrow = T)
rownames(FM) <- fwnames ; colnames(FM) <- fwnames
BM <- c(2500, 227.5, 2.13, 0.43) ; names(BM) <- fwnames
AE <- c(NA, 1.0, 1.0, 0.25) ; names(AE) <- fwnames
GE <- c(NA, 0.30, 0.30, 0.40) ; names(GE) <- fwnames
JM <- matrix(c(0,
               - FM[1,2] / BM[2],
               - FM[1,3] / BM[3],
               - FM[1,4] / BM[4],

               AE[2] * GE[2] * FM[1,2] / BM[1],
               0,
               0,
               - FM[2,4] / BM[4],

               AE[3] * GE[3] * FM[1,3] / BM[1],
               0,
               0,
               - FM[3,4] / BM[4],

               AE[4] * GE[4] * FM[1,4] / BM[1],
               AE[4] * GE[4] * FM[2,4] / BM[2],
               AE[4] * GE[4] * FM[3,4] / BM[3],
               0
), nrow = 4, ncol = 4)
rownames(JM) <- fwnames ; colnames(JM) <- fwnames
Fsum <- sum(FM)
Fm <- mean(FM[which(FM > 0)])

H <- -(
    FM[1,2]/Fsum*log(FM[1,2]/Fsum) +
    FM[1,3]/Fsum*log(FM[1,3]/Fsum) +
    FM[1,4]/Fsum*log(FM[1,4]/Fsum) +
    FM[2,4]/Fsum*log(FM[2,4]/Fsum) +
    FM[3,4]/Fsum*log(FM[3,4]/Fsum))

A <- FM[1,2]/Fsum*log((FM[1,2]*Fsum)/(sum(FM[1,])*sum(FM[,2]))) +
      FM[1,3]/Fsum*log((FM[1,3]*Fsum)/(sum(FM[1,])*sum(FM[,3]))) +
      FM[1,4]/Fsum*log((FM[1,4]*Fsum)/(sum(FM[1,])*sum(FM[,4]))) +
      FM[2,4]/Fsum*log((FM[2,4]*Fsum)/(sum(FM[2,])*sum(FM[,4]))) +
      FM[3,4]/Fsum*log((FM[3,4]*Fsum)/(sum(FM[3,])*sum(FM[,4])))

Cw <- exp((H-A)/2)/length(fwnames)

X <- 2*((length(which(FM > 0)))^2)*mean(FM[which(FM > 0)])
gini <- (abs(FM[1,2] - FM[1,3]) +
  abs(FM[1,2] - FM[1,4]) +
  abs(FM[1,2] - FM[2,4]) +
  abs(FM[1,2] - FM[3,4]) +

  abs(FM[1,3] - FM[1,2]) +
  abs(FM[1,3] - FM[1,4]) +
  abs(FM[1,3] - FM[2,4]) +
  abs(FM[1,3] - FM[3,4]) +

  abs(FM[1,4] - FM[1,2]) +
  abs(FM[1,4] - FM[1,3]) +
  abs(FM[1,4] - FM[2,4]) +
  abs(FM[1,4] - FM[3,4]) +

  abs(FM[2,4] - FM[1,2]) +
  abs(FM[2,4] - FM[1,3]) +
  abs(FM[2,4] - FM[1,4]) +
  abs(FM[2,4] - FM[3,4]) +

  abs(FM[3,4] - FM[1,2]) +
  abs(FM[3,4] - FM[1,3]) +
  abs(FM[3,4] - FM[1,4]) +
  abs(FM[3,4] - FM[2,4]))/X

test_that("The getWConn functions provide correct answers.", {
  expect_equal(fluxSizeDiversity(FM), H)
  expect_equal(averageMutualInfo(FM), A)
  expect_equal(getCw(FM), Cw)
  expect_equal(Gini(as.vector(FM[which(FM > 0)])), gini)
})

FM[which(FM > 0)]

test_that("The total number of possible loops is calculated correctly", {
  expect_equal(maxNrLoops(4), 60)
  expect_equal(maxNrLoops(5), 320)
  expect_equal(maxNrLoops(5, 5), 120)
  expect_equal(maxNrLoops(5, 2), 20)
})
