context("Assessing stabilizing characteristics")
require(ineq)

# Example data and answers
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

# Tests
test_that("The getWConn functions provide correct answers.", {
  expect_equal(fluxSizeDiversity(FM), H)
  expect_equal(averageMutualInfo(FM), A)
  expect_equal(getCw(FM), Cw)
  expect_equal(Gini(as.vector(FM[which(FM > 0)])), gini)
})

s_ini <- max(Re(eigen(JM)$values))
s_new <- max(Re(eigen(JM[-1, -1])$values))
test_that("assessComps provides correct dataframe", {
  expect_equal(ncol(assessComps(JM)), 5)
  expect_equal(nrow(assessComps(JM)), dim(JM)[1])
  expect_equal(as.character(assessComps(JM)[1,1]), rownames(JM)[1])
  expect_equal(assessComps(JM)[1,2], s_ini)
  expect_equal(assessComps(JM)[1,3], s_new)
  expect_equal(assessComps(JM)[1,4], s_new - s_ini)
  expect_equal(assessComps(JM)[1,5], (s_new - s_ini)/s_ini)

})

s_ini <- max(Re(eigen(JM)$values))
JM2 <- JM ; JM2[1] <- JM[1] * 2
s_new2 <- max(Re(eigen(JM2)$values))
test_fun <- function(x){return(x / 10)}
JM3 <- JM ; JM3[1] <- test_fun(JM[1])
s_new3 <- max(Re(eigen(JM3)$values))
test_that("assessLinksFixed provides correct dataframe", {
  expect_equal(ncol(assessLinksFixed(JM)), 6)
  expect_equal(nrow(assessLinksFixed(JM)), dim(JM)[1]*dim(JM)[2])
  expect_equal(as.character(unique(assessLinksFixed(JM)$eff.of)), rownames(JM))
  expect_equal(assessLinksFixed(JM)[1,3], s_ini)
  expect_equal(assessLinksFixed(JM)[1,4], s_new2)
  expect_equal(assessLinksFixed(JM)[1,5], s_new2 - s_ini)
  expect_equal(assessLinksFixed(JM)[1,6], (s_new2 - s_ini)/s_ini)
  expect_equal(assessLinksFixed(JM, func = test_fun)[1,3], s_ini)
  expect_equal(assessLinksFixed(JM, func = test_fun)[1,4], s_new3)
  expect_equal(assessLinksFixed(JM, func = test_fun)[1,5], s_new3 - s_ini)
  expect_equal(assessLinksFixed(JM, func = test_fun)[1,6], (s_new3 - s_ini)/s_ini)
})

test_that("assessLinksPerm provides correct dataframe", {
  expect_equal(ncol(assessLinksPerm(JM)), 3)
  expect_equal(nrow(assessLinksPerm(JM)),
               factorial(dim(JM)[1])/(factorial(2)*factorial(dim(JM)[1]-2)))
  expect_equal(all(assessLinksPerm(JM)$Pinstabiliy < 1 & assessLinksPerm(JM)$Pinstabiliy > 0), TRUE)
})

m <- matrix(1:9, 3, 3)
H <- -(
  ((m[1]/sum(m)) * log(m[1]/sum(m))) +
  ((m[2]/sum(m)) * log(m[2]/sum(m))) +
  ((m[3]/sum(m)) * log(m[3]/sum(m))) +
  ((m[4]/sum(m)) * log(m[4]/sum(m))) +
  ((m[5]/sum(m)) * log(m[5]/sum(m))) +
  ((m[6]/sum(m)) * log(m[6]/sum(m))) +
  ((m[7]/sum(m)) * log(m[7]/sum(m))) +
  ((m[8]/sum(m)) * log(m[8]/sum(m))) +
  ((m[9]/sum(m)) * log(m[9]/sum(m)))
  )

A <- ((m[1]/sum(m)) * log((m[1]*sum(m))/(sum(m[1,])*sum(m[,1])))) +
  ((m[2]/sum(m)) * log((m[2]*sum(m))/(sum(m[2,])*sum(m[,1])))) +
  ((m[3]/sum(m)) * log((m[3]*sum(m))/(sum(m[3,])*sum(m[,1])))) +
  ((m[4]/sum(m)) * log((m[4]*sum(m))/(sum(m[1,])*sum(m[,2])))) +
  ((m[5]/sum(m)) * log((m[5]*sum(m))/(sum(m[2,])*sum(m[,2])))) +
  ((m[6]/sum(m)) * log((m[6]*sum(m))/(sum(m[3,])*sum(m[,2])))) +
  ((m[7]/sum(m)) * log((m[7]*sum(m))/(sum(m[1,])*sum(m[,3])))) +
  ((m[8]/sum(m)) * log((m[8]*sum(m))/(sum(m[2,])*sum(m[,3])))) +
  ((m[9]/sum(m)) * log((m[9]*sum(m))/(sum(m[3,])*sum(m[,3]))))
test_that("The fluxSizeDiversity provides the correct answer", {
  expect_equal(fluxSizeDiversity(m), H)
})
test_that("The averageMutualInfo provides the correct answer", {
  expect_equal(averageMutualInfo(m), A)
})
test_that("The getConnPerNode provides the correct answer", {
  expect_equal(getConnPerNode(m), exp((H-A)/2))
})
test_that("The getCw provides the correct answer", {
  expect_equal(getCw(m), (exp((H-A)/2))/dim(m)[1])
})

loop <- c(-1, 2, 3)
d <- c(0.5, 0.2, 0.1)
d2 <- d[1:2]
d3 <- c(-0.5, 0.2, 0.1)
d4 <- c(NA, 0.2, 0.1)
test_that("The getLoopWeight provides the correct answer", {
  expect_equal(getLoopWeight(loop), abs(prod(loop))^(1/length(loop)))
  expect_equal(getLoopWeight(loop, d), abs(prod(loop)/prod(d))^(1/length(loop)))
  expect_error(getLoopWeight(loop, d2), "Vector IS and vector d should be same length.")
  expect_error(getLoopWeight(loop, d3), "Death rates should all be positive.")
  expect_equal(getLoopWeight(loop, d4), abs(prod(loop)/prod(d4, na.rm = T))^(1/length(loop)))
})
test_that("The getFeedback provides the correct answer", {
  expect_equal(getFeedback(loop), loop[1]*loop[2]*loop[3])
})

test_that("The total number of possible loops is calculated correctly", {
  expect_equal(maxNrLoops(4), 60)
  expect_equal(maxNrLoops(5), 320)
  expect_equal(maxNrLoops(5, 5), 120)
  expect_equal(maxNrLoops(5, 2), 20)
})

# assessFeedback(JM, findLoops = T, k = 3)
