# Automatically include lim model into the getJacobian matrix.
# Format needs to be so that the function can automatically find:
# - flow matrix
# - biomasses
# - assimilation efficiencies
# - growth efficiencies
# - mortality rates
# - dead compartments
# - defecation flows
# - parallel flows to defecation into the same compartment
# - externals
# - getXratioMatrix

library(LIM)
# readLIM allows you to easily access all data from the input file
readLIM <- Read(system.file("extdata", "foodweb2.lim", package = "fwstability"))
# Setup can directly use the input file, or the readLIM variable.
# It transforms the model into a set of linear equations to be solved.
lim <- Setup(readLIM)
# The Least Distance optimization is used to quantify all unknown flows.
lim_solved <- Ldei(lim)

# Expected answers
fwnames <- c("LABILE", "REFRAC", "MEIO", "MACRO", "IN", "OUT")
nC = lim$NExternal + lim$NComponents
FMe <- matrix(c(
0,      0, 12.25, 0, 0, 0,
0,      0,     0, 0, 0, 1.0375,
4.9625, 1.0375,0, 5, 0, 1.25,
4.4,    0,     0, 0, 0, 0.6,
2.8875, 0,     0, 0, 0, 0,
0,      0,     0, 0, 0, 0
), nrow = nC, ncol = nC, byrow = TRUE)
rownames(FMe) <- fwnames ; colnames(FMe) <- fwnames
variables <- c(lim_solved$X["meioDefLab"] + lim_solved$X["meioDefRefrac"],
               lim_solved$X["macroDefLab"] + lim_solved$X["macroDefRefrac"],
               lim_solved$X["meioGrazDet"] - lim_solved$X["meioDefLab"] - lim_solved$X["meioDefRefrac"],
               lim_solved$X["macroPredMeio"] - lim_solved$X["macroDefLab"] - lim_solved$X["macroDefRefrac"],
               lim_solved$X["meioGrazDet"] - lim_solved$X["meioDefLab"] - lim_solved$X["meioDefRefrac"] - lim_solved$X["meioResp"],
               lim_solved$X["macroPredMeio"] - lim_solved$X["macroDefLab"] - lim_solved$X["macroDefRefrac"] - lim_solved$X["macroResp"]
               )
names(variables) <- lim$Variables

# Test getFlowMatrix function
test_that("the Flowmatrix function works with parallel flows", {
  #expect_equal(Flowmatrix(lim, web = lim_solved$X), FMe)
  expect_equal(getFlowMatrix(readLIM), FMe)
})

# Test getVariables function
test_that("the getVariables function gives right answer", {
  expect_equal(getVariables(readLIM, web = lim_solved$X), variables)
})






#BM <- c(30, 20, 10, 5) ; names(BM) <- fwnames
#AE <- c(NA, NA, 0.3, 0.3) ; names(AE) <- fwnames
#GE <- c(NA, NA, 0.3, 0.3) ; names(GE) <- fwnames
#DM <- matrix(1, nrow = 4, ncol = 4)
#rownames(DM) <- fwnames
#colnames(DM) <- fwnames
#DM["MEIO", "LABILE"] <- 2/3
#DM["MACRO", "LABILE"] <- 1/2
#FDM <- FM * DM
#dead <- list(names = c("LABILE", "REFRAC"), def = c("Def", "Def"), frac = DM)
#model <- list(
#  type = "EF", FM = FM, BM = BM, AE = AE, GE = GE, dead = dead
#)
## Expected answer
#JM <- matrix(c(0,
 #              0,
#               (FM[3,1] - FM[1,3] + FM[3,4]*(1-AE[4])*(FDM[4,1]/(FDM[4,1]+FDM[4,2]))) / BM[3],
#               (FM[4,1] - FM[1,4]) / BM[4],#
#
#               0,
#               0,
#               (FM[3,2] - FM[2,3] + FM[3,4]*(1-AE[4])*(FDM[4,2]/(FDM[4,1]+FDM[4,2]))) / BM[3],
#               (FM[4,2] - FM[2,4]) / BM[4],
#
 #              AE[3] * GE[3] * FM[1,3] / BM[1],
#               AE[3] * GE[3] * FM[2,3] / BM[2],
#               0,
#               -FM[3,4] / BM[4],

 #              AE[4] * GE[4] * FM[1,4] / BM[1],
#               AE[4] * GE[4] * FM[2,4] / BM[2],
#               AE[4] * GE[4] * FM[3,4] / BM[3],
#               0
#), nrow = 4, ncol = 4)
#rownames(JM) <- fwnames ; colnames(JM) <- fwnames

#test_that("the function works with defecation and mortality into same compartment", {
#  expect_equal(getJacobian(model),
 #              JM)
#})


# use solved lim model in getJacobian
# write tests for all individual functions (getFlowMatrix, getAEmatrix, etc.)
# getFM has bug.


# test against expected outcome.
