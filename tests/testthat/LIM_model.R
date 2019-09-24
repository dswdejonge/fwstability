### Import model ###
library(LIM)
readLIM <- Read(system.file("extdata", "foodweb2.lim", package = "fwstability"))
# Original FM
model <- list(
  type = "LIM",
  LIM = readLIM
)
# Use netto FM
model1 <- list(
  type = "LIM",
  LIM = readLIM,
  netto = TRUE
)

### Expected answers ###
lim <- Setup(readLIM)
lim_solved <- Ldei(lim)
fwnames <- c("DEADLABILE", "DEADREFRAC", "MEIO", "MACRO", "IN", "OUT")
nC = lim$NExternal + lim$NComponents
# Original FM
FM <- matrix(c(
  0,      0, 12.25, 0, 0, 0,
  0,      0,     0, 0, 0, 1.0375,
  4.9625, 1.0375,0, 5, 0, 1.25,
  4.4,    0,     0, 0, 0, 0.6,
  2.8875, 0,     0, 0, 0, 0,
  0,      0,     0, 0, 0, 0
), nrow = nC, ncol = nC, byrow = TRUE)
rownames(FM) <- fwnames ; colnames(FM) <- fwnames
# Use netto FM
FM1 <- FM - t(FM)
FM1[which(FM1 < 0)] <- 0

variables <- c(lim_solved$X["meioDefLab"] + lim_solved$X["meioDefRefrac"],
               lim_solved$X["macroDefLab"] + lim_solved$X["macroDefRefrac"],
               lim_solved$X["meioRespOne"] + lim_solved$X["meioRespTwo"],
               lim_solved$X["meioMoort"],
               lim_solved$X["macroOne"] + lim_solved$X["macroTwo"],
               lim_solved$X["macroPredMeio"] - lim_solved$X["macroDefLab"] - lim_solved$X["macroDefRefrac"],
               lim_solved$X["macroPredMeio"] - lim_solved$X["macroDefLab"] - lim_solved$X["macroDefRefrac"],
               lim_solved$X["meioGrazDet"] - lim_solved$X["meioDefLab"] - lim_solved$X["meioDefRefrac"],
               lim_solved$X["meioGrazDet"] - lim_solved$X["meioDefLab"] - lim_solved$X["meioDefRefrac"]
               - lim_solved$X["meioRespOne"] - lim_solved$X["meioRespTwo"],
               lim_solved$X["macroPredMeio"] - lim_solved$X["macroDefLab"] - lim_solved$X["macroDefRefrac"]
               - lim_solved$X["macroResp"]
)
names(variables) <- lim$Variables

AE <- c(NA, NA,
        variables["meioAss"]/sum(FM[,"MEIO"], na.rm = T),
        variables["macroAss"]/sum(FM[,"MACRO"], na.rm = T)
)
GE <- c(NA, NA,
        variables["meioGrowth"]/variables["meioAss"],
        variables["macroGrowth"]/variables["macroAss"]
)
names(AE) <- lim$Components$name ; names(GE) <- lim$Components$name

BM <- c(30, 20, 10, 5) ; names(BM) <- fwnames[1:4]
DM <- matrix(0, nrow = 4, ncol = 4)
rownames(DM) <- fwnames[1:4]
colnames(DM) <- fwnames[1:4]
DM["MEIO", "DEADLABILE"] <- lim_solved$X["meioDefLab"] / FM["MEIO", "DEADLABILE"]
DM["MACRO", "DEADLABILE"] <- lim_solved$X["macroDefLab"] / FM["MACRO", "DEADLABILE"]
FDM <- FM[1:4,1:4] * DM
# Original FM
JM <- matrix(c(0,
               0,
               (FM[3,1] - FM[1,3] + FM[3,4]*(1-AE[4])*(FDM[4,1]/(FDM[4,1]+FDM[4,2]))) / BM[3],
               (FM[4,1] - FM[1,4]) / BM[4],#

               0,
               0,
               (FM[3,2] - FM[2,3] + FM[3,4]*(1-AE[4])*(FDM[4,2]/(FDM[4,1]+FDM[4,2]))) / BM[3],
               (FM[4,2] - FM[2,4]) / BM[4],

               AE[3] * GE[3] * FM[1,3] / BM[1],
               AE[3] * GE[3] * FM[2,3] / BM[2],
               0,
               -FM[3,4] / BM[4],

               AE[4] * GE[4] * FM[1,4] / BM[1],
               AE[4] * GE[4] * FM[2,4] / BM[2],
               AE[4] * GE[4] * FM[3,4] / BM[3],
               0
), nrow = 4, ncol = 4)
rownames(JM) <- fwnames[1:4] ; colnames(JM) <- fwnames[1:4]
# Use netto FM
JM1 <- matrix(c(0,
               0,
               (FM1[3,1] - FM1[1,3] + FM[3,4]*(1-AE[4])*(FDM[4,1]/(FDM[4,1]+FDM[4,2]))) / BM[3],
               (FM1[4,1] - FM1[1,4]) / BM[4],#

               0,
               0,
               (FM1[3,2] - FM1[2,3] + FM[3,4]*(1-AE[4])*(FDM[4,2]/(FDM[4,1]+FDM[4,2]))) / BM[3],
               (FM1[4,2] - FM1[2,4]) / BM[4],

               AE[3] * GE[3] * FM1[1,3] / BM[1],
               AE[3] * GE[3] * FM1[2,3] / BM[2],
               0,
               -FM1[3,4] / BM[4],

               AE[4] * GE[4] * FM1[1,4] / BM[1],
               AE[4] * GE[4] * FM1[2,4] / BM[2],
               AE[4] * GE[4] * FM1[3,4] / BM[3],
               0
), nrow = 4, ncol = 4)
rownames(JM1) <- fwnames[1:4] ; colnames(JM1) <- fwnames[1:4]
