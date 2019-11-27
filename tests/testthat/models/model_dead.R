# Dead compartment, no externals
# Plants and animals both take up detritus and deposit into detritus.
# Animals also eat plants.
#
#       /---------4----------v
# DETRITUS -7-> PLANT -4-> ANIMAL
#  ^______3______/ <____8___/
#
## Proper data format
fwnames <- c("DETRITUS", "PLANT", "ANIMAL")
FM <- matrix(c(0, 0, 0, 10, 0, 0, 4, 4, 0), nrow = 3, ncol = 3)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
BM <- c(30, 20, 10) ; names(BM) <- fwnames
AE <- c(NA, 0.6, 0.3) ; names(AE) <- fwnames
GE <- c(NA, 0.8, 0.3) ; names(GE) <- fwnames
# Detritus feedback
plant_defecation <- sum(FM[,"PLANT"])*(1-AE["PLANT"])
animal_defecation <- sum(FM[,"ANIMAL"])*(1-AE["ANIMAL"])
plant_mortality <- sum(FM[,"PLANT"])*AE["PLANT"]*GE["PLANT"] - sum(FM["PLANT",])
animal_mortality <- sum(FM[,"ANIMAL"])*AE["ANIMAL"]*GE["ANIMAL"] - sum(FM["ANIMAL",])
FM["PLANT","DETRITUS"] <- plant_defecation + plant_mortality
FM["ANIMAL", "DETRITUS"] <- animal_defecation + animal_mortality
# Get fraction matrix: fraction of flow that comprises defecation
frac <- matrix(
  c(0, plant_defecation/FM[2,1], animal_defecation/FM[3,1], rep(0, 6)),
  ncol = 3, nrow = 3)
rownames(frac) <- fwnames ; colnames(frac) <- fwnames
dead <- list(names = "DETRITUS", frac = frac)
# Combine to model
model <- list(
  type = "EF", FM = FM, BM = BM, AE = AE, GE = GE, dead = dead
)
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
rownames(JM) <- fwnames ; colnames(JM) <- fwnames
####

### Inclusion of diagonal ###
DIAG <- c(-1,-2,-3)
model2 <- model
model2$diagonal <- DIAG
JM2 <- JM; JM2[c(1,5,9)] <- DIAG
###


### Inclusion of diagonal calculated from model ###
MR <- c(NA, 10/BM[2], 5/BM[3]) ; names(MR) <- fwnames
model3 <- model
model3$diagonal <- "model"
model3$MR <- MR
JM3 <- JM
JM3[c(1,5,9)] <- c(
  -1/BM[1] * (AE[2] * FM[1,2] + AE[3] * FM[1,3]),
  -MR[2],
  -MR[3]
)
###


### Netto FM
FM4 <- getNettoFM(FM, deadnames = "DETRITUS")
JM4 <- matrix(c(0,
               (FM4[2,1] - FM4[1,2] + FM[2,3]*(1-AE[3])) / BM[2],
               (FM[3,1] - FM[1,3]) / BM[3],

               AE[2] * GE[2] * FM[1,2] / BM[1],
               0,
               -FM4[2,3] / BM[3],

               AE[3] * GE[3] * FM[1,3] / BM[1],
               AE[3] * GE[3] * FM4[2,3] / BM[2],
               0
), nrow = 3, ncol = 3)
rownames(JM4) <- fwnames ; colnames(JM4) <- fwnames
model4 <- model
model4$FM <- FM4
