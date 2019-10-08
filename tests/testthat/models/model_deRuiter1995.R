# Proper data format, from De Ruiter et al. 1995
foodweb_data <- list(
  # Compartment names
  compartment = c(
    "Detritus",
    "Roots",
    "Fungi",
    "Bacteria",
    "Phytophagous_nematodes",
    "Collembola",
    "Cryptostigmatic_mites",
    "Noncryptostigmatic_mites",
    "Fungivorous_nematodes",
    "Enchytraeids",
    "Bacteriophagous_nematodes",
    "Flagellates",
    "Bacteriophagous_mites",
    "Amoebae",
    "Predatory_nematodes",
    "Nematophagous_mites",
    "Predatory_collembola",
    "Predatory_mites"
  ),
  # Biomasses(kg ha-1)
  BM = c(
    2500,
    300,
    2.13,
    227.5,
    0.19,
    0.47,
    0.01,
    0.02,
    0.08,
    0.43,
    0.30,
    0.53,
    0.001,
    11.53,
    0.06,
    0.004,
    0.03,
    0.0635
  ),
  MR = c(
    NA,
    1.00,
    1.20,
    1.20,
    1.08,
    1.84,
    1.20,
    1.84,
    1.92,
    5.00,
    2.68,
    6.00,
    1.84,
    6.00,
    3.00,
    1.84,
    1.84,
    1.84
  ),
  AE = c(
    NA,
    NA,
    1.00,
    1.00,
    0.25,
    0.50,
    0.50,
    0.50,
    0.38,
    0.25,
    0.60,
    0.95,
    0.50,
    0.95,
    0.50,
    0.90,
    0.50,
    0.60
  ),
  GE = c(
    NA,
    NA,
    0.30,
    0.30,
    0.37,
    0.35,
    0.35,
    0.35,
    0.37,
    0.40,
    0.37,
    0.40,
    0.35,
    0.40,
    0.37,
    0.35,
    0.35,
    0.35
  )
)
names(foodweb_data$BM) <- foodweb_data$compartment
names(foodweb_data$MR) <- foodweb_data$compartment
names(foodweb_data$AE) <- foodweb_data$compartment
names(foodweb_data$GE) <- foodweb_data$compartment

FM <- matrix(
  0,
  nrow = length(foodweb_data$compartment),
  ncol = length(foodweb_data$compartment),
  byrow = T
)
if(T){
  FM[1,3] <- 40.8
  FM[1,4] <- 1606
  FM[2,5] <- 7.44
  FM[3,6] <- 6.65
  FM[3,7] <- 0.08
  FM[3,8] <- 0.34
  FM[3,9] <- 2.51
  FM[3,10] <- 0.09
  FM[1,10] <- 11.1
  FM[4,10] <- 10.1
  FM[4,11] <- 6.89
  FM[4,12] <- 9.46
  FM[4,13] <- 0.01
  FM[4,14] <- 182
  FM[11,15] <- 0.42
  FM[4,15] <- 0.33
  FM[12,15] <- 0.008
  FM[12,14] <- 0.42
  FM[14,15] <- 0.17
  FM[5,18] <- 0.06
  FM[5,17] <- 0.13
  FM[5,16] <- 0.01
  FM[5,15] <- 0.28
  FM[6,18] <- 0.31
  FM[7,18] <- 0.005
  FM[8,18] <- 0.015
  FM[9,18] <- 0.03
  FM[9,17] <- 0.05
  FM[9,16] <- 0.004
  FM[9,15] <- 0.11
  FM[13,18] <- 0.003
  FM[11,18] <- 0.10
  FM[11,17] <- 0.20
  FM[11,16] <- 0.014
  FM[15,18] <- 0.02
  FM[15,17] <- 0.04
  FM[15,16] <- 0.003
  FM[16,18] <- 0.002
  FM[17,18] <- 0.02
}

rownames(FM) <- foodweb_data$compartment
colnames(FM) <- foodweb_data$compartment

# Model
# The authors state that intragroup interferences could not be established empirically.
# Therefore, the diagonal is set to zero.
# In the original paper the interaction strengths related to the detritus compartment were not
# calculated with a function specific for detritus.
# Therefore, the results from the paper can be duplicated by NOT defining
# detritus as a dead compartment.
model <- list(
  type = "EF",
  FM = FM,
  BM = foodweb_data$BM,
  AE = foodweb_data$AE,
  GE = foodweb_data$GE
)
