# equation 13; diagonal for species
# is the negative mortality rate of species
#diagonalSpecies <- function(flow_solutions, BM, dead){
#  result <- -getMortalityRate(flow_solutions, BM, dead)
#  result <- result[-dead]
#  return(result)
#}

# equation 14; diagonal for detritus
# the total assimilated detritus is all consumers,
# divided by the biomass of detritus
#diagonalDetritus <- function(FM, BM, AE, dead){
#  detritus_ingestion <- colSums(t(FM[dead,-dead])*AE[-dead])
#  result <- -detritus_ingestion/BM[dead]
#}

#diagonal <- diag(FM)
#aii <- diagonalSpecies(flow_solutions = pars$X, BM, dead)
#add <- diagonalDetritus(FM[-externals, -externals], BM, AE, dead)
#diagonal[names(aii)] <- aii
#diagonal[names(add)] <- add
#diag(JM2) <- diagonal[-externals]
#JM2[is.na(JM2)] <- 0
