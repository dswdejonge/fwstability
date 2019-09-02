#################
# findVariables #
#################
# Find value of variables based on the flow solution values.
getsVariables <- function(LIM = NULL,
                          readLIM = NULL,
                          flowvalues = NULL) {
  # Create vector to store data in
  vals <-  numeric(LIM$NVariables)
  # Get all needed data
  # the vector with the variable values (varvec) is refreshed after every loop
  vareq <- readLIM$vars       # contains matrix which defines the variable equations
  parvec <- readLIM$pars$val  # vector with all parameter values in the right order
  flowvec <- flowvalues       # vector with all parsimonious flow values (also contains levels)

  for (i in 1:LIM$NVariables) {
    # refresh vector with variable parsimonious values
    varvec <- vals # vector with all variables values in right order: needs to be refreshed every time

    # get subset with the same equation nr
    subset <- vareq[vareq$nr == i,]

    # takes parameter, variable or flow nr from subset to use as index
    # to find the corresponding parsimonious values and
    # multiplies them by 'val' which is the coefficient (like 1 or -1) and sums them
    sum <-
      sum(parvec[subset$par1]*subset$val, na.rm = TRUE) +
      sum(parvec[subset$par2]*subset$val, na.rm = TRUE) +
      sum(parvec[subset$par3]*subset$val, na.rm = TRUE) +
      sum(parvec[subset$par4]*subset$val, na.rm = TRUE) +
      sum(varvec[subset$var]*subset$val, na.rm = TRUE) +
      sum(flowvec[subset$flow]*subset$val, na.rm = TRUE)

    # replaces 0 in data frame to actual sum value
    vals[i] = sum
  }
  vals <- data.frame(LIM$Variables, vals)
  return(vals)
}
