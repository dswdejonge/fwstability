# Check data format of matrix FM / JM
checkMformat <- function(M) {
  if(dim(M)[1] != dim(M)[2]) {
    stop("Input matrix is not square")
  } else if(is.null(colnames(M)) | is.null(rownames(M))) {
    stop("Input matrix must have named rows and columns.")
  } else if(!all(colnames(M) == rownames(M))) {
    stop("Input matrix must have same names in rows and columns.")
  } else if(!is.numeric(M)) {
    stop("Input matrix must be numeric")
  } else if(length(which(M %in% NaN | M %in% NA)) > 0) {
    stop("NAs and NaNs not allowed in matrix")
  }
}

checkExternalsFormat <- function(externals, M) {
  if((FALSE %in% (externals %in% rownames(M))) |
     (FALSE %in% (externals %in% colnames(M)))) {
    stop("the names of the external compartments are unknown")
  }
}

# input is list
checkNamingFormat <- function(matrices = NULL, vectors = NULL) {
  n <- names(vectors[[1]])
  for(m in matrices) {
    if(is.null(m)) {next}
    if(is.null(rownames(m)) | is.null(colnames(m))){
      stop("All required matrices must be named.")
    }
    if(!all(n %in% rownames(m) & rownames(m) %in% n)) {
      stop("The names must be equal in all named vectors and matrices.")
    }
  }
  for(v in vectors) {
    if(is.null(v)) {next}
    if(is.null(names(v))) {
      stop("All required vectors must be named.")
    }
    if(!all(n %in% names(v) & names(v) %in% n )) {
      stop("The names must be equal in all named vectors and matrices.")
    }
  }
}

checkBMformat <- function(BM) {
  if((TRUE %in% is.na(BM)) | (TRUE %in% (BM <= 0)) | (!is.numeric(BM))) {
    stop("biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  }
}

checkDiagonalFormat <- function(diagonal, correct_length) {
  if(!is.numeric(diagonal) & all(diagonal != "model")) {
    stop("given diagonal not numeric or set to \"model\"")
  } else if(length(diagonal) != 1 & length(diagonal) != correct_length) {
    stop("given diagonal has incorrect length")
  }
}

# input is list
checkCEformat <- function(CE) {
  for(e in CE) {
    if(any(e > 1 | e < 0, na.rm = TRUE)) {
      stop("assimilation and growth efficiencies must lie between 0 and 1")
    }
  }
}

checkDeadFormat <- function(dead, FM) {
  if(!is.null(dead)) {
    if(!is.list(dead) | is.null(names(dead))) {
      stop("argument \"dead\" must be a named list")
    } else if(is.null(dead$names)) {
      stop("the \"names\" element is required in the \"dead\" list")
    } else if(is.null(dead$frac)) {
      stop("the \"frac\" element is required in the \"dead\" list")
    } else if(length(dead) > 2) {
      stop(paste("the list \"dead\" should have 2 elements at most"))
    } else if(FALSE %in% (dead$names %in% colnames(FM))) {
      stop("the names of the dead compartments are unknown")
    } else if(FALSE %in% (dim(FM) == dim(dead$frac))){
      stop("the FM and dead$frac matrix must have the same dimensions (check for external compartments in your dead$frac matrix)")
    } else {
      checkMformat(dead$frac)
    }
  }
}

checkMortalityFormat <- function(MR, dead) {
  if(!is.null(MR)) {
    if(TRUE %in% is.na(MR) && (is.null(dead) | !all(names(which(is.na(MR))) %in% dead))) {
      stop("Mortality rates of non-dead compartments cannot be NA.")
    } else if(min(MR, na.rm = T) <= 0 | !is.numeric(MR)) {
      stop("the MR vector contains values equal or smaller than zero, or is non-numeric")
    }
  }
}

checkStabilityMethod <- function(method, JM, MR) {
  if(method != "eigenvalue" & method != "scalar" & method != "initial") {
    stop("unknown method chosen")
  } else if(method == "scalar" & is.null(MR)) {
    stop("MR vector required for the scalar method")
  }
}
