# Check data format of matrix FM / JM
checkMformat <- function(M) {
  if(dim(M)[1] != dim(M)[2]) {
    stop("Input matrix is not square")
  } else if(is.null(colnames(M)) | is.null(rownames(M))) {
    stop("Input matrix must have named rows and columns.")
  } else if(!all(colnames(M) == rownames(M))) {
    stop("Input matrix must have same names in rows and columns.")
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
    if(is.null(rownames(m)) | is.null(colnames(m))){
      stop("All required matrices must be named.")
    }
    if(!all(n == rownames(m)) | !all(n == colnames(m))) {
      stop("The names and their order must be equal in all named vectors and matrices.")
    }
  }
  for(v in vectors) {
    if(is.null(names(v))) {
      stop("All required vectors must be named.")
    }
    if(!all(n == names(v))) {
      stop("The names and their order must be equal in all named vectors and matrices.")
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

checkDeadFormat <- function(dead) {
  if(!is.null(dead)) {
    if(!is.list(dead) | is.null(names(dead))) {
      stop("argument \"dead\" must be a named list")
    } else if(is.null(dead$names)) {
      stop("\"names\" element is required in the \"dead\" list")
    } else if(length(dead) > 3) {
      stop(paste("the list \"dead\" should have 3 elements at most"))
    } else if(!is.null(dead$def) &&
              length(dead$names) !=
              length(which(dead$def == "Def" | dead$def == "noDef"))) {
      stop("the second element of the list \"dead\" may only contain the strings \"Def\" and \"noDef\"")
    }
  }
}
