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
checkNamingFormat <- function(matrices, vectors) {
  for(m in matrices) {
    if(is.null(rownames(m)) | is.null(colnames(m))){
      stop("All required matrices must be named.")
    }
  }
  for(v in vectors) {
    if(is.null(names(v))) {
      stop("All required vectors must be named.")
    }
  }
}
