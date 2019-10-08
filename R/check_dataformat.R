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
