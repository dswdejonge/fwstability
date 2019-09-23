#' Depth-First-Search all loops
#'
#' This recursive Depth-First-Search function finds all loops in the foodweb starting from
#' a certain node and writes them to a text file.
#' @param AM Adjacency matrix.
#' @param node Current node.
#' @param visited Logical vector which stores if a node is visited before.
#' @param pathway Integer vector that tracks the pathway to current node.
#' @param started Logical vector which stores which nodes have been used
#' before as starting node.
#' @param output String with filename to store loops.
#' @param start Is integer with the current start node.
#' @references \itemize{
#' \item{
#' Reference?
#' }
#' }
#' @return No returned value, writes .txt file to working directory.
#' @export
dfsall <- function(AM, node, visited, pathway, started, output, start) {
  if(visited[node]){
    if(node == start){
      cycle <- c(pathway, start)
      write(cycle,
            file = output,
            append = TRUE,
            ncolumns = length(cycle))
    }
    return()
  }
  visited[node] <- TRUE
  pathway <- c(pathway, node)
  children <- which(AM[node,] == 1)
  for(child in children){
    dfsall(AM, child, visited, pathway, started, output, start)
  }
  visited[node] <- FALSE
  pathway[-length(pathway)]
}

#' Depth-First-Search loops of length k
#'
#' This recursive Depth-First-Search function finds all loops of length k starting at a
#' certain node in the foodweb and writes them to a text file.
#' @param AM Adjacency matrix.
#' @param node Current node.
#' @param visited Logical vector which stores if a node is visited before.
#' @param pathway Integer vector that tracks the pathway to current node.
#' @param k Integer denoting the length of loops searched for.
#' @param started Logical vector which stores which nodes have been used
#' before as starting node.
#' @param output String with filename to store loops.
#' @param start Is integer with the current start node.
#' @references \itemize{
#' \item{
#' Reference?
#' }
#' }
#' @return No returned value, writes .txt file to working directory.
#' @export
dfsk <- function(AM, node, visited, pathway, k, started, output, start){
  if(started[node]){
    return()
  }
  if(visited[node]){
    if(node == start && length(pathway) == k){
      cycle <- c(pathway, start)
      write(cycle,
            file = output,
            append = TRUE,
            ncolumns = length(cycle))
      return()
    }
    return()
  }
  visited[node] <- TRUE
  pathway <- c(pathway, node)
  children <- which(AM[node,] == 1)
  if(length(pathway) == k){
    if(!(start %in% children)){
      return()
    }else{
      dfsk(AM, start, visited, pathway, k, started, output, start)
    }
  }else{
    for(child in children){
      dfsk(AM, child, visited, pathway, k, started, output, start)
    }
  }
  visited[node] <- FALSE
  pathway[-length(pathway)]
}

#' Depth-First-Search wrapper
#'
#' This function is a wrapper for two recursive Depth-First-Search functions.
#' The wrapper checks data input, redirects to the right algorithm, and checks if no
#' file will be overwritten.
#' @param AM (required) Adjacency matrix
#' @param k (optional) Integer of length loop to search for.
#' @param output (optional) Filename to store loops. Default is "allLoops.txt"
#' or "allLoops_k=<k>.txt".
#' @param verbose (optional) Default is TRUE. Set to FALSE to hide messages.
#' @references \itemize{
#' \item{
#' Reference?
#' }
#' }
#' @return Returns file name.
#' @export
dfs <- function(AM, k = NULL, output = NULL, verbose = T){
  if(dim(AM)[1] != dim(AM)[2]){
    stop("Adjacency matrix must be square.")
  } else if(!all(AM == 1 | AM == 0)) {
    stop("Adjancency matrix can only contain 0 and 1.")
  } else if(!is.null(k) & !is.double(k)) {
    stop("k must be an integer.")
  } else if(!is.null(output) & !is.character(output)) {
    stop("output should contain string with filename for output.")
  }

  visited <- logical(length = dim(AM)[1])
  pathway <- numeric()
  started <- logical(length = dim(AM)[1])

  if(is.null(k) & is.null(output)) {
    output <- "allLoops.txt"
  } else if(!is.null(k) & is.null(output)) {
    output <- paste0("allLoops_k=",k,".txt")
  }

  if(file.exists(output)) {
    stop(paste0("file ",output," already exists. Rename or remove to avoid overwriting."))
  } else if(verbose) {
    message(paste0("Loops stored as ",output))
  }

  if(is.null(k)) {
    for(startnode in 1:dim(AM)[1]){
      dfsall(AM, node = startnode, visited, pathway, started, output, start = startnode)
      started[startnode] <- TRUE
    }
  } else {
    for(startnode in 1:dim(AM)[1]){
      dfsk(AM, node = startnode, visited, pathway, k, started, output, start = startnode)
      started[startnode] <- TRUE
    }
  }
  return(output)
}
