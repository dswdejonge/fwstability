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
#' @param file String with name for textfile that will store loops.
#' @param start Is integer with the current start node.
#' @return No returned value, writes .txt file to working directory.
#' @export
dfsall <- function(AM, node, visited, pathway, started, file, start) {
  if(visited[node]){
    if(node == start){
      cycle <- c(pathway, start)
      write(cycle,
            file = file,
            append = TRUE,
            ncolumns = length(cycle))
    }
    return()
  }
  visited[node] <- TRUE
  pathway <- c(pathway, node)
  children <- which(AM[node,] == 1)
  for(child in children){
    dfsall(AM, node = child, visited, pathway, started, file, start)
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
#' @param file String with filename to store loops.
#' @param start Is integer with the current start node.
#' @return No returned value, writes .txt file to working directory.
#' @export
dfsk <- function(AM, node, visited, pathway, k, started, file, start){
  if(started[node]){
    return()
  }
  if(visited[node]){
    if(node == start && length(pathway) == k){
      cycle <- c(pathway, start)
      write(cycle,
            file = file,
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
      dfsk(AM, node = start, visited, pathway, k, started, file, start)
    }
  }else{
    for(child in children){
      dfsk(AM, node = child, visited, pathway, k, started, file, start)
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
#' @param output (required) Name for textfile to store loops.
#' The data will be stored as "<output>.txt" if all loops are searched for, or as
#' "<output>_k=<k>.txt". Default is "allLoops", so the text file will be "allLoops.txt"
#' or "allLoops_k=<k>.txt".
#' @param k (optional) Default is NULL. Integer of length loop to search for.
#' @param verbose (optional) Default is TRUE. Set to FALSE to hide messages.
#' @return Returns file name.
#' @references
#' Tarjan, R. (1972). Depth-First Search and Linear Graph Algorithms.
#' SIAM Journal on Computing, 1(2), 146–160.
#' https://doi.org/10.1137/0201010.
#' @export
dfs <- function(AM, output = "allLoops", k = NULL, verbose = T){
  # Check input data
  if(verbose){message("Checking input data.")}
  if(dim(AM)[1] != dim(AM)[2]){
    stop("Adjacency matrix must be square.")
  } else if(!all(AM == 1 | AM == 0)) {
    stop("Adjancency matrix can only contain 0 and 1.")
  } else if(!is.null(k) & !is.double(k)) {
    stop("k must be an integer.")
  } else if(is.null(output) | !is.character(output)) {
    stop("output should contain string with filename for output.")
  }

  # Initialize variables for dfs algorithm
  visited <- logical(length = dim(AM)[1])
  pathway <- numeric()
  started <- logical(length = dim(AM)[1])

  # Determine file name
  if(is.null(k)) {
    file <- paste0(output,".txt")
  } else {
    file <- paste0(output,"_k=",k,".txt")
  }
  if(file.exists(file)) {
    stop(paste0("The file ",file," already exists. Rename or remove to avoid overwriting."))
  } else if(verbose) {
    message(paste0("Loop output will be stored in ",file))
  }

  if(is.null(k)) {
    if(verbose){message("Start dfs to find all loops.")}
    for(startnode in 1:dim(AM)[1]){
      dfsall(AM, node = startnode,
             visited = visited, pathway = pathway, started = started,
             file = file, start = startnode)
      started[startnode] <- TRUE
    }
  } else {
    if(verbose){message(paste0("Start dfs to find loops of length k=",k))}
    for(startnode in 1:dim(AM)[1]){
      dfsk(AM, node = startnode,
           visited = visited, pathway = pathway, k = k, started = started,
           file = file, start = startnode)
      started[startnode] <- TRUE
    }
  }
  return(file)
}
