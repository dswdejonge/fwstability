% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_characteristics.R
\name{maxNrLoops}
\alias{maxNrLoops}
\title{Maximum number of loops}
\usage{
maxNrLoops(n, k = NULL)
}
\arguments{
\item{n}{(required) Integer. Number of network compartments.}

\item{k}{(optional) Integer. Length of loop. Default NULL searches of loops of all length
(k=2 to k = n).}
}
\value{
Returns the maximum number of possible loops.
}
\description{
This function calculates the maximum number of loops, assuming full connectance.
}
\details{
A fully connected network is assumed. The total number of loops of length k
can be found as n! / (n - k)!. If k is not given the total number of possible loops (i.e.
of length k = 2 to k = n) is found.
}
