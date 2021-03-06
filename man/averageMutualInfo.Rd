% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_characteristics.R
\name{averageMutualInfo}
\alias{averageMutualInfo}
\title{Get average mutual information.}
\usage{
averageMutualInfo(FM)
}
\arguments{
\item{FM}{(required) A square flow matrix with flows from source in rows to sink in columns.}
}
\value{
Returns a double.
}
\description{
Calculates average mutual information.
}
\details{
Cannibalistic flows (from i to i) are included in calculations.
}
\references{
\itemize{
\item{
van Altena, C., Hemerik, L., de Ruiter, P.C., 2016. Food web stability and weighted
connectance: the complexity-stability debate revisited. Theor. Ecol. 9, 49–58.
https://doi.org/10.1007/s12080-015-0291-7
}
\item{
Boit, A.; Gaedke, U., 2014. Benchmarking successional progress in a quantitative food web.
Plos One 9(2):e90404. doi:10.1371/journal.pone.0090404
}
\item{Ulanowicz, R.E., 1997. Limitations on the connectivity of ecosystem flow networks. In:
Rinaldo, A.; Marani, A. (eds) Biological models. Instituto Veneto de Scienze, Lettre ed Arti,
Venica, pp 125-143.}
}
}
