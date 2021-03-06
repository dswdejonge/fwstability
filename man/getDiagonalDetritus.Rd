% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_diagonal.R
\name{getDiagonalDetritus}
\alias{getDiagonalDetritus}
\title{Calculate diagonal for dead (detritus) compartments}
\usage{
getDiagonalDetritus(FM, BM, AE, dead_names)
}
\arguments{
\item{FM}{A named square flowmatrix, source compartments as rows,
sink compartments as columns. Should NOT contain external compartments. (required)}

\item{BM}{A named numeric vector with biomasses of all compartments, must be in the same
order as the flow matrix after externals are excluded. (required)}

\item{AE}{A named numeric vector with assimilation efficiencies of all
compartments, must be in the same order as the flow matrix after externals
are excluded. Must be a fraction i.e. between 0 and 1 (required)}

\item{dead_names}{Character vector with all names of detritus and nutrient
compartments (everything that is not fauna). (required)}
}
\value{
This function returns a named numeric vector with the diagonal values for the
detritus compartments in the food web (per unit time). It is important to review the units of the
input data. If the FM is biomass per unit time then BM must be just biomass. If FM is
biomass per unit time per surface area then BM must be biomass per surface area.
}
\description{
This function finds diagonal values for detritus compartments to be used in the Jacobian
matrix. The calculation is based on the total amount of assimilated detritus in all
consumers of detritus and the biomass of the detritus compartments, and is an
implementation of the equation from Neutel & Thorne (2014).
}
\references{
Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles
and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
651–661. https://doi.org/10.1111/ele.12266
}
