% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specresponses.R
\name{specresponses}
\alias{specresponses}
\title{Draw multiple species response curves}
\usage{
specresponses(matrix, var, main, xlab, method = "env", axis = 1, df = 5,
  bw = FALSE)
}
\arguments{
\item{matrix}{Community data, a matrix-like object with samples in rows and species in columns. Response curves are drawn for all selected columns (species).}

\item{var}{Vector containing environmental variable (per plot) \strong{OR} \code{vegan} ordination result object if \code{method = "ord"}.}

\item{main}{Optional: Main title.}

\item{xlab}{Optional: Label of x-axis.}

\item{method}{Method defining the type of variable. Default \code{method = "env"} fits a response curve to environmental variables. Alternatively \code{method = "ord"} fits a response along ordination axes.}

\item{axis}{Ordination axis (only if \code{method = "ord"}).}

\item{df}{Desired equivalent number of degrees of freedom (trace of the smoother matrix).}

\item{bw}{If set on \code{TRUE} the lines will be drawn in black/white with different line types instead of colours.}
}
\description{
This function draws multiple species response curves for rough interpretation of species responses to environmental variables or ordination axes.
It is based on \code{\link[stats]{smooth.spline}} which fits a cubic smoothing spline to the supplied data.
In contrast to \code{\link{specresponse}} this function can draw multiple curves into one plot, but will not draw points.
}
\section{Details}{

For response curves based on environmental variables the argument \code{var} takes a single vector containing the variable corresponding to the species abundances.

For a response to ordination axes (\code{method = "ord"}) the argument \code{var} requires a \code{vegan} ordination result object (e.g. from \code{\link[vegan]{decorana}}, \code{\link[vegan]{cca}}, \code{\link[vegan]{rda}} or \code{\link[vegan]{metaMDS}}).
First axis is used as default.

If you plot a response curve for only one species the use of \code{\link{specresponse}} is recommended for correct labels and the display of points.
}
\examples{
## Species (columns) need to be selected; call names() to get column numbers
names(meadows)
## Draw multiple species response curves on environmental variable in black/white
specresponses(meadows[ ,c(2,44,67,132)], env$pH, bw=T)

## Draw multiple species response curves on environmental variable with
## custom x-axis label and lower df
specresponses(meadows[ ,c(2,4,7)], env$pH, df=3, xlab = "Soil pH")

## Draw multiple species response curves on ordination axes
specresponses(meadows[ ,c(10:20)], meadows.dca, method="ord")
specresponses(meadows[ ,c(10:20)], meadows.dca, method="ord", axis = 2)

## Plot with manually log-transformed abundances
specresponses(log(meadows[ ,c(2,4,7)]+1), env$pH)

}
\author{
Friedemann Goral \email{fgoral@gwdg.de}
}
\seealso{
\code{\link{specresponse}}
}

