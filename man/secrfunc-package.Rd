\name{secrfunc-package}
\alias{secrfunc-package}
\alias{secrfunc}

\docType{package}

\title{Helper Functions for Package \pkg{secr} }

\description{

This package provides some functions used internally by the spatial 
  capture-recapture package \pkg{secr} (from version 5.4.0). The idea is to speed up
  the installation of \pkg{secr}, and possibly reduce its size. Initially the 
  functions are those for area and transect search that use numerical integration 
  code from RcppNumerical (Qiu et al. 2023) and RcppEigen. 
  
  The functions are not intended to be user-friendly and require considerable 
  pre-processing of data. Examples and tests use data pre-formatted and saved 
  in the 'extdata' folder.
  
}

\details{

\tabular{ll}{
Package: \tab secrfunc\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2025-11-26\cr
License: \tab GNU General Public License Version 2 or later\cr
}

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. (2011) Estimation of population density by spatially
explicit capture--recapture with area searches. \emph{Ecology}
\bold{92}, 2202--2207.

Efford, M. G. (2025) The SECR book. A handbook of spatially explicit 
capture--recapture methods. Version 1.0.1. Zenodo \doi{10.5281/zenodo.15109938}.
Online at \url{https://murrayefford.github.io/SECRbook/}. 

Qiu, Y., Balan, S., Beall, M., Sauder, M., Okazaki, N., and Hahn, T. (2023). RcppNumerical: 'Rcpp'
  Integration for Numerical Computing Libraries. R package version 0.6-0.
  \url{https://CRAN.R-project.org/package=RcppNumerical}
}

\keyword{ package }

\seealso{

\code{\link{hdotpolycpp}} \cr 
\code{\link{makegkPolygoncpp}}  \cr
\code{\link{polygonhistoriescpp}}  \cr
\code{\link{polygonfxicpp}} 

}
