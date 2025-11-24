\name{hdotpolycpp}
\alias{hdotpolycpp}

\title{Net Detection Hazard for Polygon and Transect Detectors)
}

\description{

Compute spatially explicit net probability of detection for individual(s) at 
given coordinates. Used by secr function \code{pdot}.

}

\usage{

hdotpolycpp(xy, traps, Tsk, markocc, cumk, detectfn, gsb, convex, dim, grain, ncores)

}

\arguments{

\item{xy}{2-column matrix of coordinates}

\item{traps}{2-column matrix; coordinates of detector polygon vertices}

\item{Tsk}{detector usage (numeric matrix detectors x occasions)}

\item{markocc}{integer vector of marking occasion codes (1 = marking)}

\item{cumk}{integer vector of start indices of vertices for each detector}

\item{detectfn}{integer detector code (see secr)}

\item{gsb}{numeric vector of real parameter values (lambda0, sigma)}

\item{convex}{logical; if TRUE all detector polygons are convex (non-convex untested)}

\item{dim}{integer; 1 = transect, 2 = polygon}

\item{grain}{integer grain argument of RcppParallel::parallelFor}

\item{ncores}{integer number of threads}

}

\details{

Uses 'integrate' function from RcppNumerical.
Sighting occasions are ignored.

}

\value{

Rcpp NumericVector with one hazard for each row in 'xy'.

}

\examples{

set.seed(123)

# 50 random points
xy <- matrix(runif(100), ncol = 2) * 2 - 1

# centred unit polygon 
traps <- matrix(c(0,0,1,1,0,0,1,1,0,0), ncol = 2) - 0.5

Tsk <- matrix(1, nrow = 1, ncol = 5)
markocc <- rep(1,5)
cumk <- c(0,5,0)   # zero-terminated
detectfn <- 14   # hazard half-normal
gsb <- c(lambda0 = 0.2, sigma = 0.5)

h <- hdotpolycpp(xy, traps, Tsk, markocc, cumk, detectfn, gsb, 
    convex = TRUE, dim = 2, grain = 1, ncores = 2)

}