\name{makegkPolygoncpp}
\alias{makegkPolygoncpp}

\title{Pre-compute Detection Hazard}

\description{

Hazards and probabilities are computed for each parameter combination, detector and mask point.
Used by \pkg{secr} internal functions \code{generalsecrloglik()} and \code{secr_makegk()}.

}

\usage{
makegkPolygoncpp(detectfn, dim, convex, grain, ncores, gsbval, cumk, traps, mask)
}

\arguments{

\item{detectfn}{integer detector code (see \pkg{secr})}
\item{dim}{integer; 1 = transect, 2 = polygon}
\item{convex}{logical; if TRUE all detector polygons are convex}

\item{grain}{integer grain argument of RcppParallel::parallelFor}
\item{ncores}{integer number of threads}
\item{gsbval}{matrix of real parameter values, one row per combination}

\item{cumk}{integer vector of start indices of vertices for each detector}
\item{traps}{2-column matrix; coordinates of detector polygon vertices}
\item{mask}{2-column matrix; coordinates of mask points}

}

\value{

List with components
\item{H}{unbounded integrated hazard from radial function}
\item{gk}{probability vector for elements of cc x detector x mask array}
\item{hk}{hazard vector for elements of cc x detector x mask array}

}

\examples{

detectfn <- 14   # hazard half-normal
dim <- 2
convex <- TRUE
grain <- 1
ncores <- 2
# for a single parameter combination
gsbval <- matrix(c(0.2,1), nrow = 1)
cumk <- c(0,5)   
# centred unit polygon 
traps <- matrix(c(0,0,1,1,0,0,1,1,0,0), ncol = 2) - 0.5
mask <- as.matrix(expand.grid(x=-2:2, y=-2:2))

makegkPolygoncpp(detectfn, dim, convex, grain, ncores, gsbval, cumk, traps, mask)

}