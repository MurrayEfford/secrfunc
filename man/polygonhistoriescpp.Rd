\name{polygonhistoriescpp}
\alias{polygonhistoriescpp}
\alias{polygonfxicpp}

\title{Probability of Polygon Capture Histories}

\description{

\code{polygonhistoriescpp} is the core function for the contribution of 
individual histories to likelihood. It is called by internal \pkg{secr} function 
\code{generalsecrloglik()}.

\code{polygonfxicpp} is the core code for 'fxi' calculations in \pkg{secr} when 
the detector type is polygon or transect based. It is called by the internal 
\pkg{secr} function \code{allhistpolygonfxi}.

}

\usage{

polygonhistoriescpp(nc, detectfn, grain, ncores, minp, binomN, w, xy, start, 
    group, hk, H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, debug)
    
polygonfxicpp(nc, detectfn, grain, ncores, minp, binomN, w, xy, start, 
    group, hk, H, gsbval, pID, mask, density, PIA, Tsk, h, hindex)

}

\arguments{

\item{nc}{integer number of individuals}
\item{detectfn}{integer detector code (see secr)}
\item{grain}{integer grain argument of RcppParallel::parallelFor}
\item{ncores}{integer number of threads}
\item{minp}{double minimum probability}

\item{binomN}{integer vector binomial code}
\item{w}{3-D integer array of capture histories (individual x occasion x detector)}
\item{xy}{2-column matrix of detection location coordinates}
\item{start}{integer vector of index of row in xy for each first detection}
\item{group}{integer vector; code for group of each individual (zero-based)}

\item{hk}{numeric array of detection hazards from makegkPolygoncpp}
\item{H}{numeric unbounded integrated hazard from radial function from makegkPolygoncpp}
\item{gsbval}{matrix of real parameter values, one row per combination}
\item{pID}{numeric matrix occasions x nmix. NOT USED?}
\item{mask}{2-column matrix of mask coordinates}

\item{density}{matrix relative density at each mask point (row) for each group? (column)}
\item{PIA}{integer vector for the slice of the parameter index array for a given mixture component and session}
\item{Tsk}{detector usage (numeric matrix detectors x occasions)}
\item{h}{for exclusive detectors, pre-computed hazard; otherwise -1}
\item{hindex}{for exclusive detectors, index to precomputed hazard; otherwise -1}
\item{debug}{logical; if TRUE...}

}

\details{

The capture history 'w' merely records detection/nondetection on a particular 
occasion. 'xy' records the location of each positive detection.

}

\value{

For \code{polygonhistoriescpp} --   
Numeric vector of probabilities, one per individual.

For \code{polygonfxicpp} --  
Numeric matrix with probability density for each individual (rows) 
at each mask point (columns).

}

\seealso{

\code{\link{makegkPolygoncpp}}

}

\examples{

datafilename <- system.file("extdata/testdata.RData", package = "secrfunc")
load(datafilename)  # loads w, dimw, xy, start, traps, mask, gsbval, cumk, binomN

detectfn <- 14   # hazard half-normal
dim <- 2
convex <- TRUE
grain <- 1
ncores <- 2

gkhk <- makegkPolygoncpp(detectfn, dim, convex, grain, ncores, gsbval, cumk, traps, mask)
# see tests for expected result

nc      <- dimw[1]  # individuals
S       <- dimw[2]  # occasions
K       <- dimw[3]  # single polygon
minp    <- 1e-200
group   <- rep(0,nc)
pID     <- matrix(1, nrow = S, ncol = 1)
density <- matrix(1/nrow(mask), nrow(mask), 1)
PIA     <- as.integer(array (1, dim = c(nc, S, K)))
Tsk     <- matrix(1, nrow = K, ncol = S)
h       <- matrix(-1)
hindex  <- matrix(-1)
debug   <- FALSE

prw <- polygonhistoriescpp(nc, detectfn, grain, ncores, minp, binomN, w, xy, start, 
    group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, debug)
# see tests for expected result

fxi <- polygonfxicpp(nc, detectfn, grain, ncores, minp, binomN, w, xy, start, 
    group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex)
# see tests for expected result

}