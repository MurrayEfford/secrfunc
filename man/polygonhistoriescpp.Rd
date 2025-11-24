\name{polygonhistoriescpp}
\alias{polygonhistoriescpp}

\title{Probability of Polygon Capture Histories}

\description{

Core function for the contribution of individual histories to likelihood.
Called by internal secr function 'generalsecrloglik'.

}

\usage{

polygonhistoriescpp(nc, detectfn, grain, ncores, minp, binomN, w, xy, start, 
    group, hk, H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, debug)

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
\item{group}{integer vector; code for group of each individual}

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

Numeric vector of probabilities, one per individual.

}

\seealso{

\code{\link{makegkPolygoncpp}}

}

\examples{

detectfn <- 14   # hazard half-normal
grain <- 1
ncores <- 2
gsbval <- matrix(c(0.2,1), nrow = 1)
cumk <- c(0,5)   
traps <- matrix(c(0,0,1,1,0,0,1,1,0,0), ncol = 2) - 0.5
mask <- as.matrix(expand.grid(x=-2:2, y=-2:2))

gkhk <- makegkPolygoncpp(detectfn, dim, TRUE, grain, ncores, gsbval, cumk, traps, mask)


nc <- 0
minp <- 1e-200
binomN <- 1
group <- rep(1,nc)

#w
#xy
#start

#polygonhistoriescpp(nc, detectfn, grain, ncores, minp, binomN, w, xy, start, 
#    group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, debug)

}