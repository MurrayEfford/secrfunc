\name{makegkPolygoncpp}
\alias{makegkPolygoncpp}
\title{makegkPolygoncpp (List return type)}
\usage{
makegkPolygoncpp(
  detectfn,
  dim,
  convex,
  grain,
  ncores,
  gsbval,
  cumk,
  traps,
  mask
)
}
\arguments{
\item{detectfn}{integer detector code (see secr)}

\item{dim}{integer}

\item{convex}{logical; if TRUE all detector polygons are convex}

\item{grain}{integer}
\item{ncores}{integer number of threads}
\item{gsbval}{real parameter values}
\item{cumk}{integer vector}
\item{traps}{2-column matrix; coordinates of detector polygon vertices}
\item{mask}{2-column matrix; coordinates of mask points}
}
\value{
An R list (based on C++ return type: List).
}
\description{
A longer description goes here.
A longer description goes here.
}
