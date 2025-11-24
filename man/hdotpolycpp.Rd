\name{hdotpolycpp}
\alias{hdotpolycpp}
\title{hdotpolycpp (NumericVector return type)}
\usage{
hdotpolycpp(
  xy,
  traps,
  Tsk,
  markocc,
  cumk,
  detectfn,
  gsb,
  convex,
  dim,
  grain,
  ncores
)
}
\arguments{
\item{xy}{2-column matrix of coordinates}

\item{traps}{coordinates of detector polygon vertices}

\item{Tsk}{detector usage}

\item{markocc}{integer vector of marking occasion codes}

\item{cumk}{integer vector}

\item{detectfn}{integer detector code (see secr)}

\item{gsb}{real parameter values}

\item{convex}{logical; if TRUE all detector polygons are convex}

\item{dim}{integer vector}

\item{grain}{integer}

\item{ncores}{ integer number of threads}
}

\value{
An R list (based on C++ return type: NumericVector).
}
\description{
A longer description goes here.
}
