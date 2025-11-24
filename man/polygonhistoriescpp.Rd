\name{polygonhistoriescpp}
\alias{polygonhistoriescpp}
\title{polygonhistoriescpp (NumericVector return type)}
\usage{
polygonhistoriescpp(
  nc,
  detectfn,
  grain,
  ncores,
  minp,
  binomN,
  w,
  xy,
  start,
  group,
  hk,
  H,
  gsbval,
  pID,
  mask,
  density,
  PIA,
  Tsk,
  h,
  hindex,
  debug
)
}
\arguments{

\item{nc}{integer number of individuals}
\item{detectfn}{integer detector code (see secr)}
\item{grain}{The input argument (type: const int grain).}
\item{ncores}{}
\item{minp}{double minimum probability}

\item{binomN}{integer vector binomial code}
\item{w}{3-D integer array of capture histories}
\item{xy}{The input argument (type: const NumericMatrix xy).}
\item{start}{}
\item{group}{integer vector; code for group of each individual}

\item{hk}{}
\item{H}{}
\item{gsbval}{}
\item{pID}{}
\item{mask}{}

\item{density}{}
\item{PIA}{}
\item{Tsk}{}
\item{h}{}
\item{hindex}{}

\item{debug}{}

}

\value{
An R list (based on C++ return type: NumericVector).
}

\description{
A longer description goes here.
}
