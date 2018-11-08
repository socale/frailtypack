
\name{plot.Diffepoce}
\alias{plot.Diffepoce}

\title{Plot difference of EPOCE estimators between two joint frailty models.}

\description{
  Plots values of the difference of two Cross-Validated Prognosis Observed Loss (CVPOL) computed with two joint frailty models. Confidence intervals are allowed.
}

\usage{
\method{plot}{Diffepoce}(x, conf.bands=TRUE, Xlab = "Time", Ylab = "EPOCE difference"
, ...)
}

\arguments{
  \item{x}{An object inheriting from \code{Diffepoce} class.}
  \item{conf.bands}{Logical value. Determines whether confidence intervals will be plotted. The default is FALSE.}
  \item{Xlab}{Label of x-axis. Default is '"Time"'}
  \item{Ylab}{Label of y-axis. Default is '"EPOCE difference"'}
  \item{\dots}{Other unused arguments.}
}

\value{
  Print one plot with one curve and its confidence interval.
}
\seealso{
 \code{\link{Diffepoce}}
}
\keyword{file}
