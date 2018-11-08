
\name{print.longiPenal}
\Rdversion{1.1}

\alias{print.longiPenal}

\title{Print a Summary of parameter estimates of a joint model for longitudinal data and a terminal event}
\description{
  Prints a short summary of parameter estimates of a joint model for longitudinal data and a terminal event, an object inheriting from class 'longiPenal'.
}

\usage{

\method{print}{longiPenal}(x, digits = max(options()$digits - 4, 6), ...)

}

\arguments{
  \item{x}{an object inheriting from \code{longiPenal} class}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unused arguments }
}

\value{

Print, separately for each part of the model (longitudinal and terminal) the parameter estimates and details on the estimation.

}


\seealso{
 \code{\link{longiPenal}}
}

\keyword{methods}
