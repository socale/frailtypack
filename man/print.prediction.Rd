
\name{print.prediction}
\alias{print.predFrailty}
\alias{print.predJoint}
\alias{print.predLongi}

\title{
Print a short summary of results of prediction function.
}
\description{
Print a short summary of results of prediction function.
}
\usage{
\method{print}{predFrailty}(x, digits = 3, ...)
\method{print}{predJoint}(x, digits = 3, ...)
\method{print}{predLongi}(x, digits = 3, ...)
}

\arguments{
  \item{x}{An object from the 'prediction' function, objects inheriting from \code{predFrailty}, \code{predJoint} and \code{predLongi} classes.}
  \item{digits}{Number of digits to print}
  \item{\dots}{Other unused arguments}
}

\value{

Print the probabilities estimated.

}

\seealso{
 \code{\link{prediction}}
}

\keyword{methods}

