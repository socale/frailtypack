\name{print.frailtyPenal}

\Rdversion{1.1}
\alias{print.frailtyPenal}

\title{Print a Short Summary of parameter estimates of a shared frailty model}

\description{
  Prints a short summary of parameter estimates of a 'frailtyPenal' object 
}

\usage{
\method{print}{frailtyPenal}(x, digits = max(options()$digits - 4, 6), ...)
}

\arguments{
  \item{x}{the result of a call to the frailtyPenal function.}
  \item{digits}{number of digits to print.}
  \item{\dots}{other unused arguments.}
}

\value{

Print the parameter estimates of the survival or hazard functions.

}

\seealso{
 \code{\link{frailtyPenal}}
}

\keyword{methods}
