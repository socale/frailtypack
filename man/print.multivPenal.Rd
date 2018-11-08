\name{print.multivPenal}
\Rdversion{1.1}
\alias{print.multivPenal}

\title{
Print a Short Summary of parameter estimates of a multivariate frailty model
}
\description{
  Prints a short summary of parameter estimates of a multivariate frailty model, or more generally an object of class 'multivPenal'.
}
\usage{
\method{print}{multivPenal}(x, digits = max(options()$digits - 4, 6), ...)
}

\arguments{
  \item{x}{the result of a call to the multivPenal function}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unused arguments }
}

\value{

Print, separately for each type of event (recurrent1, recurrent2 and terminal), the parameter estimates of the survival or hazard functions.

}

\seealso{
 \code{\link{multivPenal}}
}


\keyword{methods}
\keyword{multiv}
