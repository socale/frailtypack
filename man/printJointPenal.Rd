
\name{print.jointPenal}
\Rdversion{1.1}

\alias{print.jointPenal}

\title{Print a Short Summary of parameter estimates of a joint frailty model}
\description{
  Prints a short summary of parameter estimates of a joint frailty model, or more generally an object of class 'frailtyPenal' for joint frailty models. 
}

\usage{

\method{print}{jointPenal}(x, digits = max(options()$digits - 4, 6), ...)

}

\arguments{
  \item{x}{the result of a call to the jointPenal function}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unused arguments }
}

\value{

Print, separately for each type of event (recurrent and terminal), the parameter estimates of the survival or hazard functions.

}


\seealso{
 \code{\link{frailtyPenal}}
}

\keyword{methods}
