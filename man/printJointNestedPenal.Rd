\name{print.jointNestedPenal}
\Rdversion{1.1}

\alias{print.jointNestedPenal}

\title{Print a Short Summary of parameter estimates of a joint nested frailty model}
\description{
  Prints a short summary of parameter estimates of a joint nested frailty model, or more generally an object of class 'jointNestedPenal' for joint nested frailty models. 
}

\usage{
\method{print}{jointNestedPenal}(x, digits = max(options()$digits - 4, 6), ...)
}

\arguments{
  \item{x}{the result of a call to the jointNestedPenal function}
  \item{digits}{number of digits to print}
  \item{\dots}{other unused arguments}
}

\value{

Print, separately for each type of event (recurrent and terminal), the parameter estimates of the survival or hazard functions.

}

\seealso{
 \code{\link{frailtyPenal}}
}

\keyword{methods}
