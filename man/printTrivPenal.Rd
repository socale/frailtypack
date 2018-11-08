
\name{print.trivPenal}
\Rdversion{1.1}

\alias{print.trivPenal}

\title{Print a Summary of parameter estimates of a joint model for longitudinal data, recurrent events and a terminal event}
\description{
  Prints a short summary of parameter estimates of a joint model for longitudinal data, recurrent events and a terminal event,  an object inheriting from class 'trivPenal'.
}

\usage{

\method{print}{trivPenal}(x, digits = max(options()$digits - 4, 6), ...)

}

\arguments{
  \item{x}{an object inheriting from \code{trivPenal} class}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unused arguments }
}

\value{

Print, separately for each part of the model (longitudinal, recurrent and terminal) the parameter estimates and details on the estimation.

}


\seealso{
 \code{\link{trivPenal}}
}

\keyword{methods}
