
\name{print.trivPenalNL}
\Rdversion{1.1}

\alias{print.trivPenalNL}

\title{Print a Summary of parameter estimates of a non-linear trivariate joint model for longitudinal data, recurrent events and a terminal event}
\description{
  Prints a short summary of parameter estimates of a non-linear trivariate joint model for longitudinal data, recurrent events and a terminal event,  an object inheriting from class 'trivPenalNL'.
}

\usage{

\method{print}{trivPenalNL}(x, digits = max(options()$digits - 4, 6), ...)

}

\arguments{
  \item{x}{an object inheriting from \code{trivPenaNLl} class}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unused arguments }
}

\value{

Print, separately for each part of the model (biomarker growth, biomarker decline, recurrent events and terminal event) the parameter estimates and details on the estimation.

}


\seealso{
 \code{\link{trivPenalNL}}
}

\keyword{methods}
