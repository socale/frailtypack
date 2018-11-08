\name{summary.multivPenal}
\Rdversion{1.1}
\alias{summary.multivPenal}

\title{
summary of parameter estimates of a multivariate frailty model.
}
\description{
  This function returns hazard ratio (HR) and its confidence intervals.
}
\usage{
\method{summary}{multivPenal}(object, level = 0.95, len = 6, d = 2, lab = "hr", ...)
}

\arguments{
  \item{object}{output from a call to multivPenal for joint multivariate models} 
  \item{level}{significance level of confidence interval. Default is 95\%.}
  \item{d}{the desired number of digits after the decimal point. Default of 6 digits is used.}
  \item{len}{the total field width. Default is 6.}
  \item{lab}{label of printed results.}
  \item{\dots}{other unused arguments.}
}

\value{
  Prints HR and its confidence intervals for each covariate. Confidence level is allowed (level argument) 
}

\seealso{
 \code{\link{multivPenal}}
}

\keyword{methods}
\keyword{multiv}
