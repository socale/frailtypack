
\name{summary.additivePenal}
\Rdversion{1.1}
\alias{summary.additivePenal}
\alias{print.summary.additivePenal}

\title{summary of parameter estimates of an additive frailty model}

\description{
  This function returns hazard ratios (HR) and its confidence intervals
}
\usage{
\method{summary}{additivePenal}(object, level = 0.95, len = 6, d = 2, lab="hr", \dots)
}

\arguments{
  \item{object}{output from a call to additivePenal.} 
  \item{level}{significance level of confidence interval. Default is 95\%.}
  \item{d}{the desired number of digits after the decimal point. Default of 6 digits is used.}
  \item{len}{the total field width. Default is 6.}
  \item{lab}{label of printed results.}
  \item{\dots}{ other unused arguments.}
}


\value{
  Prints HR and its confidence intervals for each covariate. Confidence level is allowed (level argument) 
}


\seealso{
 \code{\link{additivePenal}}
}

\examples{

\dontrun{

data(dataAdditive)

modAdd <- additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
correlation=TRUE,data=dataAdditive,n.knots=8,kappa=862,hazard="Splines")

#- 'var1' is boolean as a treatment variable.

summary(modAdd)

}

}


\keyword{methods}
