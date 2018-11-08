% --- Source file: summarynestedPenal.Rd ---
\name{summary.nestedPenal}
\Rdversion{1.1}
\alias{summary.nestedPenal}
\alias{print.summary.nestedPenal}


\title{summary of regression coefficient estimates of a nested frailty model}

\description{
  This function returns hazard rations (HR) and its confidence intervals for each regression coefficient.
}
\usage{
\method{summary}{nestedPenal}(object, level = 0.95, len = 6, d = 2, lab="hr", ...)
}

\arguments{
  \item{object}{output from a call to nestedPenal.} 
  \item{level}{significance level of confidence interval. Default is 95\%.}
  \item{d}{the desired number of digits after the decimal point. Default of 6 digits is used.}
  \item{len}{the total field width. Default is 6.}
  \item{lab}{label of printed results.}
  \item{\dots}{ other unused arguments.}
}


\value{
  Prints HR and its confidence intervals for each regression coefficient. Confidence level is allowed (level argument).
}


\seealso{
 \code{\link{frailtyPenal}}
}

\examples{

\dontrun{

data(dataNested)

modNested <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
subcluster(subgroup)+cov1+cov2,data=dataNested,
n.knots=8,kappa=c(50000,50000),hazard="Splines")

#- It takes 90 minutes to converge (depends on processor)

summary(modNested)

}

}

\keyword{methods}
