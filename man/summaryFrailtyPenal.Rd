
\name{summary.frailtyPenal}
\Rdversion{1.1}

\alias{summary.frailtyPenal}
\alias{print.summary.frailtyPenal}


\title{summary of parameter estimates of a shared frailty model}

\description{
  This function returns hazard rations (HR) and its confidence intervals
}
\usage{
\method{summary}{frailtyPenal}(object, level = 0.95, len = 6, d = 2, lab="hr", ...)
}

\arguments{
  \item{object}{output from a call to frailtyPenal.} 
  \item{level}{significance level of confidence interval. Default is 95\%.}
  \item{d}{the desired number of digits after the decimal point. Default of 6 digits is used.}
  \item{len}{the total field width. Default is 6.}
  \item{lab}{label of printed results.}
  \item{\dots}{ other unused arguments.}
}


\value{
  Prints HR and its confidence intervals. Confidence level is allowed (level argument).
}


\seealso{
 \code{\link{frailtyPenal}}
}

\examples{

\dontrun{

data(kidney)

##-- Shared frailty model --##

modSha <- frailtyPenal(Surv(time,status)~age+sex+cluster(id),
n.knots=8,kappa=10000,data=kidney,hazard="Splines")

##-- Cox proportional hazard model --##

modCox <- frailtyPenal(Surv(time,status)~age+sex,
n.knots=8,kappa=10000,data=kidney,hazard="Splines")

#-- confidence interval at 95% level (default)

summary(modSha)
summary(modCox)

#-- confidence interval at 99% level

summary(modSha,level=0.99)
summary(modCox,level=0.99)

}

}

\keyword{methods}

