% --- Source file: summaryjointPenal.Rd ---
\name{summary.jointPenal}
\Rdversion{1.1}
\alias{summary.jointPenal}
\alias{print.summary.jointPenal}


\title{summary of parameter estimates of a joint frailty model}

\description{
  This function returns hazard rations (HR) and its confidence intervals.
}
\usage{
\method{summary}{jointPenal}(object, level = 0.95, len = 6, d = 2, lab="hr", ...)
}

\arguments{
  \item{object}{output from a call to frailtyPenal for joint models} 
  \item{level}{significance level of confidence interval. Default is 95\%.}
  \item{d}{the desired number of digits after the decimal point. Default of 6 digits is used.}
  \item{len}{the total field width. Default is 6.}
  \item{lab}{label of printed results.}
  \item{\dots}{ other unused arguments.}
}


\value{
  Prints HR and its confidence intervals for each covariate. Confidence level is allowed (level argument).
}


\seealso{
 \code{\link{frailtyPenal}}
}

\examples{

\dontrun{

data(readmission)

#-- gap-time
modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+
charlson+terminal(death),formula.terminalEvent=~sex+dukes+charlson,
data=readmission,n.knots=14,kappa=c(9.55e+9,1.41e+12))

#-- calendar time
modJoint.calendar <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+
sex+dukes+charlson+terminal(death),formula.terminalEvent=~sex+dukes+charlson,
data=readmission,n.knots=10,kappa=c(9.55e+9,1.41e+12),recurrentAG=TRUE)

#-- It takes around 1 minute to converge

summary(modJoint.gap)
summary(modJoint.calendar)

}

}

\keyword{methods}
