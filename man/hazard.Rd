\name{hazard}
\alias{hazard}

\title{
Hazard function.
}
\description{
Let t be a continuous variable, we determine the value of the hazard function to t after run fit.
}
\usage{
hazard(t, ObjFrailty)
}

\arguments{
  \item{t}{
time for hazard function.
}
  \item{ObjFrailty}{
an object from the frailtypack fit.
}
}

\value{
return the value of hazard function in t.
}

\examples{

\dontrun{

#-- a fit Shared
data(readmission)
fit.shared <- frailtyPenal(Surv(time,event)~dukes+cluster(id)+
strata(sex),n.knots=10,kappa=c(10000,10000),data=readmission)

#-- calling survival
hazard(20,fit.shared)

}

}