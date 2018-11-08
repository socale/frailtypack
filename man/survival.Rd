\name{survival}
\alias{survival}

\title{
Survival function
}
\description{
Let t be a continuous variable, we determine the value of the survival function to t after run fit.
}
\usage{
survival(t, ObjFrailty)
}

\arguments{
  \item{t}{
time for survival function.
}
  \item{ObjFrailty}{
an object from the frailtypack fit.
}
}

\value{
return the value of survival function in t.
}

\examples{

\dontrun{

#-- a fit Shared
data(readmission)

fit.shared <- frailtyPenal(Surv(time,event)~dukes+cluster(id)+
strata(sex),n.knots=10,kappa=c(10000,10000),data=readmission)

#-- calling survival
survival(20,fit.shared)

}

}
