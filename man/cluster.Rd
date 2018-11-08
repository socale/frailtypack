
\name{cluster}
\Rdversion{1.1}
\alias{cluster}
\title{
Identify clusters
}
\description{
This is a special function used in the context of the models for grouped data. It identifies correlated groups of observations defined by using 'cluster' function, and is used of 'frailtyPenal' formula for fitting univariate and joint models.
}

\usage{
cluster(x)
}

\arguments{
\item{x}{A character, factor, or numeric variable which is supposed to indicate 
the variable group}
}

\value{
\item{x}{A variable identified as a cluster }
}

\seealso{
 \code{\link{frailtyPenal}}
}

\examples{

\dontrun{

data(readmission)
modSha <- frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id),
n.knots=10,kappa=10000,data=readmission,hazard="Splines")

print(modSha)

}

}

\keyword{misc}
