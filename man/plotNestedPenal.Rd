% --- Source file: plot.nestedPenal.Rd ---
\name{plot.nestedPenal}
\Rdversion{1.1}
\alias{plot.nestedPenal}
\alias{lines.nestedPenal}


\title{Plot Method for a Nested frailty model.}

\description{
  Plots estimated baseline survival and hazard functions (output from an object of class 'NestedPenal' for nested frailty models). Confidence bands are allowed. 
}

\usage{

\method{plot}{nestedPenal}(x, type.plot="Hazard", conf.bands=TRUE,
	pos.legend="topright", cex.legend=0.7, main, 
	color=2, Xlab = "Time", Ylab = "Hazard function", ...)

}

\arguments{
  \item{x}{ A nested model, i.e. an object of class \code{frailtyPenal} for Nested frailty models (output from calling \code{frailtyPenal}
    function).}
  \item{type.plot}{ a character string specifying the type of curve. Possible
                   value are "Hazard", or "Survival". The default
                   is "Hazard". Only the first words are required, e.g "Haz", "Su"}
  \item{conf.bands}{ logical value. Determines whether confidence bands will be plotted. The default is to do so.}
  \item{pos.legend}{The location of the legend can be specified by setting this argument to a single keyword from the list '"bottomright"', '"bottom"', '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
     '"center"'. The default is '"topright"'}
  \item{cex.legend}{character expansion factor *relative* to current
          'par("cex")'. Default is 0.7}
  \item{main}{plot title}
  \item{color}{curve color (integer)} 
  \item{Xlab}{Label of x-axis. Default is '"Time"'}  
  \item{Ylab}{Label of y-axis. Default is '"Hazard function"'}
  \item{\dots}{ Other graphical parameters like those in \code{\link{plot.frailtyPenal}}}
}

\value{
      Print a plot of the baseline survival or hazard functions with the confidence bands or not (conf.bands argument)
}

\seealso{
 \code{\link{frailtyPenal}}
}


\examples{

\dontrun{

data(dataNested)
modNested <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
subcluster(subgroup)+cov1+cov2,data=dataNested,n.knots=8,
kappa=50000,hazard="Splines")

plot(modNested,conf.bands=FALSE)

}

}

\keyword{methods}
