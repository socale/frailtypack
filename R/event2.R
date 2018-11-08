#' Identify event2 indicator
#' 
#' This is a special function used in the context of multivariate frailty model
#' with two types of recurrent events and a terminal event (e.g., censoring
#' variable related to both recurrent events). It contains the indicator of the
#' recurrent event of type 2, normally 0=no event, 1=event, and is used on the
#' right hand side of a formula of a 'multivPenal' object.  Using
#' \code{event2()} in a formula implies that a multivariate frailty model for
#' two types of recurrent events and a terminal event is fitted.
#' 
#' 
#' @usage event2(x)
#' @param x A numeric variable but should be a boolean which equals 1 if the
#' subject has experienced an event of type 2 and 0 if not.
#' @return \item{x}{an indicator for an event of type 2}
#' @seealso \code{\link{multivPenal}}
#' @keywords misc
#' @export
event2<-function(x)
 {
  x
 }
