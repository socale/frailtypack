##' @export
print.predJointNested <- function(x, digits = 3, ...) 
{
	if(class(x)!="predJointNested") stop("The object x must be a class predJoint.")
		
	if (!is.null(cl <- x$call)){
		cat("Call:\n")
		dput(cl)
		cat("\n")
	}	
	cat("Prediction of a terminal event given the history of the patient\n")
	cat("---------------------------------------------------------------------\n")
	cat("\n")
	print(x$pred,row.names=F,digits=digits)	
}
