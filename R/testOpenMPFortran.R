
# Used to teste pactice use of OpenMP with Fortran code in an R package.
# 
# In this function we test a call of a fortran subroutine in a R package. we therefore compute a sum of two numbers 
# using a fortran subroutine. In addition, we test the use of OpenMP by compute an arithmetic Sequences with first term 1 
# and the commun difference 1. 
#
# @aliases testOpenMPFortran
# @usage testOpenMPFortran(a,b,nboot,nbThread,stat=T)
# @param a The first term of the sum
# @param b The second term of the sum
# @param nboot The maximum term of the sequence 
# @param nbThread The number of thread to used in the OpenMP procedure
# @param stat A logical specifying if you want to print the result of the sum on screen or not
#
# @keywords Fortran OpenMP R
# @author Casimir Ledoux SOFEU <scl.ledoux@gmail.com>
# 
# @seealso \code{\link{evalOpenMPFortran}}
# @examples 
# 
# \dontrun{
# frailtypack:::testOpenMPFortran(2,4,10000,10,stat=T)
# 
# }
testOpenMPFortran = function(a,b,nboot,nbThread,stat=T){
  if(stat) print(paste("voila ce que vous avez passe en parametre: a=", a, "b=", b,sep=""))
  vOut=c(0,0,0)
  ans <- .Fortran(C_somme,
                  as.double(c(a,b)),
                  s=as.double(0),
                  as.integer(nboot),
                  as.integer(nbThread),
                  vOut=as.double(vOut),
                  PACKAGE="frailtypack"
                  )
  if(stat) {
    print(paste("la somme de vos nombre vaut",ans$s))
    #print(paste("voila vOut",c(ans$vOut[1],ans$vOut[2],ans$vOut[3])))
  }
  return(ans)
}