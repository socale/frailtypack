
# Evaluated computer time occording to the number of threads used to compute an arithmetic sequence. 
# We used a Fortran code  combined with OpenMP
# 
# The main objective of this function was to compare computer burden when increase the number of thread. 
# In the example, we observe inthe associated plot that when the maximum term of the sequence is very large, 
# we do not decrease anymore the computer time with a number of threads set greater than the number of core 
# 
# @aliases evalOpenMPFortran
# @usage
# evalOpenMPFortran(nbrthread = 5, maxSum=100000000)
# 
# @param nbrthread The mximum number ot threads to considered
# @param maxSum The the maximum term of the suite
#
# @return a dataframe containing the computer times associated to each number of thread
# @keywords Fortran OpenMP R
# @author 
# Casimir Ledoux SOFEU <scl.ledoux@gmail.com>
# @export
#
# @seealso \code{\link{testOpenMPFortran}}
# @examples 
# 
# \dontrun{
# 
# ###- with default parameters
# evalOpenMPFortran()
# 
# ###- with given parameters
# evalOpenMPFortran(5,100000000)
# 
# }
# 
evalOpenMPFortran <- function(nbrthread = 5, maxSum=100000000){
  nbr=nbrthread
  result=data.frame(matrix(0,nrow = nbr,ncol=2))
  names(result)=c("Number of threads","computer time(s)")
  for(i in 1:nbr){
    result[i,1]=i
    T1<-Sys.time() 
    a=testOpenMPFortran(2,4,maxSum,i,F)
    T2<-Sys.time()
    result[i,2]=as.double(T2-T1)
  }
  plot(result,type="b")
  return(result)
}
  