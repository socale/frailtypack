#data(gastadj)
#som=NULL
#a=frailtypack:::test(gastadj,nrow(gastadj),ncol(gastadj),som)

test = function(donnee,nbrow,nbcol,som){
  # ans <- .Fortran(C_test,
  #                 as.double(donnee[,1]),
  #                 as.double(donnee[,2]),
  #                 as.double(donnee[,3]),
  #                 as.double(donnee[,4]),
  #                 as.double(donnee[,5]),
  #                 as.double(donnee[,6]),
  #                 as.integer(nbrow),
  #                 som=rep(0,nbrow),
  #                 PACKAGE="frailtypack"
  #       )
  som=rep(0,nbrow)
  ans <- .Fortran(C_test,
                  as.matrix(donnee[,1:6]),
                  as.integer(nbrow),
                  som=som, # important d'initialiser les sortie dans l'appel de la fonction
                  PACKAGE="frailtypack"
  )
  return(ans$som)
}