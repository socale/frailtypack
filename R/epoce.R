# nobs - number of recurrent event
# nsujet - number of individuals
# nsujety - number of repeated measurements




#' Estimators of the Expected Prognostic Observed Cross-Entropy (EPOCE) for
#' evaluating predictive accuracy of joint models.
#' 
#' 
#' This function computes estimators of the Expected Prognostic Observed
#' Cross-Entropy (EPOCE) for evaluating the predictive accuracy of joint models
#' using \code{frailtyPenal}, \code{longiPenal}, \code{trivPenal} or
#' \code{trivPenalNL}. On the same data as used for estimation of the joint
#' model, this function computes both the Mean Prognosis Observed Loss (MPOL)
#' and the Cross-Validated Prognosis Observed Loss (CVPOL), two estimators of
#' EPOCE. The latter corrects the MPOL estimate for over-optimism by
#' approximated cross-validation. On external, this function only computes
#' MPOL.
#' 
#' 
#' @usage epoce(fit, pred.times, newdata = NULL, newdata.Longi = NULL)
#' @param fit A jointPenal, longiPenal, trivPenal or trivPenalNL object.
#' @param pred.times Time or vector of times to compute epoce.
#' @param newdata Optional. In case of joint models obtained with
#' \code{frailtyPenal}, \code{trivPenal} or \code{trivPenalNL}. For models
#' inheriting from \code{trivPenal} or \code{trivPenalNL} class, if
#' \code{newdata} is given, \code{newdata.Longi} must be given as well.  When
#' missing, the data used for estimating the fit are used, and CVPOL and MPOL
#' are computed (internal validation). When \code{newdata} is specified, only
#' MPOL is computed on this new dataset (external validation). The new dataset
#' and the dataset used in the estimation must have the same covariates with
#' the same coding without missing data.
#' @param newdata.Longi Optional. In case of joint models obtained with
#' \code{longiPenal}, \code{trivPenal} or \code{trivPenalNL}. For models
#' inheriting from \code{longiPenal}, if the \code{newdata.Longi} is given,
#' \code{newdata} must be \code{NULL}, but for models from \code{trivPenal} or
#' \code{trivPenalNL} class, if \code{newdata.Longi} is given, \code{newdata}
#' must be provided as well. The two datasets newdata and newdata.Longi must
#' include the information concerning the same patients with the same
#' characteristics and the appropriate data on follow up (recurrences for
#' \code{newdata} and longitudinal measurements for \code{newdata.Longi}).
#' @return \item{data}{name of the data used to compute epoce}
#' \item{new.data}{a boolean which is FALSE if computation is done on the same
#' data as for estimation, and TRUE otherwise} \item{pred.times}{time or vector
#' of times used in the function} \item{mpol}{values of MPOL for each
#' pred.times} \item{cvpol}{values of CVPOL for each pred.times}
#' \item{IndivContrib}{all the contributions to the log-likelihood for each
#' pred.times} \item{AtRisk}{number of subject still at risk for each
#' pred.times}
#' @references D. Commenges, B. Liquet, C. Proust-Lima (2012). Choice of
#' prognostic estimators in joint models by estimating differences of expected
#' conditional Kullback-Leibler risks. \emph{Biometrics}, \bold{68(2)},
#' 380-387.
#' @keywords misc
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' ########################################
#' #### EPOCE on a joint frailty model ####
#' ########################################
#' 
#' data(readmission)
#' 
#' modJoint.gap <- frailtyPenal(Surv(t.start,t.stop,event)~ cluster(id) +
#'   dukes + charlson + sex + chemo + terminal(death),
#'   formula.terminalEvent = ~ dukes + charlson + sex + chemo ,
#'   data = readmission, n.knots = 8, kappa =c(2.11e+08,9.53e+11),
#'   recurrentAG=TRUE)
#' 
#' # computation on the same dataset
#' temps <- c(200,500,800,1100)
#' epoce <- epoce(modJoint.gap,temps)
#' 
#' print(epoce)
#' plot(epoce,type = "cvpol")
#' 
#' # computation on a new dataset
#' # here a sample of readmission with the first 50 subjects
#' s <- readmission[1:100,]
#' epoce <- epoce(modJoint.gap,temps,newdata=s)
#' 
#' print(epoce)
#' plot(epoce,type = "cvpol")
#' 
#' #################################################
#' #### EPOCE on a joint  model for a biomarker ####
#' #########   and a terminal event  ###############
#' #################################################
#' 
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' # Survival data preparation - only terminal events 
#' colorectalSurv <- subset(colorectal, new.lesions == 0)
#' 
#' modLongi <- longiPenal(Surv(time0, time1, state) ~ age +
#' treatment + who.PS, tumor.size ~  year*treatment + age +
#' who.PS, colorectalSurv, data.Longi =colorectalLongi,
#' random = c("1", "year"),  id = "id", link = "Random-effects", 
#' left.censoring = -3.33, hazard = "Weibull", 
#' method.GH = "Pseudo-adaptive")
#' 
#' # computation on the same dataset
#' time <- c(1, 1.5, 2, 2.5)
#' epoce <- epoce(modLongi,time)
#' 
#' print(epoce)
#' plot(epoce, type = "cvpol")
#' 
#' # computation on a new dataset
#' # here a sample of colorectal data with the first 50 subjects
#' s <-  subset(colorectal, new.lesions == 0 & id%in%1:50)
#' s.Longi <- subset(colorectalLongi, id%in%1:50)
#' epoce <- epoce(modLongi, time, newdata = s, newdata.Longi = s.Longi)
#' 
#' print(epoce)
#' plot(epoce, type = "cvpol")
#' 
#' 
#' ###################################################
#' #### EPOCE on a joint model for a biomarker, ######
#' #### recurrent events and a terminal event   ######
#' ###################################################
#' 
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' # Linear model for the biomarker
#' # (computation takes around 30 minutes)
#' model.trivPenalNL <-trivPenal(Surv(gap.time, new.lesions) ~ cluster(id)
#' + age + treatment + who.PS + prev.resection + terminal(state),
#' formula.terminalEvent =~ age + treatment + who.PS + prev.resection, 
#' tumor.size ~ year * treatment + age + who.PS, data = colorectal,
#' data.Longi = colorectalLongi, random = c("1", "year"), id = "id", 
#' link = "Random-effects", left.censoring = -3.33, recurrentAG = FALSE,
#' hazard = "Weibull", method.GH="Pseudo-adaptive", n.nodes=7)
#' 
#' # computation on the same dataset
#' time <- c(1, 1.5, 2, 2.5)
#' 
#' # (computation takes around 10 minutes)
#' epoce <- epoce(model.trivPenalNL,time)
#' print(epoce)
#' plot(epoce, type = "cvpol")
#' 
#' # computation on a new dataset
#' # here a sample of colorectal data with the first 100 subjects
#' s <-  subset(colorectal,  id%in%1:100)
#' s.Longi <- subset(colorectalLongi, id%in%1:100)
#' # (computation takes around 10 minutes)
#' epoce <- epoce(model.trivPenalNL, time, newdata = s, newdata.Longi = s.Longi)
#' 
#' print(epoce)
#' plot(epoce, type = "cvpol")
#' 
#' 
#' 
#' # Non-linear model for the biomarker
#' 
#' # No information on dose - creation of a dummy variable 
#' colorectalLongi$dose <- 1
#' 
#' # (computation can take around 40 minutes)
#' model.trivPenalNL <- trivPenalNL(Surv(time0, time1, new.lesions) ~ cluster(id) + age + treatment
#'  + terminal(state), formula.terminalEvent =~ age + treatment, biomarker = "tumor.size",
#'  formula.KG ~ 1, formula.KD ~ treatment, dose = "dose", time.biomarker = "year", 
#'  data = colorectal, data.Longi =colorectalLongi, random = c("y0", "KG"), id = "id", 
#'  init.B = c(-0.22, -0.16, -0.35, -0.19, 0.04, -0.41, 0.23), init.Alpha = 1.86,
#'  init.Eta = c(0.5, 0.57, 0.5, 2.34), init.Biomarker = c(1.24, 0.81, 1.07, -1.53),
#'  recurrentAG = TRUE, n.knots = 5, kappa = c(0.01, 2), method.GH = "Pseudo-adaptive")
#' 
#' # computation on the same dataset
#' time <- c(1, 1.5, 2, 2.5)
#' 
#' epoce <- epoce(model.trivPenalNL, time)
#' 
#' 
#' }
#' 
#' 
epoce <- function(fit, pred.times, newdata = NULL, newdata.Longi = NULL){

        if (missing(fit)) stop("The argument fit must be specified")
        if (class(fit)!="jointPenal" & class(fit)!="longiPenal" & class(fit)!="trivPenal" & class(fit)!="trivPenalNL") stop("The argument fit must be a class 'jointPenal', 'longiPenal' or 'trivPenal'")
        if (missing(pred.times)) stop("The argument pred.times must be specified")
        if (class(pred.times)!="numeric") stop("pred.times must contain numerical values")

        if(!missing(newdata) & (class(newdata)!="data.frame")) stop("The argument newdata must be a 'data.frame'")
        if(!missing(newdata.Longi) & (class(newdata.Longi)!="data.frame")) stop("The argument newdata must be a 'data.frame'")

  if(class(fit)== "jointPenal" & !missing(newdata.Longi))warning("The argument newdata.Longi is not required and thus ignored")
        if(class(fit)== "longiPenal" & !missing(newdata.Longi) & missing(newdata))warning("For an object of class 'longiPenal' both datasets should be given")
        if(class(fit)== "trivPenal" & !missing(newdata.Longi) & missing(newdata))warning("For an object of class 'trivPenal' both datasets should be given")
        if(class(fit)== "trivPenalNL" & missing(newdata.Longi) & !missing(newdata))warning("For an object of class 'trivPenalNL' both datasets should be given")


        nt <- length(pred.times)
        vopt <- fit$varHtotal
       
        b <- fit$b
        np <- length(fit$b)
        typeof <- fit$typeof
        nva <- fit$nvar

        if (typeof == 0){
                nz <- fit$n.knots
                zi <- fit$zi
        }else{
                nz <- 0
                zi <- 0
        }

        if (typeof == 1){
                nbintervR <- fit$nbintervR
                nbintervDC <- fit$nbintervDC
                ttt <- fit$time
                tttdc <- fit$timedc
        }else{
                nbintervR <- 0
                nbintervDC <- 0
                ttt <- 0
                tttdc <- 0
        }

        # recuperation des profils d'individus
        m <- fit$call
        m0 <- match.call()

        if (!missing(newdata)){
                if (length(colnames(eval(m$data)))!=length(colnames(eval(m0$newdata)))) stop("Your new dataset must have the same number of columns than the dataset used in the 'fit'")
                if (any(colnames(eval(m$data))!=colnames(eval(m0$newdata)))) stop("Your new dataset must have the very same variables than the dataset used in the 'fit'")
        }

        if (!missing(newdata.Longi)){
          if (length(colnames(eval(m$data.Longi)))!=length(colnames(eval(m0$newdata.Longi)))) stop("Your new dataset for longitudinal data must have the same number of columns than the dataset used in the 'fit'")
          if (any(colnames(eval(m$data.Longi))!=colnames(eval(m0$newdata.Longi)))) stop("Your new dataset for longitudinal data must have the very same variables than the dataset used in the 'fit'")
        }

        if (is.null(m$recurrentAG)) recurrentAG <- FALSE
        else recurrentAG <- TRUE

  if(class(fit)=="jointPenal" | class(fit)=="trivPenal" | class(fit) == "trivPenalNL"){
        #m$formula.terminalEvent <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa <- m$maxit <- m$hazard <- m$nb.int1 <-m$nb.int2 <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$init.Theta <- m$init.Alpha <- m$Alpha <- m$... <- NULL
        m$formula.LongitudinalData <- m$formula.terminalEvent <- m$recurrentAG <- m$data.Longi <- m$n.knots <- m$random <- m$link  <- m$id <- m$kappa <- m$maxit <- m$hazard <- m$nb.int  <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$left.censoring <- m$print.times <- m$init.Random <- m$init.Eta <- m$init.Alpha <- m$method.GH <- m$intercept <- m$n.nodes <- m$biomarker <- m$formula.KG <- m$formula.KD <- m$dose <- m$time.biomarker <- m$BoxCox  <-  m$init.Biomarker <- m$... <- m$RandDist <- NULL

        m[[1]] <- as.name("model.frame")
        if (!missing(newdata)) m[[3]] <- as.name(m0$newdata) # nouveau dataset

        dataset <- eval(m, sys.parent())

        typeofY <- attr(model.extract(dataset, "response"),"type")
        Y <- model.extract(dataset, "response")

        if (typeofY=="right"){
                tt0 <- rep(0,dim(dataset)[1])
                tt1 <- Y[,1]
                c <- Y[,2]
        }else{
                tt0 <- Y[,1]
                tt1 <- Y[,2]
                c <- Y[,3]
        }
        tt0 <- as.numeric(tt0)
        tt1 <- as.numeric(tt1)
        c <- as.numeric(c)

        class(m$formula) <- "formula"
        special <- c("strata", "cluster", "subcluster", "terminal", "num.id", "timedep")

        Terms <- terms(m$formula, special)#, data = m$data)

        m$formula <- Terms

        dropx <- NULL

        tempc <- untangle.specials(Terms, "cluster", 1:10)
        cluster <- strata(dataset[, tempc$vars], shortlabel = TRUE)
#       numbers <- table(cluster)[order(unique(cluster))]
#       newCluster <- rep(1:nsujet,numbers)
        dropx <- c(dropx,tempc$terms)

        tempterm <- untangle.specials(Terms, "terminal", 1:10)
        terminal <- strata(dataset[, tempterm$vars], shortlabel = TRUE)
        terminal <- as.numeric(as.character(terminal))
        dropx <- c(dropx,tempterm$terms)

        if (!is.null(dropx)) newTerms <- Terms[-dropx]
        else newTerms <- Terms
      
        X <- model.matrix(newTerms, dataset)
        if (ncol(X) > 1) X <- X[, -1, drop = FALSE]
        nva1 <- ncol(X)

        if (!missing(newdata)){
                nobs <- nrow(newdata)
                nsujet <- length(unique(cluster))
        }else{
                nobs <- fit$n
                nsujet <- fit$groups
        }

        if (!recurrentAG){
                tt0dc <- aggregate(tt1,by=list(cluster),FUN=sum)[,2]#rep(0,nsujet)
                tt1dc <- aggregate(tt1,by=list(cluster),FUN=sum)[,2]
        }else{
                tt0dc <- rep(0,nsujet)
                tt1dc <- aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
        }
        cdc <- aggregate(terminal,by=list(cluster),FUN=function(x) x[length(x)])[,2]

       
        m2 <- fit$call

        m2$formula.LongitudinalData <- m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int1 <-m2$nb.int2 <- m2$RandDist <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$init.Theta <- m2$init.Alpha <- m2$Alpha <- m2$method.GH <- m2$intercept <- m2$init.Eta <- m2$data.Longi <- m2$init.Random <- m2$left.censoring <- m2$random <- m2$link <- m2$id <- m2$n.nodes <- m2$biomarker <- m2$formula.KG <- m2$formula.KD <- m2$dose <- m2$time.biomarker <- m2$BoxCox  <-  m2$init.Biomarker <- m2$... <- NULL

        m2$formula[[3]] <- m2$formula.terminalEvent[[2]]
        m2$formula.terminalEvent <- NULL
        m2[[1]] <- as.name("model.frame")

        if (!missing(newdata)) m2[[3]] <- as.name(m0$newdata) # nouveau dataset

        datasetdc <- eval(m2, sys.parent())

        class(m2$formula) <- "formula"
        special2 <- c("strata", "timedep")
        Terms2 <- terms(m2$formula, special2)#, data = m3$data)

        X2 <- model.matrix(Terms2, datasetdc)
        if (ncol(X2) > 1) X2 <- X2[, -1, drop = FALSE]
        nva2 <- ncol(X2)

        if (!is.null(ncol(X2))){
                Xdc <- aggregate(X2[,1],by=list(cluster), FUN=function(x) x[length(x)])[,2]
                if (ncol(X2)>1){
                        for (i in 2:ncol(X2)){
                                Xdc.i <- aggregate(X2[,i],by=list(cluster), FUN=function(x) x[length(x)])[,2]
                                Xdc <- cbind(Xdc,Xdc.i)
                        }
                }
        }else{
                Xdc <- aggregate(X2,by=list(cluster), FUN=function(x) x[length(x)])[,2]
        }

        if (!missing(newdata) & length(fit$coef[1:(fit$nvarEnd+fit$nvarRec)])!=(ncol(X)+ncol(X2))) stop("Different covariates in model and newdata. Verify your dataset, be careful to the factor variables.")
  }
    if(class(fit) == "trivPenal" | class(fit) == "longiPenal"){

      m2 <- fit$call
      m2$formula <- m2$formula.terminalEvent <- m2$data <- m2$random <- m2$id <- m2$link <- m2$n.knots <- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$left.censoring <- m2$init.Random <- m2$init.Eta <- m2$method.GH <- m2$intercept <- m2$... <- NULL
      if (!missing(newdata.Longi)){m2[[3]] <- as.name(m0$newdata.Longi) # nouveau dataset
                                  data.Longi <- newdata.Longi
      }else{data.Longi <- eval(m2$data.Longi) }

      special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep")

      #========= Longitudinal Data preparation =========================
      class(m2$formula.LongitudinalData) <- "formula"

      TermsY <- terms(m2$formula.LongitudinalData, special, data = data.Longi)

      llY <- attr(TermsY, "term.labels")#liste des variables explicatives
      ord <- attr(TermsY, "order")

      #=========================================================>

      name.Y <- as.character(attr(TermsY, "variables")[[2]])
      yy <- data.Longi[,which(names(data.Longi)==name.Y)]

      # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
      ind.placeY <- which(llY%in%names(which(lapply(data.Longi[,which(names(data.Longi)%in%llY)],function(x) length(levels(x)))>2)))
      vec.factorY <- NULL
      vec.factorY <- c(vec.factorY,llY[ind.placeY])


      mat.factorY <- matrix(vec.factorY,ncol=1,nrow=length(vec.factorY))

      # Fonction servant a prendre les termes entre "as.factor"
      vec.factorY <-apply(mat.factorY,MARGIN=1,FUN=function(x){
        if (length(grep("as.factor",x))>0){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          return(substr(x,start=pos1,stop=pos2))
        }else{
          return(x)
        }})

      ind.placeY <- grep(paste(vec.factorY,collapse="|"),llY)

      if(is.factor(data.Longi[,names(data.Longi)==llY[1]]))X_L<- as.numeric(data.Longi[,names(data.Longi)==llY[1]])-1
      else X_L<- data.Longi[,names(data.Longi)==llY[1]]

      if(length(llY)>1){
        for(i in 2:length(llY)){
          if(is.factor(data.Longi[,names(data.Longi)==llY[i]]))X_L<- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llY[i]])-1)
          else X_L<- cbind(X_L,data.Longi[,names(data.Longi)==llY[i]])
        }}
      #X_L<- data.Longi[,names(data.Longi)%in%(llY)]

      if(sum(ord)>length(ord)){
        for(i in 1:length(ord)){
          if(ord[i]>1){
            v1 <- strsplit(as.character(llY[i]),":")[[1]][1]
            v2 <- strsplit(as.character(llY[i]),":")[[1]][2]

            if(is.factor(data.Longi[,names(data.Longi)==v1]) && length(levels(data.Longi[,names(data.Longi)==v1]))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            if(is.factor(data.Longi[,names(data.Longi)==v2]) && length(levels(data.Longi[,names(data.Longi)==v2]))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            if(is.factor(data.Longi[,names(data.Longi)==v1]) || !is.factor(data.Longi[,names(data.Longi)==v2])){
              X_L <- cbind(X_L,(as.numeric(data.Longi[,names(data.Longi)==v1])-1)*data.Longi[,names(data.Longi)==v2])
              llY[i]<-paste(llY[i],levels(data.Longi[,names(data.Longi)==v1])[2],sep="")
            }else if(!is.factor(data.Longi[,names(data.Longi)==v1]) || is.factor(data.Longi[,names(data.Longi)==v2])){
              X_L <- cbind(X_L,data.Longi[,names(data.Longi)==v1]*(as.numeric(data.Longi[,names(data.Longi)==v2])-1))
              llY[i]<-paste(llY[i],levels(data.Longi[,names(data.Longi)==v2])[2],sep="")
            }else{
              X_L <- cbind(X_L,data.Longi[,names(data.Longi)==v1]*data.Longi[,names(data.Longi)==v2])
            }

          }
        }
      }
     if(dim(X_L)[2]!=length(llY))stop("The variables in the longitudinal part must be in the data.Longi")
      X_L <- as.data.frame(X_L)
      names(X_L) <- llY

      Intercept <- rep(1,dim(X_L)[1])

      if(fit$intercept)X_L <- cbind(Intercept,X_L)

      X_Lall<- X_L
      "%+%"<- function(x,y) paste(x,y,sep="")

      if(length(vec.factorY) > 0){
        for(i in 1:length(vec.factorY)){
          X_L <- cbind(X_L[,-(which(names(X_L)==vec.factorY[i]))],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorY[i], collapse= "+")), data.Longi)[,-1])
        }


        vect.factY<-names(X_L)[which(!(names(X_L)%in%llY))]

        occurY <- rep(0,length(vec.factorY))

        for(i in 1:length(vec.factorY)){
          #occur[i] <- sum(vec.factor[i] == vect.fact)
          occurY[i] <- length(grep(vec.factorY[i],vect.factY))
        }
      }

      if (ncol(X_L) == 0){
        noVarY <- 1
      }else{
        noVarY <- 0
      }
      #=========================================================>
  
      clusterY <- data.Longi$id
      maxy_rep <- max(table(clusterY))
   
      uni.cluster<-as.factor(unique(clusterY))
      npred <- length(uni.cluster)


      nva3<-ncol(X_L)


      varY <- as.matrix(sapply(X_L, as.numeric))


      #=======================================>
      #======= Construction du vecteur des indicatrice
      if(length(vec.factorY) > 0){
        #               ind.place <- ind.place -1
        k <- 0
        for(i in 1:length(vec.factorY)){
          ind.placeY[i] <- ind.placeY[i]+k
          k <- k + occurY[i]-1
        }
      }

      if(fit$link=="Random-effects")link <- 1
      if(fit$link=="Current-level") link <- 2

      matzy <- NULL
      names.matzy <- fit$names.re

      matzy <- data.matrix(X_Lall[,which(names(X_Lall)%in%names.matzy)])


      if(fit$leftCensoring==FALSE){s_cag_id = 0
                                   s_cag = 0}else{
                                     s_cag_id = 1
                                     s_cag = fit$leftCensoring.threshold
                                   }

    }
  if(class(fit)== "longiPenal"){
    #m$formula.terminalEvent <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa <- m$maxit <- m$hazard <- m$nb.int1 <-m$nb.int2 <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$init.Theta <- m$init.Alpha <- m$Alpha <- m$... <- NULL
    m$formula.LongitudinalData <- m$formula.terminalEvent <- m$recurrentAG <- m$data.Longi <- m$n.knots <- m$random <- m$link  <- m$id <- m$kappa <- m$maxit <- m$hazard <- m$nb.int  <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$left.censoring <- m$print.times <- m$init.Random <- m$init.Eta <- m$init.Alpha <- m$method.GH <- m$intercept <- m$n.nodes <- m$... <- NULL

    m[[1]] <- as.name("model.frame")
    if (!missing(newdata)) m[[3]] <- as.name(m0$newdata) # nouveau dataset

    dataset <- eval(m, sys.parent())

    typeofY <- attr(model.extract(dataset, "response"),"type")
    Y <- model.extract(dataset, "response")


    if (typeofY=="right"){
      tt0dc <- rep(0,dim(dataset)[1])
      tt1dc <- Y[,1]
      cdc <- Y[,2]
    } else {
      tt0dc <- Y[,1]
      tt1dc <- Y[,2]
      cdc <- Y[,3]
    }

    tt0dc <- as.numeric(tt0dc)
    tt1dc <- as.numeric(tt1dc)
    cdc <- as.numeric(cdc)

    class(m$formula) <- "formula"
    special <- c("strata", "cluster", "subcluster", "num.id", "timedep")

    Terms <- terms(m$formula, special)#, data = m$data)

    m$formula <- Terms

    dropx <- NULL

    newTerms <- Terms

    Xdc <- model.matrix(newTerms, dataset)
    if (ncol(Xdc) > 1) Xdc <- Xdc[, -1, drop = FALSE]
    nva2 <- ncol(Xdc)


    if (!missing(newdata)){
      nsujet <- dim(newdata)[1]
    }else{
    #  nobs <- fit$n
      nsujet <- fit$groups
    }

  }
  if(class(fit) == "trivPenalNL"){
    
    
      nva3 <- fit$nvarKG
      nva4 <- fit$nvarKD

    m3 <- fit$call # longitudinal (KG)
    m3$formula <- m3$formula.terminalEvent <- m3$biomarker <- m3$formula.KD <- m3$dose <- m3$data <- m3$recurrentAG <- m3$random <- m3$id <- m3$link <- m3$n.knots <- m3$kappa <- m3$maxit <- m3$hazard <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$left.censoring <- m3$init.Random <- m3$init.Eta <- m3$init.Alpha <- m3$method.GH <- m3$n.nodes  <- m3$init.GH <- m3$time.biomarker <- m3$BoxCox <- m3$... <- NULL

    Names.data.Longi <- m3$data.Longi
    if (!missing(newdata.Longi)){m3[[3]] <- as.name(m0$newdata.Longi) # nouveau dataset
    data.Longi <- newdata.Longi
    }else{data.Longi <- eval(m3$data.Longi) }
    formula.KG <- fit$formula.KG
 
    
    m4 <- fit$call # longitudinal (KD)
    m4$formula <- m4$formula.terminalEvent <- m4$biomarker <- m4$formula.KG <- m4$dose <- m4$data <- m4$recurrentAG <- m4$random <- m4$id <- m4$link <- m4$n.knots <- m4$kappa <- m4$maxit <- m4$hazard <- m4$init.B <- m4$LIMparam <- m4$LIMlogl <- m4$LIMderiv <- m4$print.times <- m4$left.censoring <- m4$init.Random <- m4$init.Eta <- m4$init.Alpha <- m4$method.GH <- m4$n.nodes <- m4$init.GH <- m4$time.biomarker <- m4$BoxCox <- m4$... <- NULL
    
    Y <- data.Longi[,which(names(data.Longi)==fit$biomarker)]
    
    if(!is.null(formula.KG[3]) && formula.KG[3] != "1()"){
      
      TermsKG <- if (missing(data.Longi)){
        terms(formula.KG, special)
      }else{
        terms(formula.KG, special, data = data.Longi)
      }
      
      ord <- attr(TermsKG, "order") # longueur de ord=nbre de var.expli
      
      #si pas vide tous si il ya au moins un qui vaut 1 on arrete
      
      m2$formula.KG <- TermsKG
      
      
      m2[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait
      
      
      if (NROW(m3) == 0)stop("No (non-missing) observations") #nombre ligne different de 0
      
      llKG <- attr(TermsKG, "term.labels")#liste des variables explicatives
      
      
      #=========================================================>
      
      name.KG <- as.character(attr(TermsKG, "variables")[[2]])
      KG <- data.Longi[,which(names(data.Longi)==name.KG)]
      
      
      # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
      ind.placeKG <- which(lapply(as.data.frame(data.Longi[,which(names(data.Longi)%in%llKG)]),function(x) length(levels(x)))>2)
      
      defined.factor <- llKG[grep("factor",llKG)]
      
      vec.factorKG.tmp <- NULL
      if(length(defined.factor)>0){
        mat.factorKG.tmp <- matrix(defined.factor,ncol=1,nrow=length(defined.factor))
        
        # Fonction servant a prendre les termes entre "as.factor"
        vec.factorKG.tmp <-apply(mat.factorKG.tmp,MARGIN=1,FUN=function(x){
          if (length(grep("factor",x))>0){
            if(length(grep(":",x))>0){
              if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
                
                pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
                pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
                pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
                pos4 <- length(unlist(strsplit(x,split="")))
                if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2)return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
                else return(NaN)
              }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
                pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
                pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
                pos4 <- length(unlist(strsplit(x,split="")))-1
                if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos3,stop=pos4))])))>2)return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
                else return(NaN)
              }else{#both factors
                pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
                pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
                pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
                pos4 <- length(unlist(strsplit(x,split="")))-1
                if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2 || length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos3,stop=pos4))])))>2)return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
                else return(NaN)
              }
            }else{
              pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
              pos2 <- length(unlist(strsplit(x,split="")))-1
              if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2)return(substr(x,start=pos1,stop=pos2))
              else return(NaN)
            }
          }else{
            return(x)
          }})
        
        vec.factorKG.tmp <- vec.factorKG.tmp[which(vec.factorKG.tmp!="NaN")]
        
        if(length(vec.factorKG.tmp)>0){
          for(i in 1:length(vec.factorKG.tmp)){
            if(length(grep(":",vec.factorKG.tmp[i]))==0){
              if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==vec.factorKG.tmp[i])])))>2)ind.placeKG <- c(ind.placeKG,which(llKG%in%paste("as.factor(",vec.factorKG.tmp[i],")",sep="")))
            }
            
          }}
        
      }
      
      ind.placeKG <- sort(ind.placeKG)
      
      
      
      mat.factorKG2 <- matrix(llKG,ncol=1,nrow=length(llKG))
      
      # Fonction servant a prendre les termes entre "as.factor"
      llKG2 <-apply(mat.factorKG2,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          x<-substr(x,start=pos1,stop=pos2)
          return(paste(x,levels(as.factor(data.Longi[,which(names(data.Longi)==x)]))[2],sep=""))
        }else{
          return(x)
        }})
      
      # Fonction servant a prendre les termes entre "as.factor" - without the name of the level
      llKG3 <-apply(mat.factorKG2,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          return(substr(x,start=pos1,stop=pos2))
        }else{
          return(x)
        }})
      
      llKG.real.names <- llKG3  
      llKG3 <- llKG3[!llKG2%in%llKG]
      
      
      
      if(is.factor(data.Longi[,names(data.Longi)==llKG.real.names[1]]))X_L<- as.numeric(data.Longi[,names(data.Longi)==llKG.real.names[1]])-1
      else X_L<- data.Longi[,names(data.Longi)==llKG.real.names[1]]
      
      if(length(llKG) == 1)X_L <- as.data.frame(X_L)
      
      
      if(length(llKG)>1){
        
        
        for(i in 2:length(llKG.real.names)){
          
          if(is.factor(data.Longi[,names(data.Longi)==llKG.real.names[i]]))X_L<- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llKG.real.names[i]])-1)
          else X_L<- cbind(X_L,data.Longi[,names(data.Longi)==llKG.real.names[i]])
        }
      }
      #X_L<- data.Longi[,names(data.Longi)%in%(llY)]
      
      llKG.fin <- llKG.real.names
      llKG <- llKG.real.names
      
      if(sum(ord)>length(ord)){
        
        for(i in 1:length(ord)){
          if(ord[i]>1){
            
            name_v1 <- strsplit(as.character(llKG[i]),":")[[1]][1]
            name_v2 <- strsplit(as.character(llKG[i]),":")[[1]][2]
            
            if(length(grep("factor",name_v1))>0){name_v1<-substring(name_v1,11,nchar(name_v1)-1)
            v1 <- as.factor(data.Longi[,names(data.Longi)==name_v1])}
            else{v1 <- data.Longi[,names(data.Longi)==name_v1]}
            if(length(grep("factor",name_v2))>0){name_v2<-substring(name_v2,11,nchar(name_v2)-1)
            v2 <- as.factor(data.Longi[,names(data.Longi)==name_v2])}
            else{v2 <- data.Longi[,names(data.Longi)==name_v2]}
            
            llKG[i] <- paste(name_v1,":",name_v2,sep="")
            #   if(is.factor(v1) && length(levels(v1))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            #   if(is.factor(v2) && length(levels(v2))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            if(is.factor(v1) && !is.factor(v2)){
              
              dummy <- model.matrix( ~ v1 - 1)
              # if(length(levels(v1)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v1))){
                X_L <- cbind(X_L,dummy[,j]*v2)
                if(i>1 && i<length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKG.fin[(i+1+j-2):length(llKG.fin)])
                else if(i==length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""))
                else llKG.fin <- c(paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKG.fin[(2+j-2):length(llKG.fin)])
              }
              
            }else if(!is.factor(v1) && is.factor(v2)){
              
              dummy <- model.matrix( ~ v2 - 1)
              #  if(length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v2))){
                
                X_L <- cbind(X_L,dummy[,j]*v1)
                
                if(i>1 && i<length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llKG.fin[(i+1+j-2):length(llKG.fin)])
                else if(i==length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""))
                else llKG.fin <- c(paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llKG.fin[(2+j-2):length(llKG.fin)])
              }
            }else if(is.factor(v1) && is.factor(v2)){
              
              
              dummy1 <- model.matrix( ~ v1 - 1)
              dummy2 <- model.matrix( ~ v2 - 1)
              #   if(length(levels(v1)>2) || length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v1))){
                for(k in 2:length(levels(v2))){
                  
                  X_L <- cbind(X_L,dummy1[,j]*dummy2[,k])
                  if(i>1 && i<length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKG.fin[(i+1+j-2+k-2):length(llKG.fin)])
                  else if(i==length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""))
                  else llKG.fin <- c(paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKG.fin[(2+j-2+k-2):length(llKG.fin)])
                }
              } 
            }else{
              
              X_L <- cbind(X_L,v1*v2)
            }
            
          }
        }
      }
      
      
      if(length(grep(":",llKG))>0){
        for(i in 1:length(grep(":",llKG))){
          if(length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llKG[grep(":",llKG)[i]],":")[[1]])[1]]))>2 || length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llKG[grep(":",llKG)[i]],":")[[1]])[2]]))>2){
            ind.placeKG <- c(ind.placeKG,grep(":",llKG)[i])
            #     vec.factorY <- c(vec.factorY,llY[grep(":",llY)[i]])
          }
        }
      }
      vec.factorKG <- NULL
      
      if(length(vec.factorKG.tmp)>0)vec.factorKG <- c(llKG[ind.placeKG],vec.factorKG.tmp)
      else vec.factorKG <- c(vec.factorKG,llKG[ind.placeKG])
      
      vec.factorKG <- unique(vec.factorKG)
      
      
      
      mat.factorKG <- matrix(vec.factorKG,ncol=1,nrow=length(vec.factorKG))
      # Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions
      vec.factorKG <-apply(mat.factorKG,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0){
          if(length(grep(":",x))>0){
            if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
              
              pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
              pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
              pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
              pos4 <- length(unlist(strsplit(x,split="")))
              return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
            }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
              pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
              pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
              pos4 <- length(unlist(strsplit(x,split="")))-1
              return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
            }else{#both factors
              pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
              pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
              pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
              pos4 <- length(unlist(strsplit(x,split="")))-1
              return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
            }
          }else{
            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- length(unlist(strsplit(x,split="")))-1
            return(substr(x,start=pos1,stop=pos2))}
        }else{
          return(x)
        }})
      
      
      for(i in 1:length(llKG.fin)){
        
        if(sum(names(data.Longi)==llKG.fin[i])>0){
          if(is.factor(data.Longi[,names(data.Longi)==llKG.fin[i]]) && length(levels(data.Longi[,names(data.Longi)==llKG.fin[i]]))==2){
            llKG.fin[i] <- paste(llKG.fin[i],levels(data.Longi[,names(data.Longi)==llKG.fin[i]])[2],sep="")}
        }
      }
      
      #  llY <- llY.fin
      X_L <- as.data.frame(X_L)
      if(dim(X_L)[2]!=length(llKG.fin))stop("The variables in the longitudinal part must be in the data.Longi")
      
      names(X_L) <- llKG.fin
      
      
      X_Lall<- X_L
      "%+%"<- function(x,y) paste(x,y,sep="")
      if(length(vec.factorKG) > 0){
        for(i in 1:length(vec.factorKG)){
          if(length(grep(":",vec.factorKG[i]))==0){
            
            factor.spot <- which(names(X_L)==vec.factorKG[i])
            
            if(factor.spot<ncol(X_L))  X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKG[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1],X_L[(factor.spot+1):ncol(X_L)])
            else X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKG[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1])
            
          } }
        
      }
      varKG <- as.matrix(sapply(X_L, as.numeric))
      
      nsujety<-nrow(X_L)
      
    }
    
    
    #=========================================================>
    
    clusterY <- data.Longi$id
    
    max_rep <- max(table(clusterY))
    uni.clusterY<-as.factor(unique(clusterY))
    
   
    
    if(is.null(formula.KG[3]) | formula.KG[3] == "1()"){
       varKG <- c()#rep(0, dim(data.Longi)[1])
      nsujety <- length(Y)
    }
    
    formula.KD <- fit$formula.KD
    if(!is.null(formula.KD[3]) && formula.KD[3] != "1()"){
      
      TermsKD <- if (missing(data.Longi)){
        terms(formula.KD, special)
      }else{
        terms(formula.KD, special, data = data.Longi)
      }
      
      ord <- attr(TermsKD, "order") # longueur de ord=nbre de var.expli
      
      #si pas vide tous si il ya au moins un qui vaut 1 on arrete
      
      m4$formula.KD <- TermsKD
      
      
      m4[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait
      
      
      if (NROW(m4) == 0)stop("No (non-missing) observations") #nombre ligne different de 0
      
      llKD <- attr(TermsKD, "term.labels")#liste des variables explicatives
      
      
      #=========================================================>
      
      name.KD <- as.character(attr(TermsKD, "variables")[[2]])
      KD <- data.Longi[,which(names(data.Longi)==name.KD)]
      
      
      # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
      
      ind.placeKD <- which(lapply(as.data.frame(data.Longi[,which(names(data.Longi)%in%llKD)]),function(x) length(levels(x)))>2)
      
      defined.factor <- llKD[grep("factor",llKD)]
      
      vec.factorKD.tmp <- NULL
      if(length(defined.factor)>0){
        mat.factorKD.tmp <- matrix(defined.factor,ncol=1,nrow=length(defined.factor))
        
        # Fonction servant a prendre les termes entre "as.factor"
        vec.factorKD.tmp <-apply(mat.factorKG.tmp,MARGIN=1,FUN=function(x){
          if (length(grep("factor",x))>0){
            if(length(grep(":",x))>0){
              if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
                
                pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
                pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
                pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
                pos4 <- length(unlist(strsplit(x,split="")))
                if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2)return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
                else return(NaN)
              }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
                pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
                pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
                pos4 <- length(unlist(strsplit(x,split="")))-1
                if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos3,stop=pos4))])))>2)return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
                else return(NaN)
              }else{#both factors
                pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
                pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
                pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
                pos4 <- length(unlist(strsplit(x,split="")))-1
                if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2 || length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos3,stop=pos4))])))>2)return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
                else return(NaN)
              }
            }else{
              pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
              pos2 <- length(unlist(strsplit(x,split="")))-1
              if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2)return(substr(x,start=pos1,stop=pos2))
              else return(NaN)
            }
          }else{
            return(x)
          }})
        
        vec.factorKD.tmp <- vec.factorKD.tmp[which(vec.factorKD.tmp!="NaN")]
        
        if(length(vec.factorKD.tmp)>0){
          for(i in 1:length(vec.factorKD.tmp)){
            if(length(grep(":",vec.factorKD.tmp[i]))==0){
              if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==vec.factorKD.tmp[i])])))>2)ind.placeKD <- c(ind.placeKD,which(llKD%in%paste("as.factor(",vec.factorKD.tmp[i],")",sep="")))
            }
            
          }}
        
      }
      
      ind.placeKD <- sort(ind.placeKD)
      
      mat.factorKD2 <- matrix(llKD,ncol=1,nrow=length(llKD))
      
      # Fonction servant a prendre les termes entre "as.factor"
      llKD2 <-apply(mat.factorKD2,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          x<-substr(x,start=pos1,stop=pos2)
          return(paste(x,levels(as.factor(data.Longi[,which(names(data.Longi)==x)]))[2],sep=""))
        }else{
          return(x)
        }})
      
      # Fonction servant a prendre les termes entre "as.factor" - without the name of the level
      llKD3 <-apply(mat.factorKD2,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          return(substr(x,start=pos1,stop=pos2))
        }else{
          return(x)
        }})
      
      llKD.real.names <- llKD3  
      llKD3 <- llKD3[!llKD2%in%llKD]
      
      
      if(is.factor(data.Longi[,names(data.Longi)==llKD.real.names[1]]))X_L<- as.numeric(data.Longi[,names(data.Longi)==llKD.real.names[1]])-1
      else X_L<- data.Longi[,names(data.Longi)==llKD.real.names[1]]
      
      
      
      
      if(length(llKD)>1){
        for(i in 2:length(llKD.real.names)){
          
          if(is.factor(data.Longi[,names(data.Longi)==llKD.real.names[i]]))X_L<- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llKD.real.names[i]])-1)
          else X_L<- cbind(X_L,data.Longi[,names(data.Longi)==llKD.real.names[i]])
        }}
      
      #X_L<- data.Longi[,names(data.Longi)%in%(llY)]
      
      llKD.fin <- llKD.real.names
      llKD <- llKD.real.names
      
      if(sum(ord)>length(ord)){
        
        for(i in 1:length(ord)){
          if(ord[i]>1){
            
            name_v1 <- strsplit(as.character(llKD[i]),":")[[1]][1]
            name_v2 <- strsplit(as.character(llKD[i]),":")[[1]][2]
            
            if(length(grep("factor",name_v1))>0){name_v1<-substring(name_v1,11,nchar(name_v1)-1)
            v1 <- as.factor(data.Longi[,names(data.Longi)==name_v1])}
            else{v1 <- data.Longi[,names(data.Longi)==name_v1]}
            if(length(grep("factor",name_v2))>0){name_v2<-substring(name_v2,11,nchar(name_v2)-1)
            v2 <- as.factor(data.Longi[,names(data.Longi)==name_v2])}
            else{v2 <- data.Longi[,names(data.Longi)==name_v2]}
            
            llKD[i] <- paste(name_v1,":",name_v2,sep="")
            #   if(is.factor(v1) && length(levels(v1))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            #   if(is.factor(v2) && length(levels(v2))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            if(is.factor(v1) && !is.factor(v2)){
              
              dummy <- model.matrix( ~ v1 - 1)
              # if(length(levels(v1)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v1))){
                X_L <- cbind(X_L,dummy[,j]*v2)
                if(i>1 && i<length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKD.fin[(i+1+j-2):length(llKD.fin)])
                else if(i==length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""))
                else llKD.fin <- c(paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKD.fin[(2+j-2):length(llKD.fin)])
              }
              
            }else if(!is.factor(v1) && is.factor(v2)){
              
              dummy <- model.matrix( ~ v2 - 1)
              #  if(length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v2))){
                
                X_L <- cbind(X_L,dummy[,j]*v1)
                
                if(i>1 && i<length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llKD.fin[(i+1+j-2):length(llKD.fin)])
                else if(i==length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""))
                else llKD.fin <- c(paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llKD.fin[(2+j-2):length(llKD.fin)])
              }
            }else if(is.factor(v1) && is.factor(v2)){
              
              
              dummy1 <- model.matrix( ~ v1 - 1)
              dummy2 <- model.matrix( ~ v2 - 1)
              #   if(length(levels(v1)>2) || length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v1))){
                for(k in 2:length(levels(v2))){
                  
                  X_L <- cbind(X_L,dummy1[,j]*dummy2[,k])
                  if(i>1 && i<length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKG.fin[(i+1+j-2+k-2):length(llKG.fin)])
                  else if(i==length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""))
                  else llKD.fin <- c(paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKD.fin[(2+j-2+k-2):length(llKD.fin)])
                }
              } 
            }else{
              
              X_L <- cbind(X_L,v1*v2)
            }
            
          }
        }
      }
      
      
      
      if(length(grep(":",llKD))>0){
        for(i in 1:length(grep(":",llKD))){
          if(length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llKD[grep(":",llKD)[i]],":")[[1]])[1]]))>2 || length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llKD[grep(":",llKD)[i]],":")[[1]])[2]]))>2){
            ind.placeKD <- c(ind.placeKD,grep(":",llKD)[i])
            #     vec.factorY <- c(vec.factorY,llY[grep(":",llY)[i]])
          }
        }
      }
      
      vec.factorKD <- NULL
      
      if(length(vec.factorKD.tmp)>0)vec.factorKD <- c(llKD[ind.placeKD],vec.factorKD.tmp)
      else vec.factorKD <- c(vec.factorKD,llKD[ind.placeKD])
      
      vec.factorKD <- unique(vec.factorKD)
      
      
      
      mat.factorKD <- matrix(vec.factorKD,ncol=1,nrow=length(vec.factorKD))
      # Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions
      vec.factorKD <-apply(mat.factorKD,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0){
          if(length(grep(":",x))>0){
            if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
              
              pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
              pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
              pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
              pos4 <- length(unlist(strsplit(x,split="")))
              return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
            }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
              pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
              pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
              pos4 <- length(unlist(strsplit(x,split="")))-1
              return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
            }else{#both factors
              pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
              pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
              pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
              pos4 <- length(unlist(strsplit(x,split="")))-1
              return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
            }
          }else{
            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- length(unlist(strsplit(x,split="")))-1
            return(substr(x,start=pos1,stop=pos2))}
        }else{
          return(x)
        }})
      
      
      for(i in 1:length(llKD.fin)){
        
        if(sum(names(data.Longi)==llKD.fin[i])>0){
          if(is.factor(data.Longi[,names(data.Longi)==llKD.fin[i]]) && length(levels(data.Longi[,names(data.Longi)==llKD.fin[i]]))==2){
            llKD.fin[i] <- paste(llKD.fin[i],levels(data.Longi[,names(data.Longi)==llKD.fin[i]])[2],sep="")}
        }
      }
      
      #  llY <- llY.fin
      X_L <- as.data.frame(X_L)
      if(dim(X_L)[2]!=length(llKD.fin))stop("The variables in the longitudinal part must be in the data.Longi")
      
      names(X_L) <- llKD.fin
      
      
      X_Lall<- X_L
      "%+%"<- function(x,y) paste(x,y,sep="")
      if(length(vec.factorKD) > 0){
        for(i in 1:length(vec.factorKD)){
          if(length(grep(":",vec.factorKD[i]))==0){
            
            factor.spot <- which(names(X_L)==vec.factorKD[i])
            
            if(factor.spot<ncol(X_L))  X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKD[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1],X_L[(factor.spot+1):ncol(X_L)])
            else X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKD[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1])
            
          } }
        
      }
      
      varKD <- as.matrix(sapply(X_L, as.numeric))
      
      
    }else{
      varKD <- c()#rep(0, dim(data.Longi)[1])
    }
    
    matzy <- NULL       # here matzy is for time and dose
    matzy <- cbind(data.Longi[,which(colnames(data.Longi)==fit$time.biomarker)],data.Longi[,which(colnames(data.Longi)==fit$dose)])
    
    if(dim(matzy)[2] != 2)stop("Both information on the biomarker measurement times and dose must be included in the data.")
    
    
    matzy <- as.matrix(matzy)
    
    if(fit$link=="Random-effects")link <- 1
    if(fit$link=="Current-level") link <- 2
    if(fit$leftCensoring==FALSE){
      s_cag_id = 0
      s_cag = 0
    }else{
      s_cag_id = 1
      s_cag = fit$leftCensoring.threshold
    }	
    }
        
   
  
        cat("\n")
        cat("Calculating ... \n")
  if(class(fit)== 'jointPenal'){
        if(fit$logNormal==0){
        ans <- .Fortran(C_cvpl,
                        as.integer(nobs),
                        as.integer(nsujet),
                        as.integer(cluster),
                        as.integer(c),
                        as.integer(cdc),
                        as.integer(nva1),
                        as.integer(nva2),
                        as.double(X),
                        as.double(Xdc),
                        as.integer(typeof),
                        as.integer(nz),
                        as.double(zi),
                        as.double(ttt),
                        as.double(tttdc),
                        as.integer(nbintervR),
                        as.integer(nbintervDC),
                        as.integer(np),
                        as.double(b),
                        as.double(vopt),
                        as.double(tt0),
                        as.double(tt1),
                        as.double(tt0dc),
                        as.double(tt1dc),
                        as.integer(nt),
                        as.double(pred.times),
                        rl_cond=as.double(rep(0,nt)),
                        epoir=as.double(rep(0,nt)),
                        contribt=as.double(rep(0,nt*nsujet)),
                        atrisk=as.double(rep(0,nt))
		       )
}else{
#  cat('logn...')
  ans <- .Fortran(C_cvpl_logn,
                  as.integer(nobs),
                  as.integer(nsujet),
                  as.integer(cluster),
                  as.integer(c),
                  as.integer(cdc),
                  as.integer(nva1),
                  as.integer(nva2),
                  as.double(X),
                  as.double(Xdc),
                  as.integer(typeof),
                  as.integer(nz),
                  as.double(zi),
                  as.double(ttt),
                  as.double(tttdc),
                  as.integer(nbintervR),
                  as.integer(nbintervDC),
                  as.integer(np),
                  as.double(b),
                  as.double(vopt),
                  as.double(tt0),
                  as.double(tt1),
                  as.double(tt0dc),
                  as.double(tt1dc),
                  as.integer(nt),
                  as.double(pred.times),
                  rl_cond=as.double(rep(0,nt)),
                  epoir=as.double(rep(0,nt)),
                  contribt=as.double(rep(0,nt*nsujet)),
                  atrisk=as.double(rep(0,nt))
	         )
				  
  }}else if(class(fit) == "longiPenal"){

    ans <- .Fortran(C_cvpl_long,
                    as.integer(nsujet),
                    as.integer(1),
                    as.integer(length(clusterY)),
                    as.integer(0),
                    as.integer(clusterY),
                    as.integer(0),
                    as.integer(cdc),
                    as.double(yy),
                    as.integer(1),
                    as.integer(nva2),
                    as.integer(nva3),
                    as.integer(fit$ne_re),
                    as.integer(0),
                    as.integer(fit$netadc),
                    as.integer(link),
                    as.double(matrix(0,nrow=1,ncol=1)),
                    as.double(Xdc),
                    as.double(as.matrix(varY)),
                    as.double(matzy),
                    as.double(s_cag),
                    as.integer(s_cag_id),
                    as.integer(typeof),
                    as.integer(2),
                    as.integer(nz),
                    as.double(zi),
                    as.integer(np),
                    as.double(b),
                    as.double(vopt),
                    as.double(0),
                    as.double(0),
                    as.double(tt0dc),
                    as.double(tt1dc),
                    as.integer(nt),
                    as.double(pred.times),
                    rl_cond=as.double(rep(0,nt)),
                    epoir=as.double(rep(0,nt)),
                    contribt=as.double(rep(0,nt*nsujet)),
                    atrisk=as.double(rep(0,nt))
		   )

  }else if(class(fit) == "trivPenal"){


    ans <- .Fortran(C_cvpl_long,
                    as.integer(nsujet),
                    as.integer(nobs),
                    as.integer(length(clusterY)),
                    as.integer(cluster),
                    as.integer(clusterY),
                    as.integer(c),
                    as.integer(cdc),
                    as.double(yy),
                    as.integer(nva1),
                    as.integer(nva2),
                    as.integer(nva3),
                    as.integer(fit$ne_re),
                    as.integer(fit$netar),
                    as.integer(fit$netadc),
                    as.integer(link),
                    as.double(X),
                    as.double(Xdc),
                    as.double(as.matrix(varY)),
                    as.double(matzy),
                    as.double(s_cag),
                    as.integer(s_cag_id),
                    as.integer(typeof),
                    as.integer(3),
                    as.integer(nz),
                    as.double(zi),
                    as.integer(np),
                    as.double(b),
                    as.double(vopt),
                    as.double(tt0),
                    as.double(tt1),
                    as.double(tt0dc),
                    as.double(tt1dc),
                    as.integer(nt),
                    as.double(pred.times),
                    rl_cond=as.double(rep(0,nt)),
                    epoir=as.double(rep(0,nt)),
                    contribt=as.double(rep(0,nt*nsujet)),
                    atrisk=as.double(rep(0,nt))
		   )

  }else if(class(fit) == "trivPenalNL"){
    
    GH <- c(0,0)
    GH[1] <- ifelse(fit$methodGH == "Standard", 0, 1)
  #  GH[1] <- 0
    GH[2] <- fit$n.nodes
    
    box_cox <- c(0,1)
    box_cox[1] <- ifelse(fit$BoxCox == TRUE, 1, 0)
    if(!is.null(fit$BoxCox_parameter))box_cox[2] <- fit$BoxCox_parameter
    
    lappend <- function (lst, ...){
      lst <- c(lst, list(...))
      return(lst)
    }
    nodes <- gauss.quad(20,kind="hermite")$nodes
    weights <- gauss.quad(20,kind="hermite")$weights*exp(nodes^2)
    
    nodes2 <- gauss.quad(20,kind="hermite")$nodes
    weights2 <- gauss.quad(20,kind="hermite")$weights*exp(nodes2^2)
    
    tmp <- rep(list(nodes), fit$ne_re)
    tmp <- lappend(tmp, nodes2)
    
    nodes <- as.data.frame(do.call(expand.grid,tmp))
    
    tmp <- rep(list(weights), fit$ne_re)
    tmp <- lappend(tmp, weights2)
    weights <- as.data.frame(do.call(expand.grid,tmp))
    
    nodes <- sapply(nodes, as.double)
    weights <- sapply(weights, as.double)
 
 
    ans <- .Fortran(C_cvplnl,
                    as.integer(nsujet),
                    as.integer(nobs),
                    as.integer(length(clusterY)),
                    as.integer(cluster),
                    as.integer(clusterY),
                    as.integer(c),
                    as.integer(cdc),
                    as.double(Y),
                    as.integer(nva1),
                    as.integer(nva2),
                    as.integer(nva3),
                    as.integer(nva4),
                    as.integer(fit$ne_re),
                    as.integer(fit$random.which),
                    as.double(box_cox),
                    as.integer(fit$netar),
                    as.integer(fit$netadc),
                    as.integer(link),
                    as.double(X),
                    as.double(Xdc),
                    as.double(cbind(varKG,varKD)),
                    as.double(matzy),
                    as.double(s_cag),
                    as.integer(s_cag_id),
                    as.integer(typeof),
                    as.integer(nz),
                    as.double(zi),
                    as.integer(np),
                    as.double(b),
                    as.double(vopt),
                    as.double(tt0),
                    as.double(tt1),
                    as.double(tt0dc),
                    as.double(tt1dc),
                    as.integer(nt),
                    as.double(pred.times),
                    rl_cond=as.double(rep(0,nt)),
                    epoir=as.double(rep(0,nt)),
                    contribt=as.double(rep(0,nt*nsujet)),
                    atrisk=as.double(rep(0,nt)),
                    as.integer(GH),
                    as.double(fit$b_pred),
                    as.double(fit$weights),#weights),
                    as.double(fit$nodes),#nodes),#
                    as.integer(fit$n.nodes^fit$ne_re*20)#
    )
    
  }
        out <- NULL
         if (!missing(newdata)) out$data <- m0$newdata
        else out$data <- fit$data
        out$new.data <- !is.null(newdata)
        out$pred.times <- pred.times
        out$mpol <- ans$rl_cond
        if(!missing(newdata) && any(out$mpol<0))stop("The program stopped abnormally. This may be related to the new datasets with insufficient information")
        if (missing(newdata)) out$cvpol <- ans$epoir
        out$IndivContrib <- matrix(ans$contribt,nrow=nsujet,ncol=nt)
        out$AtRisk <- ans$atrisk

        cat("Estimators of EPOCE computed for",length(pred.times),"times \n")

        class(out) <- c("epoce")
        out
}
