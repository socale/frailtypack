! np - number of parameters , b- estimated parameters
! nz - number of knots, nbintervR and nbintervDC - for piecewise constant number of intervals
! nva1 - number of covariates for recurrent, nva2 - number of covariates for terminal
! nst - number of strata, typeof - type of baseline hazard estimation
! zi - vector of knots, HIHOut - inverse of hessian, time, timedc - time from piecewise constant
! ntimeAll - number of times for which we make predictions,npred0 - number of individual predictions
! predTime vector of times for prediction, window - vector of windows for predictions
! predtimerec
! nrec0 - number of recurrent events, vaxpred0 - covariates for recurrences
! vaxdcpred0 - covariates for death, predAll1, predAll2, predAll3 - empty matrices for sortie
! predAlllow/predAllhigh - empty matrices (for ci?)
! icproba - if 0 no MC, if 1 - MC, nsample - MC.sample, intcens - censure par intervalle ou pas
! trunctime -
! lowertime -                   uppertime -
! timeAll -t+window, model - typeJoint

! ============================================== pred iction Joint

    
      subroutine predict_tri(np,b,nz,nva10,nva20,nva30,nb_re0,nzyr,nzyd,link0,nst,typeof0,zi0,HIHOut, &
        ntimeAll,npred0,predTime,window,predtimerec,nrec0,nrecy0,yy0,vaxpred0,vaxdcpred0,vaxypred0,groupey,uniGroupey,nsujety, &
        nsujet,predAll1,predAlllow1,predAllhigh1, &
        icproba,nsample,movingwindow,timeAll,s_cag_id0,s_cag0)
        !model - type of a model 0- recur/survie, 1- longi/survie, 2-longi/recur/survie
        
        
        use donnees_indiv,only:nmescur,mu,ycurrent,z2,b1,X2cur,Z1cur,nmescurr
        use comon,only:etaydc,sigmae,netadc,s_cag_id,s_cag,ut,&
                    utt,nva,link,npp,nea,vey,nb1,netar,indic_Alpha,&
                    nva1,nva2, nva3,effet,zi,nz1,typeof,nb_re,alpha,&
                    etayr,typeJoint !it_cur
        use lois_normales
        use prediction
        
        implicit none    
    
        integer::i,ii,iii,j,k,it,kk,jj,npred0,nrec0,nsample,nrecy0
        integer,intent(in)::np,nz,nva10,nva20,nva30,nb_re0,nzyd,nzyr,nst,typeof0,ntimeAll,&
            icproba,movingwindow,nsujety,s_cag_id0,nsujet,link0
        double precision,dimension(nsujety),intent(in)::yy0
        double precision,dimension(npred0,nsujety)::yy_matrice
        integer,dimension(nsujety),intent(in)::groupey
        double precision,dimension(np),intent(in)::b
        double precision,dimension(nz+6),intent(in)::zi0
        double precision,dimension(np,np),intent(in)::HIHOut
        double precision,dimension(nsujet,nva10),intent(in)::vaxpred0
        double precision,dimension(npred0,nva20),intent(in)::vaxdcpred0
        double precision,dimension(nsujety,nva30),intent(in)::vaxypred0
        double precision,dimension(1,npred0)::XbetapredDC,XbetapredDCalea
        double precision,dimension(1,nsujety) :: XbetapredY,XbetapredYalea
        double precision,dimension(1,nsujet) :: XbetapredR,XbetapredRalea
        integer,dimension(npred0)::nreci,nreci_all,nrecyi, uniGroupey
        double precision::predTime,window,predTime2,scR,shR,scDC,shDC,&
            scRalea,shRalea,scDCalea,shDCalea,alea
        double precision::ss11,ss12,s_cag0
        double precision,dimension(npred0)::predProba1
        double precision,dimension(npred0,ntimeAll),intent(out)::predAll1
        double precision,dimension(npred0,ntimeAll),intent(out)::predAlllow1,predAllhigh1
        double precision,dimension(nz+2)::theR,theDC,theRalea,theDCalea
        double precision,dimension(2)::surv,survDCalea,lam
        double precision,dimension(npred0,nrec0+2)::survR,hazR,survRalea,hazRalea
        double precision,dimension(ntimeAll)::timeAll
        double precision,dimension(nsample,np)::balea
        double precision,dimension(nsample,npred0)::predProbaalea1
        double precision,dimension(1,nva30)::coefBetaYalea
        double precision,dimension(1,nva20)::coefBetadcalea
        double precision,dimension(1,nva10)::coefBetaRalea
        double precision,dimension(npred0,nrec0)::predtimerec
        double precision,dimension(npred0,nrec0+2)::predtimerec2
        double precision,dimension(1,nva20)::coefBetadc
        double precision,dimension(1,nva30)::coefBetay
        double precision,dimension(1,nva10)::coefBetaR
    
        integer :: ndim, restar,nf2
    
    
        typejoint = 3
        s_cag_id = s_cag_id0    
        s_cag = s_cag0    
        link = link0
        netadc = nzyd
        netar = nzyr    
        nea = nb_re0+1
        nb1 = nb_re0    
        nb_re = nb_re0+ INT((nb_re0*(nb_re0-1))/2.d0)
        typeof = typeof0
        indic_Alpha = 1
        nva1 = nva10
        nva2 = nva20
        nva3 = nva30
        effet = 1
  
        allocate(zi(-2:nz+3),b1(np))

        b1(1:np) = b(1:np)    
        zi(-2:nz+3) = zi0(1:nz+3+2)     
        nz1= nz
        npp=np    
        nva = nva1+nva2+nva3
        coefBetaR(1,:) = b((np-nva+1):(np-nva3-nva2))
        XbetapredR = matmul(coefBetaR,transpose(vaxpred0))    
        coefBetaDC(1,:) = b((np-nva3-nva2+1):(np-nva3))    
        XbetapredDC = matmul(coefBetadc,transpose(vaxdcpred0))
        coefBetaY(1,:) = b((np-nva3+1):np)    
        XbetapredY = matmul(coefBetay,transpose(vaxypred0))
        !XbetapredDC = matmul(coefBetadc,transpose(vaxdcpred0))    

        allocate(survRi(nrec0+2),hazRi(nrec0+2))
    
        ! Determine when to calculate Survival function (gap time)
        ! and the number of recurrences in the prediction (for each pred i)    
        ! timeAll(1) = predTime + window
        ! do i=2,ntimeAll
        ! timeAll(i) = timeAll(i-1) + window
        ! end do
    
        if (icproba.eq.1) then ! generation des parametres
            do j=1,nsample
                do i=1,np
                    call rnorm(b(i),sqrt(HIHOut(i,i)),alea)
                    balea(j,i) = alea
                end do
            end do
        end if
        k=1    
        ndim =  netar+1    
        restar = 0
        nf2 = 1    
        allocate(Ut(nea,nea),Utt(nea,nea))    
        ut = 0.d0
        utt = 0.d0
        allocate(ycurrent(nrecy0),mu(nrecy0,1),z2(nrecy0,nzyr),etaydc(netadc),etayr(netar))    
        allocate(vey(nsujety,nva3),X2cur(1,nrecy0),Z1cur(1,nb1))
        do i=1,nsujety
            do j=1,nva3
                vey(i,j) = vaxypred0(i,j)
            end do
        end do 
            
        do iii=1,ntimeAll
            nrecyi = 0
            nreci_all = 0
            k=1
            do i=1,npred0    
                if (movingwindow.eq.1) then
                    predtimerec2(i,1) = predTime
                else
                    predtimerec2(i,1) = timeAll(iii) - window
                endif                
                !do while ((groupey(k).eq.i).and.(k.lt.nsujety))
                do while ((groupey(k).eq.uniGroupey(i)).and.(k.lt.nsujety))
                    !write(*,*)npred0,i,k,vaxypred0(k,2),predtimerec2(i,1)
                    if(vaxypred0(k,2).le.predtimerec2(i,1)) then ! check if measurement happened before prediction time
                        nrecyi(i) = nrecyi(i)+1
                    end if
                    nreci_all(i) = nreci_all(i) + 1
                    if(k.eq.(nsujety-1)) then
                        if(vaxypred0(k+1,2).le.predtimerec2(i,1)) then ! check if measurement happened before prediction time
                            nrecyi(i) = nrecyi(i)+1
                        end if
                        nreci_all(i) = nreci_all(i) + 1
                    end if
                    k = k+1
                end do    
                nreci(i) = 0 
    
                do ii=2,nrec0+1    
                    if (predtimerec(i,ii-1).le.predtimerec2(i,1)) then ! check if relapse happened before prediction time
                        predtimerec2(i,ii) = predtimerec(i,ii-1)
                    else ! otherwise 0
                        predtimerec2(i,ii) = 0
                    endif
                    !if(i.eq.6)write(*,*)'66',timeAll(iii),i,ii,predtimerec2(i,:)
                    if ((ii.gt.2).and.(predtimerec2(i,ii-1).gt.0).and.(predtimerec2(i,ii).eq.0)) then
                        nreci(i) = ii-2
                    endif
                end do
                
                if (predtimerec2(i,nrec0+1).gt.0) nreci(i) = nrec0
                predtimerec2(i,nrec0+2) = timeAll(iii)!+window    
            end do
            
            !stop
            ! les y jusqu'? predtimerec2
            yy_matrice = 0.d0
            it = 1
            do kk=1,npred0    
                    if(nrecyi(kk).gt.0) then
                    yy_matrice(kk,1:nrecyi(kk)) = yy0(it:(it+nrecyi(kk)-1))    
                    end if
                    it = it+nreci_all(kk)    
            end do
    
                      
            ! Calcul des risques de base
            ! A chaque fois, calcul? pour :
            ! DC au temps de base (predtimerec2(1,1)) et ? l'horizon (predtimerec2(1,nrec0+2))
            ! Recurrence au temps de base et pour chaque temps de rechute entr? (predtimerec2(i,ii))
            ! pour chaque prediction demand?e  
    
            select case (typeof)
                case(0)
                    theR = b(1:(nz+2))*b(1:(nz+2))
                    theDC = b((nz+3):2*(nz+2))*b((nz+3):2*(nz+2))    
                    predTime2 = predtimerec2(1,1)    
                    call survival_frailty(predTime2,theR,theDC,nz+2,zi,surv,lam,nst)
                    survR(:,1) = surv(1)
                    hazR(:,1) = lam(1)
                    survDC(1) = surv(2)
                    
                    do i=1,npred0    
                        predtime_cm(1) = predtimerec2(i,1)
                        predtime_cm(2) = predtimerec2(i,nrec0+2)    
                        do ii=1,nrec0
                            predTime2 = predtimerec2(i,ii+1)
                            call survival_frailty(predTime2,theR,theDC,nz+2,zi,surv,lam,nst)
                            survR(i,ii+1) = surv(1)
                            hazR(i,ii+1) = lam(1)
                        end do
                    end do
                    predTime2 = predtimerec2(1,nrec0+2)
                    call survival_frailty(predTime2,theR,theDC,nz+2,zi,surv,lam,nst)
                    survR(:,nrec0+2) = surv(1)
                    hazR(:,nrec0+2) = lam(1)
                    survDC(2) = surv(2)    
                    
                case(2)    
                    scR = b(2)**2 !shapeweib(1)
                    shR = b(1)**2 !scaleweib(1)
                    scDC = b(4)**2 !shapeweib(2)
                    shDC = b(3)**2 !scaleweib(2)
                    
                    survR(:,1) = exp(-(predtimerec2(1,1)/scR)**shR)
                    hazR(:,1) = (shR/scR)*((predtimerec2(1,1)/scR)**(shR-1))
                    survDC(1) = exp(-(predtimerec2(1,1)/scDC)**shDC)
                    
                    do i=1,npred0
                        predtime_cm(1) = predtimerec2(i,1)
                        predtime_cm(2) = predtimerec2(i,nrec0+2)                    
                        do ii=1,nrec0+1
                            survR(i,ii+1) = exp(-(predtimerec2(i,ii+1)/scR)**shR)
                            hazR(i,ii+1) = (shR/scR)*((predtimerec2(i,ii+1)/scR)**(shR-1))
                        end do
                    end do
                    survDC(2) = exp(-(predtimerec2(1,nrec0+2)/scDC)**shDC)   
            end select
    
            alpha = b(np-nva-nb_re-1-netadc - netar)    
            sigmae = b(np-nva-nb_re)*b(np-nva-nb_re)    
    
            etaydc = b(np-nva-nb_re-netadc:np-nva-nb_re-1)
    
            etayr =b(np-nva-nb_re-netadc - netar:np-nva-nb_re-netadc -1) 
    
            Ut = 0.d0
            Utt = 0.d0
    
            do jj=1,nb1
                do k=1,jj
                Ut(jj,k)=b(np-nva-nb_re+k+jj*(jj-1)/2)
                Utt(k,jj)=b(np-nva-nb_re+k+jj*(jj-1)/2)    
                end do
            end do    
            Ut(nea,nea) =  b(np-nva-nb_re-1-netadc - netar-1)
            Utt(nea,nea) =  b(np-nva-nb_re-1-netadc - netar-1)    
  
            it = 1
            do i=1,npred0
                ycurrent  =0.d0
                mu = 0.d0
                z2 = 0.d0
                xbetapreddci=xbetapreddc(1,i)
                XbetapredRi = XbetapredR(1,i)    
                survRi = 0.d0
                hazRi= 0.d0
                survRi = survR(i,:)
                hazRi = hazR(i,:)
                nmescur = nrecyi(i)
                nmescurr = nreci(i)
                
                if(nmescur.gt.0) then    
                    ycurrent(1:nmescur) = yy_matrice(i,1:nrecyi(i))    
                    mu(1:nmescur,1) = matmul(vaxypred0(it:(it+nmescur-1),1:(nva3)),b((np-nva3+1):np))    
                    if(nzyr.eq.1) then    
                        z2(1:nmescur,1)= 1.d0
                    else    
                        z2(1:nrecyi(i),1) = 1.d0
                        z2(1:nrecyi(i),2) = vaxypred0(it:(it+nmescur-1),2)
                    end if    
                end if
                
                if(nb1.eq.1) call gauherPred_tri2(ss11,1)
                if(nb1.eq.2) call gauherPred_tri3(ss11,1) 
                
                if(nb1.eq.3) call gauherPred_tri4(ss11,1) 
                
                if(nb1.eq.1) call gauherPred_tri2(ss12,2)
                if(nb1.eq.2) call gauherPred_tri3(ss12,2)
                if(nb1.eq.3) call gauherPred_tri4(ss12,2)
             
                                
                predProba1(i) = ss11/ss12
            
                it = it +nreci_all(i)
            end do
         
            predAll1(:,iii) = predProba1    
    
            !=============================================
            ! Variabilite des proba predites
            ! Creation d'un vecteur balea, qui correspond au vecteur b o? chaque parametre
            ! est tir? au sort selon sa loi
            ! seProba1(:)=0.d0; seProba2(:)=0.d0; seProba3(:)=0.d0;seProba4(:)=0.d0;
            ! lowProba1(:)=0.d0; lowProba2(:)=0.d0; lowProba3(:)=0.d0;lowProba4(:)=0.d0;
            ! highProba1(:)=0.d0; highProba2(:)=0.d0; highProba3(:)=0.d0;highProba4(:)=0.d0;
            ! predProbaalea1(:,:)=0.d0;predProbaalea2(:,:)=0.d0;
            ! predProbaalea3(:,:)=0.d0;predProbaalea4(:,:)=0.d0;
    
            if (icproba.eq.1) then ! calcul de l'intervalle de confiance seulement si demande     
                do j=1,nsample
                    ss11 = 0.d0
                    ss12 = 0.d0    
                    XbetapredRalea = 0.d0
                    XbetapredYalea = 0.d0
                    XbetapredDCalea = 0.d0
                    survRalea = 0.d0
                    hazRalea = 0.d0
                    survDCalea = 0.d0 
                    
                    coefBetaRalea(1,:) = balea(j,(np-nva+1):(np-nva2-nva3))
                    coefBetaYalea(1,:) = balea(j,(np-nva3+1):np)
                    coefBetadcalea(1,:) = balea(j,(np-nva3-nva2+1):(np-nva3))
    
                    XbetapredRalea = matmul(coefBetaRalea,transpose(vaxpred0))
                    XbetapredYalea = matmul(coefBetaYalea,transpose(vaxypred0))
                    XbetapredDCalea = matmul(coefBetadcalea,transpose(vaxdcpred0))
    
                    select case (typeof)
                        case(0)    
                            theRalea = balea(j,1:(nz+2))*balea(j,1:(nz+2))
                            theDCalea = balea(j,(nz+3):2*(nz+2))*balea(j,(nz+3):2*(nz+2))
                            predTime2 = predtimerec2(1,1)
                            call survival_frailty(predTime2,theRalea,theDCalea,nz+2,zi,surv,lam,nst)
                            survRalea(:,1) = surv(1)
                            hazRalea(:,1) = lam(1)
                            survDCalea(1) = surv(2)    
                            do i=1,npred0
                                do ii=1,nrec0
                                    predTime2 = predtimerec2(i,ii+1)
                                    call survival_frailty(predTime2,theRalea,theDCalea,nz+2,zi,surv,lam,nst)
                                    survRalea(i,ii+1) = surv(1)
                                    hazRalea(i,ii+1) = lam(1)
                                end do
                            end do
                            predTime2 = predtimerec2(1,nrec0+2)
                            call survival_frailty(predTime2,theRalea,theDCalea,nz+2,zi,surv,lam,nst)
                            survRalea(:,nrec0+2) = surv(1)
                            hazRalea(:,nrec0+2) = lam(1)
                            survDCalea(2) = surv(2)    
                        case(2)
                            scRalea = balea(j,2)**2 !shapeweib(1)
                            shRalea = balea(j,1)**2 !scaleweib(1)
                            scDCalea = balea(j,4)**2 !shapeweib(2)
                            shDCalea = balea(j,3)**2 !scaleweib(2)
                            
                            survRalea(:,1) = exp(-(predtimerec2(1,1)/scRalea)**shRalea)
                            hazRalea(:,1) = (shRalea/scRalea)*((predtimerec2(1,1)/scRalea)**(shRalea-1))
                            survDCalea(1) = exp(-(predtimerec2(1,1)/scDCalea)**shDCalea)
                            
                            do i=1,npred0
                                do ii=1,nrec0+1
                                    survRalea(i,ii+1) = exp(-(predtimerec2(i,ii+1)/scRalea)**shRalea)
                                    hazRalea(i,ii+1) = (shRalea/scRalea)*((predtimerec2(i,ii+1)/scRalea)**(shRalea-1))
                                end do
                            end do
                            survDCalea(2) = exp(-(predtimerec2(1,nrec0+2)/scDCalea)**shDCalea)
                    end select
    
                    alpha = balea(j,np-nva-nb_re-1-netadc - netar)    
                    sigmae =balea(j,np-nva-nb_re)*balea(j,np-nva-nb_re)    
    
                    etaydc = balea(j,np-nva-nb_re-netadc:np-nva-nb_re-1)
                    etayr =balea(j,np-nva-nb_re-netadc - netar:np-nva-nb_re-netadc -1)
    
                    Ut = 0.d0
                    Utt = 0.d0
    
                    do jj=1,nb1
                        do k=1,jj
                            Ut(jj,k) = balea(j,np-nva-nb_re+k+jj*(jj-1)/2)
                            Utt(k,jj) = balea(j,np-nva-nb_re+k+jj*(jj-1)/2)                        
                        end do
                    end do     
                    Ut(nea,nea) =  balea(j,np-nva-nb_re-1-netadc - netar-1)
                    Utt(nea,nea) =  balea(j,np-nva-nb_re-1-netadc - netar-1)   
    
                    !theta = b(np-nva1-nva2-1)*b(np-nva1-nva2-1)
                    !alpha = b(np-nva1-nva2)
                    it = 1
                    do i=1,npred0
                        z2 = 0.d0
                        ycurrent = 0.d0
                        mu = 0.d0
                        survRi = 0.d0
                        hazRi= 0.d0
                        survRi = survR(i,:)
                        hazRi = hazR(i,:)
                        xbetapreddci=xbetapreddcalea(1,i)
                        XbetapredRi = XbetapredRalea(1,i)
                        nmescur = nrecyi(i)
                        nmescurr = nreci(i)
                        if(nmescur.gt.0) then    
                            ycurrent(1:nmescur) = yy_matrice(i,1:nrecyi(i))
                            mu(1:nmescur,1) = matmul(vaxypred0(it:(it+nmescur-1),1:(nva3)),balea(j,(np-nva3+1):np))    
                            if(nzyr.eq.1) then    
                                z2(1:nrecyi(i),1) = 1.d0
                            else    
                                z2(1:nrecyi(i),1) = 1.d0
                                z2(1:nrecyi(i),2) = vaxypred0(it:(it+nmescur-1),2)
                            end if    
                        end if    
                        if(nb1.eq.1) call gauherPred_tri2(ss11,1)
                        if(nb1.eq.2) call gauherPred_tri3(ss11,1)
                        if(nb1.eq.3) call gauherPred_tri4(ss11,1)
                        
                        if(nb1.eq.1) call gauherPred_tri2(ss12,2)
                        if(nb1.eq.2) call gauherPred_tri3(ss12,2)
                         if(nb1.eq.3) call gauherPred_tri4(ss11,1)
                        predProbaalea1(j,i) = ss11/ss12
                        it = it +nreci_all(i)   
                    end do  
                end do
    
                ! utilisation de la fonction percentile2 de aaUseFunction
                do i=1,npred0
                    call percentile2(predProbaalea1(:,i),nsample,predAlllow1(i,iii),predAllhigh1(i,iii))    
                end do    
            endif ! calcul de l'intervalle de confiance seulement si demande    
        end do
        
        !stop
        deallocate(mu,ycurrent,z2,survRi,hazRi,ut,utt)
        deallocate(vey,X2cur,Z1cur,b1,zi,etaydc,etayr)
        !   deallocate(ut)
        
        end subroutine predict_tri  
    
     !**************************************************************************
            !********* Gauss-Hermit pour la dimension 4 - mod?le trviarie b_10,b_11, b_12,v
            !**************************************************************************    
      SUBROUTINE gauherPred_tri4(ss,choix)    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only:typeof!,netadc,auxig
        use donnees_indiv,only : frailpol3!,numpat
        
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca
        external::gauherPred_tri2
        integer::j
        ss=0.d0
        if (typeof.eq.0) then
            do j=1,20    
                frailpol3 = x2(j)
                call gauherPred_tri3(auxfunca,choix)
                ss = ss+w2(j)*(auxfunca)    
            end do
        else   
            do j=1,32
                frailpol3 = x3(j)
                call gauherPred_tri3(auxfunca,choix)    
                ss = ss+w3(j)*(auxfunca)
            end do    
        endif
    
        return
    
        END SUBROUTINE gauherPred_tri4
        
            !**************************************************************************
            !********* Gauss-Hermit pour la dimension 3 - mod?le trviarie b_10,b_11, v
            !**************************************************************************    
      SUBROUTINE gauherPred_tri3(ss,choix)    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only:typeof!,netadc,auxig
        use donnees_indiv,only : frailpol2!,numpat
        
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca
        external::gauherPred_tri2
        integer::j
        ss=0.d0
        if (typeof.eq.0) then
            do j=1,20    
                frailpol2 = x2(j)
                call gauherPred_tri2(auxfunca,choix)
                ss = ss+w2(j)*(auxfunca)    
            end do
        else   
            do j=1,32
                frailpol2 = x3(j)
                call gauherPred_tri2(auxfunca,choix)    
                ss = ss+w3(j)*(auxfunca)
            end do    
        endif
    
        return
    
        END SUBROUTINE gauherPred_tri3
    
    !===============================================
      SUBROUTINE gauherPred_tri2(ss,choix)    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only:typeof!auxig
        use donnees_indiv,only : frailpol!,numpat
        
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca
        external::gauherPred_tri
        integer::j  

        ss=0.d0
        if (typeof.eq.0) then
            do j=1,20
                ! if (choix.eq.3) then
                frailpol = x2(j)
                call gauherPred_tri(auxfunca,choix)
                ss = ss+w2(j)*(auxfunca)
                ! endif
            end do
        else
            do j=1,32
                !  if (choix.eq.3) then
                frailpol = x3(j)
                call gauherPred_tri(auxfunca,choix)
                ss = ss+w3(j)*(auxfunca)
            end do    
        endif
    
        return
    
      END SUBROUTINE gauherPred_tri2
    
    
      SUBROUTINE gauherPred_tri(ss,choix)    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only:typeof,nb1!auxig
        use donnees_indiv,only : frailpol,frailpol2,frailpol3
        
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca,func1pred1GHtri,func2pred1GHtri,&
                          func1pred2GHtri,func2pred2GHtri,&
                          func1pred3GHtri,func2pred3GHtri
        external::func1pred1GHtri,func2pred1GHtri,func1pred2GHtri,func2pred2GHtri,&
        func1pred3GHtri,func2pred3GHtri
        integer::j
    
        ss=0.d0
        auxfunca = 0.d0
        if (typeof.eq.0) then
            do j=1,20
                if (nb1.eq.1) then
                    if(choix.eq.1) then
                        auxfunca=func1pred1GHtri(frailpol2,x2(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred1GHtri(frailpol2,x2(j))
                    endif                
                    ss = ss+w2(j)*(auxfunca)
                else if(nb1.eq.2)then
                    if(choix.eq.1) then
                        auxfunca=func1pred2GHtri(frailpol2,frailpol,x2(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred2GHtri(frailpol2,frailpol,x2(j))
                    endif
                    ss = ss+w2(j)*(auxfunca)
                 else if(nb1.eq.3)then
                    if(choix.eq.1) then
                        auxfunca=func1pred3GHtri(frailpol3,frailpol2,frailpol,x2(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred3GHtri(frailpol3,frailpol2,frailpol,x2(j))
                    endif
                    ss = ss+w2(j)*(auxfunca)
                end if
            end do
        else
            do j=1,32
                if (nb1.eq.1) then
                    if(choix.eq.1) then
                        auxfunca=func1pred1GHtri(frailpol2,x3(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred1GHtri(frailpol2,x3(j))
                    endif
                    ss = ss+w3(j)*(auxfunca)
                else if(nb1.eq.2) then
                    if(choix.eq.1) then
                        auxfunca=func1pred2GHtri(frailpol2,frailpol,x3(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred2GHtri(frailpol2,frailpol,x3(j))
                    endif
                    ss = ss+w3(j)*(auxfunca)
                    else if(nb1.eq.3) then
                    if(choix.eq.1) then
                        auxfunca=func1pred3GHtri(frailpol3,frailpol2,frailpol,x3(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred3GHtri(frailpol3,frailpol2,frailpol,x3(j))
                    endif
                    ss = ss+w3(j)*(auxfunca)
                end if
            end do
        endif   
    
        return
    
      END SUBROUTINE gauherPred_tri    
    
    !=========================
    ! Prediction  : 1 effet aleatoire
    !=========================
      double precision function func1pred1GHtri(frail,frail2)
        ! calcul de l integrant (numerateur de la fonction de prediction)
        use comon,only:etaydc,sigmae,etayr,nb1,& !netar
            s_cag,s_cag_id,alpha,ut,utt,link,npp
            !etaydc2,etayr2,netadc,nva3,vey
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1,nmescurr!x2cur,z1cur
        use prediction
        use optim
        implicit none   
    
        double precision,intent(in)::frail,frail2
        double precision,dimension(nmescur)::mu1
        double precision :: yscalar,eps,finddet,det,prod_cag
        integer :: j,jj,k,ier
        double precision,dimension((nb1+1)*(nb1+1+1)/2)::matv
        double precision,dimension((nb1+1),1)::  Xea2
        double precision,dimension((nb1+1)):: uii, Xea22
        double precision,dimension(1)::uiiui
        double precision,dimension((nb1+1),(nb1+1))::mat
        double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
        logical :: upper
        double precision :: alnorm
        double precision :: resultdc,resultR,abserr,resabs,resasc
    
        upper = .false.
    
        Xea2(1,1) = frail
        Xea2(2,1) = frail2
    
        Xea22(1) = frail
        Xea22(2) = frail2  
    
        mat = matmul(ut,utt)
    
        jj=0
        do j=1,(nb1+1)
            do k=j,(nb1+1)
                jj=j+k*(k-1)/2
                matv(jj)=mat(j,k)
            end do
        end do
        ier = 0
        eps = 1.d-10
    
        call dsinvj(matv,(nb1+1),eps,ier)
    
        mat=0.d0
        do j=1,(nb1+1)
            do k=1,(nb1+1)
                if (k.ge.j) then
                    mat(j,k)=matv(j+k*(k-1)/2)
                else
                    mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
        end do
    
        uii = matmul(Xea22,mat)
        det = finddet(matmul(ut,utt),(nb1+1))    
        uiiui=matmul(uii,Xea2)    
        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(1) = resultdc    
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(2),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(2) = resultdc
            survRi(1) = 0.d0
            ! do ii=it_rec,it_rec+nmescurr-1
            resultR = 0.d0
            call integrationdc(survRCM_pred,0.d0,predtime_cm(1),resultR,abserr,resabs,resasc,1,b1,npp,xea22)    
            survRi(1) = survRi(1) + resultR !c'est deja res1-res3
            ! end do
        end if    
    
        if(nmescur.gt.0) then    
            mu1(1:nmescur) = mu(1:nmescur,1) +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if    
        
        yscalar = 0.d0
        prod_cag = 1.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
                    prod_cag =  prod_cag*(1.d0-alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper))    
                else
                    yscalar = yscalar + (ycurrent(k)-mu1(k))**2
                end if
            end do
        else        
            do k=1,nmescur
                yscalar = yscalar + (ycurrent(k)-mu1(k))**2
            end do
        end if
    
        yscalar = dsqrt(yscalar)   
        func1pred1GHtri = 0.d0
        
        if(link.eq.1) then
            func1pred1GHtri = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+frail2*alpha)) &
                - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+frail2*alpha))) &
                * exp((frail2+etayr(1)*frail))**nmescurr &
                * (survRi(1)**( exp(XbetapredRi+frail2+etayr(1)*frail))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*dsqrt(2.d0*pi)
        else if(link.eq.2) then
            func1pred1GHtri =  (dexp(-survDC(1)*dexp(frail2*alpha))&
                -dexp(- survDC(2)*dexp(frail2*alpha))) &
                * exp((frail2+etayr(1)*frail))**nmescurr &
                * dexp(-survRi(1)* exp(frail2)) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*dsqrt(2.d0*pi)    
        end if    
    
        return
    
      end function func1pred1GHtri    
    
            !*************************************************
      double precision  function func2pred1GHtri(frail,frail2)
        ! calcul de l integrant (denominateur de la fonction de prediction)
        use optim
        use comon,only:etaydc,sigmae,etayr,s_cag,s_cag_id,nb1,& !netar
            alpha,ut,utt,link,npp!,nva3,vey,netadc,etaydc2,etayr2
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1,nmescurr!x2cur,z1cur,
        use prediction
        implicit none     
    
        double precision,intent(in)::frail,frail2
        double precision,dimension(nmescur)::mu1
        double precision :: yscalar,eps,finddet,det,prod_cag
        integer :: j,jj,k,ier
        double precision,dimension((nb1+1)*(nb1+1+1)/2)::matv
        double precision,dimension((nb1+1),1)::  Xea2
        double precision,dimension((nb1+1)):: uii, Xea22
        double precision,dimension(1)::uiiui
        double precision,dimension((nb1+1),(nb1+1))::mat
        double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
        logical :: upper
        double precision :: alnorm
        double precision :: resultdc,resultR,abserr,resabs,resasc
    
        upper = .false.    
        Xea2(1,1) = frail
        Xea2(2,1) = frail2
        Xea22(1) = frail
        Xea22(2) = frail2
    
        mat = matmul(ut,utt)
    
        jj=0
        do j=1,(nb1+1)
            do k=j,(nb1+1)
                jj=j+k*(k-1)/2
                matv(jj)=mat(j,k)
            end do
        end do
        ier = 0
        eps = 1.d-10    
        call dsinvj(matv,(nb1+1),eps,ier)    
        mat=0.d0
        do j=1,(nb1+1)
            do k=1,(nb1+1)
                if (k.ge.j) then
                    mat(j,k)=matv(j+k*(k-1)/2)
                else
                    mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
        end do
    
        uii = matmul(Xea22,mat)
        det = finddet(matmul(ut,utt),(nb1+1))    
        uiiui=matmul(uii,Xea2)   
        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(1) = resultdc    
            survRi(1) = 0.d0
            ! do ii=it_rec,it_rec+nmescurr-1
            resultR = 0.d0
            call integrationdc(survRCM_pred,0.d0,predtime_cm(1),resultR,abserr,resabs,resasc,1,b1,npp,xea22)    
            survRi(1) = survRi(1) + resultR !c'est deja res1-res3
            !  end do
        end if    

        if(nmescur.gt.0) then    
        mu1(1:nmescur) = mu(1:nmescur,1) +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if
    
        yscalar = 0.d0
        prod_cag = 1.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
                    prod_cag = prod_cag*(1.d0-alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper))   
                else
                    yscalar = yscalar + (ycurrent(k)-mu1(k))**2
                end if
            end do
        else        
            do k=1,nmescur
                yscalar = yscalar + (ycurrent(k)-mu1(k))**2
            end do
        end if
    
        !yscalar = norm2(ycurrent(1:nmescur) - mu1(1:nmescur))    
        yscalar = dsqrt(yscalar)    
        func2pred1GHtri = 0.d0
        if(link.eq.1) then
            func2pred1GHtri = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+frail2*alpha)) ) &
                * exp((frail2+etayr(1)*frail))**nmescurr &
                * (survRi(1)**( exp(XbetapredRi+frail2+etayr(1)*frail))) &
                * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*dsqrt(2.d0*pi)
        else if(link.eq.2) then
            func2pred1GHtri =  dexp(-survDC(1)*dexp(frail2*alpha)) &
                * exp((frail2+etayr(1)*frail))**nmescurr &
                * dexp(-survRi(1)* exp(frail2)) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*dsqrt(2.d0*pi)    
        end if
        
        return
    
      end function func2pred1GHtri    
    
    !=========================
    ! Prediction  : 2 effet aleatoire
    !=========================
      double precision function func1pred2GHtri(frail,frail2,frail3)
        ! calcul de l integrant (numerateur de la fonction de prediction)
        use comon,only:etaydc,sigmae,etayr,nb1,& !netar
            s_cag,s_cag_id,alpha,ut,utt,link,npp!netadc,nva3,vey
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1,nmescurr!x2cur,z1cur
        use prediction
        use optim
        implicit none  
    
        double precision,intent(in)::frail,frail2,frail3
        double precision,dimension(nmescur)::mu1
        double precision :: yscalar,eps,finddet,det,prod_cag
        integer :: j,jj,k,ier
        double precision,dimension((nb1+1)*(nb1+1+1)/2)::matv
        double precision,dimension((nb1+1),1)::  Xea2
        double precision,dimension((nb1+1)):: uii, Xea22
        double precision,dimension(1)::uiiui
        double precision,dimension((nb1+1),(nb1+1))::mat
        double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
        logical :: upper
        double precision :: alnorm
        double precision :: resultdc,resultR,abserr,resabs,resasc  

        upper = .false.
    
        Xea2(1,1) = frail
        Xea2(2,1) = frail2
        Xea2(3,1) = frail3
        Xea22(1) = frail
        Xea22(2) = frail2
        Xea22(3) = frail3
    
        mat = matmul(ut,utt)
    
        jj=0
        do j=1,(nb1+1)
            do k=j,(nb1+1)
                jj=j+k*(k-1)/2
                matv(jj)=mat(j,k)
            end do
        end do
        ier = 0
        eps = 1.d-10    
        
        call dsinvj(matv,(nb1+1),eps,ier)    
        mat=0.d0
        do j=1,(nb1+1)
            do k=1,(nb1+1)
                if (k.ge.j) then
                    mat(j,k)=matv(j+k*(k-1)/2)
                else
                    mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
        end do
    
        uii = matmul(Xea22,mat)
        det = finddet(matmul(ut,utt),(nb1+1))  
    
        uiiui=matmul(uii,Xea2)    
        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(1) = resultdc    
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(2),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(2) = resultdc
            survRi(1) = 0.d0
            ! do ii=it_rec,it_rec+nmescurr-1    
            resultR = 0.d0
            call integrationdc(survRCM_pred,0.d0,predtime_cm(1),resultR,abserr,resabs,resasc,1,b1,npp,xea22)    
            survRi(1) = survRi(1) + resultR !c'est deja res1-res3
            ! end do
        end if    
    
        if(nmescur.gt.0) then    
            mu1(1:nmescur) = mu(1:nmescur,1) +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if
    
        yscalar = 0.d0
        prod_cag = 1.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
                    prod_cag =  prod_cag*(1.d0-alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper))    
                else
                    yscalar = yscalar + (ycurrent(k)-mu1(k))**2
                end if
            end do
        else
            do k=1,nmescur
                yscalar = yscalar + (ycurrent(k)-mu1(k))**2
            end do
        end if
            
        yscalar = dsqrt(yscalar)
        func1pred2GHtri  = 0.d0
                
        if(link.eq.1) then
            func1pred2GHtri = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+etaydc(2)*frail2+frail3*alpha)) &
                - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+etaydc(2)*frail2+frail3*alpha))) &
                * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail+etayr(2)*frail2))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*(2.d0*pi)**(-3.d0/2.d0)
        else if(link.eq.2) then
            func1pred2GHtri =  (dexp(-survDC(1)*dexp(frail3*alpha))&
                -dexp(- survDC(2)*dexp(frail3*alpha))) &
                * exp(frail3)**nmescurr &
                * dexp(-survRi(1)* exp(frail3)) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*(2.d0*pi)**(-3.d0/2.d0)    
        end if
      

        return
    
        end function func1pred2GHtri    
    
        !*************************************************
      double precision  function func2pred2GHtri(frail,frail2,frail3)
        ! calcul de l integrant (denominateur de la fonction de prediction)
        use optim
        use comon,only:etaydc,sigmae,etayr,nb1,& !netar
            s_cag,s_cag_id,alpha,ut,utt,link,npp !,nva3,vey,netadc
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1,nmescurr!x2cur,z1cur    
        use prediction
        implicit none   
    
        double precision,intent(in)::frail,frail2,frail3
        double precision,dimension(nmescur)::mu1
        double precision :: yscalar,eps,finddet,det,prod_cag
        integer :: j,jj,k,ier
        double precision,dimension((nb1+1)*(nb1+1+1)/2)::matv
        double precision,dimension((nb1+1),1)::  Xea2
        double precision,dimension((nb1+1)):: uii, Xea22
        double precision,dimension(1)::uiiui
        double precision,dimension((nb1+1),(nb1+1))::mat
        double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
        logical :: upper
        double precision :: alnorm
        double precision :: resultdc,resultR,abserr,resabs,resasc
    
        upper = .false.
    
        Xea2(1,1) = frail
        Xea2(2,1) = frail2
        Xea2(3,1) = frail3
        Xea22(1) = frail
        Xea22(2) = frail2
        Xea22(3) = frail3
    
        mat = matmul(ut,utt)
    
        jj=0
        do j=1,(nb1+1)
            do k=j,(nb1+1)
                jj=j+k*(k-1)/2
                matv(jj)=mat(j,k)
            end do
        end do
        ier = 0
        eps = 1.d-10    
        call dsinvj(matv,(nb1+1),eps,ier)    
        mat=0.d0
        do j=1,(nb1+1)
            do k=1,(nb1+1)
                if (k.ge.j) then
                    mat(j,k)=matv(j+k*(k-1)/2)
                else
                    mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
        end do    
        uii = matmul(Xea22,mat)
        det = finddet(matmul(ut,utt),(nb1+1))    
        uiiui=matmul(uii,Xea2)   

        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(1) = resultdc   
            survRi(1) = 0.d0
            ! do ii=it_rec,it_rec+nmescurr-1
            resultR = 0.d0
            call integrationdc(survRCM_pred,0.d0,predtime_cm(1),resultR,abserr,resabs,resasc,1,b1,npp,xea22)    
            survRi(1) = survRi(1) + resultR !c'est deja res1-res3
            ! end do
        end if
    
        if(nmescur.gt.0) then    
             mu1(1:nmescur) = mu(1:nmescur,1) +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if
    
        yscalar = 0.d0
        prod_cag = 1.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
                    prod_cag = prod_cag*(1.d0-alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper))    
                else
                    yscalar = yscalar + (ycurrent(k)-mu1(k))**2
                end if
            end do
        else
            do k=1,nmescur
                yscalar = yscalar + (ycurrent(k)-mu1(k))**2
            end do
        end if
    
        !yscalar = norm2(ycurrent(1:nmescur) - mu1(1:nmescur))    
        yscalar = dsqrt(yscalar)   
        func2pred2GHtri = 0.d0
        if(link.eq.1) then
            func2pred2GHtri = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+etaydc(2)*frail2+frail3*alpha)) ) &
                * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail+etayr(2)*frail2))) &
                * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*(2.d0*pi)**(-3.d0/2.d0)
        else if(link.eq.2) then
            func2pred2GHtri =  dexp(-survDC(1)*dexp(frail3*alpha)) &
                * exp(frail3)**nmescurr &
                * dexp(-survRi(1)* exp(frail3)) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*(2.d0*pi)**(-3.d0/2.d0)    
        end if    
    
    
    
        return
    
      end function func2pred2GHtri    
    
    
        
    !=========================
    ! Prediction  : 3 effet aleatoire
    !=========================
      double precision function func1pred3GHtri(frail,frail2,frail3,frail4)
        ! calcul de l integrant (numerateur de la fonction de prediction)
        use comon,only:etaydc,sigmae,etayr,nb1,nea,& !netar
            s_cag,s_cag_id,alpha,ut,utt,link,npp!netadc,nva3,vey
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1,nmescurr!x2cur,z1cur
        use prediction
        use optim
        implicit none  
    
        double precision,intent(in)::frail,frail2,frail3,frail4
        double precision,dimension(nmescur)::mu1
        double precision :: yscalar,eps,finddet,det,prod_cag
        integer :: j,jj,k,ier
        double precision,dimension((nb1+1)*(nb1+1+1)/2)::matv
        double precision,dimension((nb1+1),1)::  Xea2
        double precision,dimension((nb1+1)):: uii, Xea22
        double precision,dimension(1)::uiiui
        double precision,dimension((nb1+1),(nb1+1))::mat
        double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
        logical :: upper
        double precision :: alnorm
        double precision :: resultdc,resultR,abserr,resabs,resasc  

        upper = .false.
    
        Xea2(1,1) = frail
        Xea2(2,1) = frail2
        Xea2(3,1) = frail3
        Xea2(4,1) = frail4
        Xea22(1) = frail
        Xea22(2) = frail2
        Xea22(3) = frail3
        Xea22(4) = frail4
    
        mat = matmul(ut,utt)
    
        jj=0
        do j=1,(nb1+1)
            do k=j,(nb1+1)
                jj=j+k*(k-1)/2
                matv(jj)=mat(j,k)
            end do
        end do
        ier = 0
        eps = 1.d-10    
        
        call dsinvj(matv,(nb1+1),eps,ier)    
        mat=0.d0
        do j=1,(nb1+1)
            do k=1,(nb1+1)
                if (k.ge.j) then
                    mat(j,k)=matv(j+k*(k-1)/2)
                else
                    mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
        end do
    
        uii = matmul(Xea22,mat)
        det = finddet(matmul(ut,utt),(nb1+1))  
    
        uiiui=matmul(uii,Xea2)    
        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(1) = resultdc    
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(2),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(2) = resultdc
            survRi(1) = 0.d0
            ! do ii=it_rec,it_rec+nmescurr-1    
            resultR = 0.d0
            call integrationdc(survRCM_pred,0.d0,predtime_cm(1),resultR,abserr,resabs,resasc,1,b1,npp,xea22)    
            survRi(1) = survRi(1) + resultR !c'est deja res1-res3
            ! end do
        end if    
    
        if(nmescur.gt.0) then    
            mu1(1:nmescur) = mu(1:nmescur,1) +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if
    
        yscalar = 0.d0
        prod_cag = 1.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
                    prod_cag =  prod_cag*(1.d0-alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper))    
                else
                    yscalar = yscalar + (ycurrent(k)-mu1(k))**2
                end if
            end do
        else
            do k=1,nmescur
                yscalar = yscalar + (ycurrent(k)-mu1(k))**2
            end do
        end if
            
        yscalar = dsqrt(yscalar)
        func1pred3GHtri  = 0.d0
                
        if(link.eq.1) then
            func1pred3GHtri = (survDC(1)**(exp(XbetapredDCi+dot_product(etaydc,Xea22(1:nb1))+Xea22(nea)*alpha)) &
                - survDC(2)**(exp(XbetapredDCi+dot_product(etaydc,Xea22(1:nb1))+Xea22(nea)*alpha))) &
                * exp(dot_product(etayr,Xea22(1:nb1))+Xea22(nea))**nmescurr &
                * (survRi(1)**( exp(XbetapredRi+dot_product(etayr,Xea22(1:nb1))+Xea22(nea)))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*(2.d0*pi)**(-3.d0/2.d0)
        else if(link.eq.2) then
            func1pred3GHtri =  (dexp(-survDC(1)*dexp(Xea22(nea)*alpha))&
                -dexp(- survDC(2)*dexp(Xea22(nea)*alpha))) &
                * exp(Xea22(nea))**nmescurr &
                * dexp(-survRi(1)* exp(Xea22(nea))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*(2.d0*pi)**(-3.d0/2.d0)    
        end if
      

        return
    
        end function func1pred3GHtri    
    
        !*************************************************
      double precision  function func2pred3GHtri(frail,frail2,frail3,frail4)
        ! calcul de l integrant (denominateur de la fonction de prediction)
        use optim
        use comon,only:etaydc,sigmae,etayr,nb1,nea,& !netar
            s_cag,s_cag_id,alpha,ut,utt,link,npp !,nva3,vey,netadc
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1,nmescurr!x2cur,z1cur    
        use prediction
        implicit none   
    
        double precision,intent(in)::frail,frail2,frail3,frail4
        double precision,dimension(nmescur)::mu1
        double precision :: yscalar,eps,finddet,det,prod_cag
        integer :: j,jj,k,ier
        double precision,dimension((nb1+1)*(nb1+1+1)/2)::matv
        double precision,dimension((nb1+1),1)::  Xea2
        double precision,dimension((nb1+1)):: uii, Xea22
        double precision,dimension(1)::uiiui
        double precision,dimension((nb1+1),(nb1+1))::mat
        double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
        logical :: upper
        double precision :: alnorm
        double precision :: resultdc,resultR,abserr,resabs,resasc
    
        upper = .false.
    
        Xea2(1,1) = frail
        Xea2(2,1) = frail2
        Xea2(3,1) = frail3
        Xea2(4,1) = frail4
        Xea22(1) = frail
        Xea22(2) = frail2
        Xea22(3) = frail3
        Xea22(4)  = frail4 
        
        mat = matmul(ut,utt)
    
        jj=0
        do j=1,(nb1+1)
            do k=j,(nb1+1)
                jj=j+k*(k-1)/2
                matv(jj)=mat(j,k)
            end do
        end do
        ier = 0
        eps = 1.d-10    
        call dsinvj(matv,(nb1+1),eps,ier)    
        mat=0.d0
        do j=1,(nb1+1)
            do k=1,(nb1+1)
                if (k.ge.j) then
                    mat(j,k)=matv(j+k*(k-1)/2)
                else
                    mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
        end do    
        uii = matmul(Xea22,mat)
        det = finddet(matmul(ut,utt),(nb1+1))    
        uiiui=matmul(uii,Xea2)   

        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(1) = resultdc   
            survRi(1) = 0.d0
            ! do ii=it_rec,it_rec+nmescurr-1
            resultR = 0.d0
            call integrationdc(survRCM_pred,0.d0,predtime_cm(1),resultR,abserr,resabs,resasc,1,b1,npp,xea22)    
            survRi(1) = survRi(1) + resultR !c'est deja res1-res3
            ! end do
        end if
    
        if(nmescur.gt.0) then    
             mu1(1:nmescur) = mu(1:nmescur,1) +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if
    
        yscalar = 0.d0
        prod_cag = 1.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
                    prod_cag = prod_cag*(1.d0-alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper))    
                else
                    yscalar = yscalar + (ycurrent(k)-mu1(k))**2
                end if
            end do
        else
            do k=1,nmescur
                yscalar = yscalar + (ycurrent(k)-mu1(k))**2
            end do
        end if
    
        !yscalar = norm2(ycurrent(1:nmescur) - mu1(1:nmescur))    
        yscalar = dsqrt(yscalar)   
        func2pred3GHtri = 0.d0
        if(link.eq.1) then
            func2pred3GHtri = ((survDC(1)**(exp(XbetapredDCi+dot_product(etaydc,Xea22(1:nb1))+Xea22(nea)*alpha)) ) &
                * exp(dot_product(etaydc,Xea22(1:nb1))+Xea22(nea)*alpha)**nmescurr &
                * (survRi(1)**( exp(XbetapredRi+dot_product(etayr,Xea22(1:nb1))+Xea22(nea)))) &
                * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*(2.d0*pi)**(-3.d0/2.d0)
        else if(link.eq.2) then
            func2pred3GHtri =  dexp(-survDC(1)*dexp(Xea22(nea)*alpha)) &
                * exp(Xea22(nea))**nmescurr &
                * dexp(-survRi(1)* exp(Xea22(nea))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                *dexp(-uiiui(1)/2.d0)/dsqrt(det)*(2.d0*pi)**(-3.d0/2.d0)    
        end if    
    
        return
    
      end function func2pred3GHtri    
    
    
    
    !====================================================================
      double precision function survRCM_pred(tps,it2,bh,np,frail)    
        use tailles
        use comon
        use betatttps
        use donnees_indiv
        use random_effect
        use prediction,only: XbetapredRi
    
        integer::j,np,k,n,it2,it1
        double precision::tps
        double precision,dimension(-2:np)::the2
        double precision::bbb,su
        double precision,dimension(np)::bh
        double precision,dimension(nea)::frail
        double precision,dimension(1)::current_m
    
        k=0
        j=0
        su=0.d0
        bbb=0.d0
        ! frail(1:nea) = re(1:nea)
        it1 = it2   
        X2cur(1,1) = 1.d0
        X2cur(1,2) = tps
        if(nva3.gt.2) then
            do k=3,nva3
                X2cur(1,k) = dble(vey(it_cur+1,k))    
            end do
        end if
    
        Z1cur(1,1) = 1.d0
        if(nb1.eq.2)  Z1cur(1,2) =tps    
        ! write(*,*)'frail',frail
        current_m = 1.d0
        if(nea.gt.1) then
            current_m =dot_product(X2cur(1,1:nva3),bh((np-nva3+1):np))&
                            +dot_product(Z1cur(1,1:nb1),frail(1:nb1))
        else
            current_m = dot_product(X2cur(1,1:nva3),bh((np-nva3+1):np))&
                            +Z1cur(1,1:nb1)*frail(1:nb1)
        end if    
    
        select case(typeof)
            case(0) ! calcul du risque splines    
                if(netar+netadc.ge.1) then
                    n = (np-nva-effet-indic_ALPHA-1-nb_re - netadc - netar)/(effet+1) !nst          !to znaczy ze dzielimy lliczbe wezlow na 2
                else
                    n = (np-nva-effet-indic_ALPHA-nb_re - netadc - netar)/(effet+1)
                endif   
                do k=1,n
                    ! the1(k-3)=(bh(k))**2.d0
                    the2(k-3)=(bh(k))**2.d0
                end do    
                call susps(tps,the2,nz,su,bbb,zi)   
            case(2) ! calcul du risque weibull    
                betaR = bh(1)**2
                etaR = bh(2)**2    
                if (tps.eq.0.d0) tps = 1d-12 ! utile car log(0) => -Inf
                !write(*,*)'betaD',betaD,etaD
                bbb = (betaR*dexp((betaR-1.d0)*dlog(tps))/(etaR**betaR))   
        end select
    
        survRCM_pred = bbb*XbetapredRi*dexp(etayr(1)*current_m(1))!+re(2))!+cdc(i)*etaydc1*current_mean(1)  
    
        return
    
      end function survRCM_pred 