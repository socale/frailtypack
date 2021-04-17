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

! ============================================== prediction Joint
    
        subroutine predict_biv(np,b,nz,nva20,nva30,nb_re0,nzyd,link0,nst,typeof0,zi0,HIHOut, &
        ntimeAll,npred0,predTime,window,nrec0,yy0,vaxdcpred0,vaxypred0,groupey,uniGroupey,nsujety, &
        predAll1,predAlllow1,predAllhigh1, &
        icproba,nsample,movingwindow,timeAll,s_cag_id0,s_cag0)
    
    !model - type of a model 0- recur/survie, 1- longi/survie, 2-longi/recur/survie
        use donnees_indiv,only:nmescur,mu,ycurrent,z2,b1,it_cur,X2cur,Z1cur
        use comon,only:etaydc,sigmae,netadc,s_cag_id,s_cag,ut,utt,nva,link,npp,&
        nea,vey,nb1,netar,indic_Alpha,nva1,nva2, nva3,effet,zi,nz1,typeof,nb_re,typeJoint
        use lois_normales
        use prediction
        implicit none    
    
        integer::i,iii,j,k,it,jj,npred0,nrec0,nsample
        integer,intent(in)::np,nz,nva20,nva30,nb_re0,nzyd,nst,typeof0,ntimeAll,&
                        nsujety,icproba,movingwindow,s_cag_id0,link0
        double precision,dimension(nsujety),intent(in)::yy0
        double precision,dimension(npred0,nrec0)::yy_matrice
        integer,dimension(nsujety),intent(in)::groupey
        double precision,dimension(np),intent(in)::b
        double precision,dimension(nz+6),intent(in)::zi0
        double precision,dimension(np,np),intent(in)::HIHOut
        double precision,dimension(npred0,nva20),intent(in)::vaxdcpred0
        double precision,dimension(nsujety,nva30),intent(in)::vaxypred0
        double precision,dimension(1,npred0)::XbetapredDC,XbetapredDCalea
        double precision,dimension(1,nsujety) :: XbetapredY,XbetapredYalea
        integer,dimension(npred0)::nreci,nreci_all,uniGroupey
        double precision::predTime,window,predTime2,scDC,shDC,&
                        scDCalea,shDCalea,alea
        double precision::ss11,ss12,s_cag0
        double precision,dimension(npred0)::predProba1
        double precision,dimension(npred0,ntimeAll),intent(out)::predAll1
        double precision,dimension(npred0,ntimeAll),intent(out)::predAlllow1,predAllhigh1
        double precision,dimension(nz+2)::theR,theDC,theRalea,theDCalea
        double precision,dimension(2)::surv,survDCalea,lam
        double precision,dimension(ntimeAll)::timeAll
        double precision,dimension(nsample,np)::balea
        double precision,dimension(nsample,npred0)::predProbaalea1!,predProbaalea2,predProbaalea3
        double precision,dimension(1,nva30)::coefBetaalea
        double precision,dimension(1,nva20)::coefBetadcalea
        double precision,dimension(npred0,nrec0+2)::predtimerec2
        double precision,dimension(1,nva20)::coefBetadc
        double precision,dimension(1,nva30)::coefBetay

        integer :: ndim, restar,nf2
        double precision:: epsabs,epsrel
        external :: func1pred_biv,func2pred_biv
        double precision,dimension(nzyd) :: xea
        
        typeJoint = 2
        link = link0
        netadc = nzyd
        nea = nb_re0
        nb1 = nb_re0
        typeof = typeof0
        netar = 0
        indic_Alpha = 0
        nva1 = 0
        nva2 = nva20
        nva3 = nva30
        effet = 0
        allocate(zi(-2:nz+3),b1(np))
        b1(1:np) = b(1:np)
        nb_re = nb_re0 + INT((nb_re0*(nb_re0-1))/2.d0)
        zi(-2:nz+3) = zi0(1:nz+3+2)      
        nz1= nz
        npp=np
                
        coefBetadc(1,:) = b((np-nva3-nva2+1):(np-nva3))    
        XbetapredDC = matmul(coefBetadc,transpose(vaxdcpred0))
        coefBetay(1,:) = b((np-nva3+1):np)
        XbetapredY=0.d0    
        XbetapredY = matmul(coefBetay,transpose(vaxypred0))
        
        s_cag_id = s_cag_id0
        s_cag = s_cag0
        nva = nva3+nva2
        !  XbetapredDC = matmul(coefBetadc,transpose(vaxdcpred0))
        
        allocate(vey(nsujety,nva3),X2cur(1,nrec0),Z1cur(1,nea))
        do i=1,nsujety
            do j=1,nva3
                vey(i,j) = vaxypred0(i,j)
            end do
        end do
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
        ndim =  netadc
        EPSABS=1.d-100
        EPSREL=1.d-100    
        restar = 0
        nf2 = 1  
    
        allocate(ut(nb1,nb1),utt(nb1,nb1))
        allocate(ycurrent(nrec0),mu(nrec0,1),z2(nrec0,nzyd),etaydc(netadc))
                
        do iii=1,ntimeAll
            nreci = 0
            nreci_all = 0
            k=1
            do i=1,npred0
                if (movingwindow.eq.1) then
                    predtimerec2(i,1) = predTime
                else
                    predtimerec2(i,1) = timeAll(iii) - window
                endif    
                predtimerec2(i,2) = timeAll(iii)
                predtime_cm(1:2) = predtimerec2(i,1:2)
                !do while ((groupey(k).eq.i).and.(k.lt.nsujety))    
                do while ((groupey(k).eq.uniGroupey(i)).and.(k.lt.nsujety))    
                    if(vaxypred0(k,2).le.predtimerec2(i,1)) then ! check if measurement happened before prediction time
                        nreci(i) = nreci(i)+1    
                    end if
                    nreci_all(i) = nreci_all(i) + 1    
                    if(k.eq.(nsujety-1)) then
                        if(vaxypred0(k+1,2).le.predtimerec2(i,1)) then ! check if measurement happened before prediction time
                            nreci(i) = nreci(i)+1    
                        end if
                        nreci_all(i) = nreci_all(i) + 1    
                    end if
                    k = k+1    
                end do    
            end do
    
            ! les y jusqu'? predtimerec2
            yy_matrice = 0.d0
            it = 1
            do k=1,npred0    
                if(nreci(k).gt.0) then
                    yy_matrice(k,1:nreci(k)) = yy0(it:(it+nreci(k)-1))    
                end if
                it = it+nreci_all(k)    
            end do    
    
            ! Calcul des risques de base
            ! A chaque fois, calcul? pour :
            ! DC au temps de base (predtimerec2(1,1)) et ? l'horizon (predtimerec2(1,nrec0+2))
            ! Recurrence au temps de base et pour chaque temps de rechute entr? (predtimerec2(i,ii))
            ! pour chaque prediction demand?e
                        
            if(link.eq.1) then
                select case (typeof)
                    case(0)    
                        theDC = b(1:(nz+2))*b(1:(nz+2))
                        predTime2 = predtimerec2(1,1)
                        call survival_frailty(predTime2,theR,theDC,nz+2,zi,surv,lam,nst)
                        survDC(1) = surv(2)    
                        predTime2 = predtimerec2(1,2)
                        call survival_frailty(predTime2,theR,theDC,nz+2,zi,surv,lam,nst)
                        survDC(2) = surv(2)    
                    case(2)    
                        scDC = b(2)**2 !shapeweib(2)
                        shDC = b(1)**2 !scaleweib(2)    
                        survDC(1) = exp(-(predtimerec2(1,1)/scDC)**shDC)
                        survDC(2) = exp(-(predtimerec2(1,2)/scDC)**shDC)
                end select    
            end if    

            sigmae = b(np-nva2-nva3-nb_re)*b(np-nva2-nva3-nb_re)    
            Ut = 0.d0
            Utt = 0.d0  

            do jj=1,nb1
                do k=1,jj
                    Ut(jj,k)=b(np-nva-nb_re+k+jj*(jj-1)/2)
                    Utt(k,jj)=b(np-nva-nb_re+k+jj*(jj-1)/2)    
                end do
            end do
    
            etaydc = b(np-nva-nb_re-netadc:np-nva-nb_re-1)
            
            it = 1
            do i=1,npred0
                it_cur = it
                ycurrent  =0.d0
                mu = 0.d0
                z2 = 0.d0
                xbetapreddci=xbetapreddc(1,i)
                nmescur = nreci(i)
                if(nmescur.gt.0) then                
                    ycurrent(1:nmescur) = yy_matrice(i,1:nreci(i))
                    mu(1:nmescur,1) = matmul(vaxypred0(it:(it+nmescur-1),1:(nva3)),b((np-nva3+1):np))                    
                    if(nzyd.eq.1) then                    
                        z2(1:nmescur,1)= 1.d0
                    else                    
                        z2(1:nreci(i),1) = 1.d0
                        z2(1:nreci(i),2) = vaxypred0(it:(it+nmescur-1),2)
                    end if                
                end if                
                xea = 0.d0   
                if(nb1.eq.1) then
                    call gauherPred_biv(ss11,1)
                    xea = 0.d0
                    call gauherPred_biv(ss12,2)    
                else if(nb1.eq.2) then
                    call gauherPred_biv2(ss11,1)    
                    xea = 0.d0    
                    call gauherPred_biv2(ss12,2)
                else if(nb1.eq.3) then
                    call gauherPred_biv3(ss11,1)    
                    xea = 0.d0    
                    call gauherPred_biv3(ss12,2)
                end if    
               
                
                predProba1(i) = ss11/ss12
                it = it +nreci_all(i)    
            end do

            predAll1(:,iii) = predProba1
    
    
            !=============================================
            ! Variabilite des proba predites
            ! Creation d'un vecteur balea, qui correspond au vecteur b o? chaque parametre
            ! est tir? au sort selon sa loi
    !        seProba1(:)=0.d0; seProba2(:)=0.d0; seProba3(:)=0.d0;seProba4(:)=0.d0;
    !        lowProba1(:)=0.d0; lowProba2(:)=0.d0; lowProba3(:)=0.d0;lowProba4(:)=0.d0;
    !        highProba1(:)=0.d0; highProba2(:)=0.d0; highProba3(:)=0.d0;highProba4(:)=0.d0;
    !        predProbaalea1(:,:)=0.d0;predProbaalea2(:,:)=0.d0;
    !        predProbaalea3(:,:)=0.d0;predProbaalea4(:,:)=0.d0;    
            if (icproba.eq.1) then ! calcul de l'intervalle de confiance seulement si demande    
                do j=1,nsample
                    ss11 = 0.d0
                    ss12 = 0.d0    
                    XbetapredYalea = 0.d0
                    XbetapredDCalea = 0.d0    
                    survDCalea = 0.d0    
                    coefBetaalea(1,:) = balea(j,(np-nva3+1):np)
                    coefBetadcalea(1,:) = balea(j,(np-nva2-nva3+1):(np-nva3))    
                    XbetapredYalea = matmul(coefBetaalea,transpose(vaxypred0))
                    XbetapredDCalea = matmul(coefBetadcalea,transpose(vaxdcpred0))
    
                    select case (typeof)
                        case(0)    
                            !theRalea = balea(j,1:(nz+2))*balea(j,1:(nz+2))
                            theDCalea = balea(j,(nz+3):2*(nz+2))*balea(j,(nz+3):2*(nz+2))
                            predTime2 = predtimerec2(1,1)
                            call survival_frailty(predTime2,theRalea,theDCalea,nz+2,zi,surv,lam,nst)

                            !survRalea(:,1) = surv(1)
                            !hazRalea(:,1) = lam(1)
                            survDCalea(1) = surv(2)    
                            predTime2 = predtimerec2(1,nrec0+2)
                            call survival_frailty(predTime2,theRalea,theDCalea,nz+2,zi,surv,lam,nst)

                            !survRalea(:,nrec0+2) = surv(1)
                            !hazRalea(:,nrec0+2) = lam(1)
                            survDCalea(2) = surv(2)   
                        case(2)    
                            scDCalea = balea(j,2)**2 !shapeweib(2)
                            shDCalea = balea(j,1)**2 !scaleweib(2)    
                            survDCalea(1) = exp(-(predtimerec2(1,1)/scDCalea)**shDCalea)   
                            survDCalea(2) = exp(-(predtimerec2(1,nrec0+2)/scDCalea)**shDCalea)
                    end select    
    
                    sigmae = balea(j,np-nva2-nva3-nb_re)*balea(j,np-nva2-nva3-nb_re)    
                    Ut = 0.d0
                    Utt = 0.d0    
                    do jj=1,nb1
                        do k=1,jj
                            Ut(jj,k)=balea(j,np-nva-nb_re+k+jj*(jj-1)/2)
                            Utt(k,jj)=balea(j,np-nva-nb_re+k+jj*(jj-1)/2)    
                        end do
                    end do
    
                    etaydc = balea(j,np-nva-nb_re-netadc:np-nva-nb_re-1) 
                    it = 1
                    do i=1,npred0
                        z2 = 0.d0
                        ycurrent = 0.d0
                        mu = 0.d0    
                        xbetapreddci=xbetapreddcalea(1,i)
                        nmescur = nreci(i)
                        if(nmescur.gt.0) then    
                            ycurrent(1:nmescur) = yy_matrice(i,1:nreci(i))
                            mu(1:nmescur,1) = matmul(vaxypred0(it:(it+nmescur),1:(nva3)),balea(j,(np-nva3+1):np))
                            if(nzyd.eq.1) then    
                                z2(1:nreci(i),1) = 1.d0
                            else    
                                z2(1:nreci(i),1) = 1.d0
                                z2(1:nreci(i),2) = vaxypred0(it:(it+nmescur-1),2)
                            end if    
                        end if  

                        if(nb1.eq.1) then
                    call gauherPred_biv(ss11,1)
                    xea = 0.d0
                    call gauherPred_biv(ss12,2)    
                else if(nb1.eq.2) then
                    call gauherPred_biv2(ss11,1)    
                    xea = 0.d0    
                    call gauherPred_biv2(ss12,2)
                else if(nb1.eq.3) then
                    call gauherPred_biv3(ss11,1)    
                    xea = 0.d0    
                    call gauherPred_biv3(ss12,2)
                end if      
                       
                        predProbaalea1(j,i)= ss11/ss12
                        it = it +nreci_all(i)   
                    end do    
                end do
    
                ! utilisation de la fonction percentile2 de aaUseFunction
                do i=1,npred0
                    call percentile2(predProbaalea1(:,i),nsample,predAlllow1(i,iii),predAllhigh1(i,iii))
                !call percentile2(predProbaalea2(:,i),nsample,predAlllow2(i,iii),predAllhigh2(i,iii))
                !call percentile2(predProbaalea3(:,i),nsample,predAlllow3(i,iii),predAllhigh3(i,iii))
                end do    
            endif ! calcul de l'intervalle de confiance seulement si demande    
        end do
    
        deallocate(mu,ycurrent,z2,ut,utt)
        deallocate(vey,X2cur,Z1cur,b1,zi,etaydc)
        !deallocate(ut)

        end subroutine predict_biv    
    
    !=========================
    ! Prediction 1 : exactement j recurrences
    !=========================
        subroutine func1pred_biv(ndim2,xea,nf2,funvls)
        ! calcul de l integrant (numerateur de la fonction de prediction)
        use comon,only:nb1,etaydc,sigmae,s_cag_id,s_cag,ut
        use donnees_indiv,only:nmescur,mu,z2,ycurrent
        use prediction
    
        implicit none    
    
        double precision,dimension(nb1):: xea,xea2,ui
        integer :: nf2,j,ndim2,k
        double precision::yscalar,funvls,vraisind,prod_cag
        double precision,dimension(nmescur)::mu1

        ndim2 = nb1
        nf2 = 1
        vraisind = 0.d0
    
        Xea2=0.d0
        do j=1,nb1
            Xea2(j)=Xea(j)
        end do  
    
        ui = MATMUL(Ut,Xea2)    
        mu1 = mu(1:nmescur,1) + MATMUL(Z2(1:nmescur,1:nb1),ui(1:nb1))
        yscalar= 0.d0
        prod_cag = 0.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
                    
                    prod_cag = prod_cag+dlog(0.5*(1+erf((-mu1(k)+s_cag)/(sigmae*dsqrt(2.d0))))) 
               
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
           vraisind = ((survDC(1)**(exp(XbetapredDCi+dot_product(etaydc,ui(1:nb1)))) &
                - survDC(2)**(exp(XbetapredDCi+dot_product(etaydc,ui(1:nb1))))) &
                * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag)    
    
    
        funvls = vraisind
    
        end subroutine func1pred_biv
    
        subroutine func2pred_biv(ndim2,xea,nf2,funvls)
        ! calcul de l integrant (denominateur de la fonction de prediction)

        use comon,only:nb1,etaydc,sigmae,s_cag_id,s_cag,ut
        use donnees_indiv,only:nmescur,mu,z2,ycurrent
        use prediction
        implicit none   
    
        double precision,dimension(nb1):: xea,xea2,ui
        integer :: nf2,j,ndim2,k
        double precision::funvls,yscalar,prod_cag
        double precision,dimension(nmescur)::mu1
    
        ndim2 =nb1
        nf2 = 1    
    
        Xea2=0.d0
        do j=1,nb1
            Xea2(j)=Xea(j)
        end do
    
        ui = MATMUL(Ut,Xea2)    
        mu1 = mu(1:nmescur,1) + MATMUL(Z2(1:nmescur,1:nb1),ui(1:nb1))   
        yscalar = 0.d0
        prod_cag = 0.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
             
                    prod_cag = prod_cag+dlog(0.5*(1+erf((-mu1(k)+s_cag)/(sigmae*dsqrt(2.d0)))))             !dlog(alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper)) !
                    ! mu1(k) = ycurrent(k)
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
    
              funvls = ((survDC(1)**(exp(XbetapredDCi+dot_product(etaydc,ui(1:nb1))))) &
                * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag)    
       
    
        end subroutine func2pred_biv
        
             !===============================================
        SUBROUTINE gauherPred_biv3(ss,choix)
    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only:typeof!auxig
        use donnees_indiv,only : frailpol2!,numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca
        external::gauherPred_biv
        integer::j   
    
        ss=0.d0
        if (typeof.eq.0) then
            do j=1,20
                !if (choix.eq.3) then
                frailpol2 = x2(j)
                call gauherPred_biv2(auxfunca,choix)
                ss = ss+w2(j)*(auxfunca)
                !endif
            end do
        else
            do j=1,32
                !  if (choix.eq.3) then
                frailpol2 = x3(j)
                call gauherPred_biv2(auxfunca,choix)
                ss = ss+w3(j)*(auxfunca)
            end do    
        endif
    
        return
    
        END SUBROUTINE gauherPred_biv3   
    
        !===============================================
        SUBROUTINE gauherPred_biv2(ss,choix)
    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only:typeof!auxig
        use donnees_indiv,only : frailpol!,numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca
        external::gauherPred_biv
        integer::j   
    
        ss=0.d0
        if (typeof.eq.0) then
            do j=1,20
                !if (choix.eq.3) then
                frailpol = x2(j)
                call gauherPred_biv(auxfunca,choix)
                ss = ss+w2(j)*(auxfunca)
                !endif
            end do
        else
            do j=1,32
                !  if (choix.eq.3) then
                frailpol = x3(j)
                call gauherPred_biv(auxfunca,choix)
                ss = ss+w3(j)*(auxfunca)
            end do    
        endif
    
        return
    
        END SUBROUTINE gauherPred_biv2    
    
        SUBROUTINE gauherPred_biv(ss,choix)
    
        use tailles
        use donnees,only:x2,w2,x3!,w3
        use comon,only:typeof,netadc!auxig
        use donnees_indiv,only : frailpol,frailpol2
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca,func1pred_bivGH1,func2pred_bivGH1,&
        func1pred_bivGH2,func2pred_bivGH2,func1pred_bivGH3,func2pred_bivGH3
        external::func1pred_bivGH1,func2pred_bivGH1,func1pred_bivGH2,func2pred_bivGH2,&
        func1pred_bivGH3,func2pred_bivGH3
        integer::j
        
        auxfunca = 0.d0
        ss=0.d0
        if (typeof.eq.0) then        
            do j=1,20
                if (netadc.eq.3) then
                    if(choix.eq.1) then                   
                        auxfunca=func1pred_bivGH3(frailpol2,frailpol,x2(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred_bivGH3(frailpol2,frailpol,x2(j))
                    endif
                    ss = ss+w2(j)*(auxfunca)
               
                else if (netadc.eq.2) then
                    if(choix.eq.1) then                   
                        auxfunca=func1pred_bivGH2(frailpol,x2(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred_bivGH2(frailpol,x2(j))
                    endif
                    ss = ss+w2(j)*(auxfunca)
                else
                    if(choix.eq.1) then
                        auxfunca=func1pred_bivGH1(x2(j))
                    else if(choix.eq.2) then    
                        auxfunca = func2pred_bivGH1(x2(j))
                    endif
                    ss = ss+w2(j)*(auxfunca)
                endif
            end do
        else
            do j=1,32
                if (netadc.eq.3) then
                    if(choix.eq.1) then    
                        auxfunca=func1pred_bivGH3(frailpol2,frailpol,x3(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred_bivGH3(frailpol2,frailpol,x3(j))
                    endif
                    ss = ss+w2(j)*(auxfunca)
                else if (netadc.eq.2) then
                    if(choix.eq.1) then    
                        auxfunca=func1pred_bivGH2(frailpol,x3(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred_bivGH2(frailpol,x3(j))
                    endif
                    ss = ss+w2(j)*(auxfunca)
                else
                    if(choix.eq.1) then
                        auxfunca=func1pred_bivGH1(x3(j))
                    else if(choix.eq.2) then
                        auxfunca = func2pred_bivGH1(x3(j))
                    endif
                    ss = ss+w2(j)*(auxfunca)
                endif
            end do
        endif
    
        !end if
    
        return
    
        END SUBROUTINE gauherPred_biv   
        
        
        !=========================
    ! Prediction  : 3 effets        aleatoires
    !=========================
        double precision function  func1pred_bivGH3(frail,frail2,frail3)
        ! calcul de l integrant (numerateur de la fonction de prediction)
        use comon,only:nb1,etaydc,sigmae,s_cag_id,s_cag,&
        ut,utt,link,npp !nva3,vey
    
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1
        use prediction
        use optim
        implicit none    
    
        !double precision::XbetapredDCi
        !double precision,dimension(2)::survDC
        double precision,intent(in)::frail,frail2,frail3
        !double precision,dimension(netadc):: xea,xea2,ui
        double precision :: eps,finddet,det,prod_cag,alnorm
        integer :: j,k
        double precision::yscalar
        double precision,dimension(nmescur)::mu1
        integer :: jj,ier
        double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1):: Xea2
        double precision,dimension(nb1):: uii, Xea22
        double precision,dimension(1)::uiiui
        double precision,dimension(nb1,nb1)::mat
        double precision,external::survdcCM_pred
        double precision :: resultdc,abserr,resabs,resasc
        logical :: upper
        double precision,parameter::pi=3.141592653589793d0    
        
        upper = .false.      
        Xea2=0.d0
        Xea2(1,1) = frail
        Xea2(2,1) = frail2
        Xea2(3,1) = frail3
        Xea22(1) = frail
        Xea22(2) = frail2 
        Xea22(3) = frail3
        mat = matmul(ut,utt)    
        jj=0
        ! jjj = 0

        do j=1,nb1
            do k=j,nb1
                jj=j+k*(k-1)/2
                !jjj = jjj +1
                matv(jj)=mat(j,k)
                
            end do
        end do
        ier = 0
        eps = 1.d-10 

        call dsinvj(matv,nb1,eps,ier)    
        mat=0.d0

        do j=1,nb1
            do k=1,nb1
                if (k.ge.j) then
                    mat(j,k)=matv(j+k*(k-1)/2)
                else
                    mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
        end do    
    
        uii = matmul(Xea22,mat)
        det = finddet(matmul(ut,utt),nb1)
          
        uiiui=matmul(uii,Xea2)

        if(nmescur.gt.0) then
            mu1(1:nmescur) = mu(1:nmescur,1) +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if   
    
        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(1) = resultdc    
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(2),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(2) = resultdc    
        end if    
        yscalar= 0.d0
        prod_cag = 1.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
                   prod_cag =prod_cag*(1.d0-alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper))
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
        if(link.eq.1)  then
            func1pred_bivGH3 = (survDC(1)**exp(XbetapredDCi+dot_product(etaydc,Xea22(1:nb1))) &
                - survDC(2)**exp(XbetapredDCi+dot_product(etaydc,Xea22(1:nb1)))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*prod_cag *dexp(-uiiui(1)/2.d0)&
                *1/(dsqrt(det)*2.d0*pi)
        else                
            func1pred_bivGH3 = (dexp(-survDC(1))-dexp(-survDC(2))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*prod_cag *dexp(-uiiui(1)/2.d0)&
                *1/(dsqrt(det)*2.d0*pi)
        end if  

        return
    
        end function func1pred_bivGH3
    
        double precision function  func2pred_bivGH3(frail,frail2,frail3)
        ! calcul de l integrant (numerateur de la fonction de prediction)
        use comon,only:nb1,etaydc,sigmae,s_cag_id,s_cag,&
                  ut,utt,link,npp !nva3,vey
    
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1
        use prediction
        use optim
        implicit none  
    
        double precision,intent(in)::frail,frail2,frail3
        !double precision,dimension(netadc):: xea,xea2,ui
        integer :: j,k
        double precision::yscalar,prod_cag,eps,det,finddet,alnorm
        double precision,dimension(nmescur)::mu1
        integer :: jj,ier
        double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22
        double precision,dimension(1)::uiiui
        double precision,dimension(nb1,nb1)::mat
        double precision,external::survdcCM_pred
        double precision :: resultdc,abserr,resabs,resasc
        logical :: upper
        double precision,parameter::pi=3.141592653589793d0

        upper = .false.    
        Xea2=0.d0
        Xea2(1,1) = frail
        Xea2(2,1) = frail2
        Xea2(3,1) = frail3
        Xea22(1) = frail
        Xea22(2) = frail2  
        Xea22(3) = frail3
        mat = matmul(ut,utt)

        jj=0
        ! jjj = 0
        do j=1,nb1
            do k=j,nb1
                jj=j+k*(k-1)/2
                !jjj = jjj +1
                matv(jj)=mat(j,k)
                !bb2vv(jjj)=bb2(j,k)   
            end do
        end do
        ier = 0
        eps = 1.d-10    
        call dsinvj(matv,nb1,eps,ier)    
        mat=0.d0
        do j=1,nb1
            do k=1,nb1
                if (k.ge.j) then
                    mat(j,k)=matv(j+k*(k-1)/2)
                else
                    mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
        end do    
        uii = matmul(Xea22,mat)
        det = finddet(matmul(ut,utt),nb1)
          
        uiiui=matmul(uii,Xea2)   

        if(nmescur.gt.0) then
             mu1(1:nmescur) = mu(1:nmescur,1) +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if
        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(1) = resultdc    
        end if 
    
        yscalar= 0.d0
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
    
        yscalar = dsqrt(yscalar)
    
        if(link.eq.1) then
            func2pred_bivGH3 = ((survDC(1)**(exp(XbetapredDCi+dot_product(etaydc,Xea22(1:nb1)))) ) &
                * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*prod_cag*dexp(-uiiui(1)/2.d0)&
                *1/(dsqrt(det)*2.d0*pi)
        else
            func2pred_bivGH3 = dexp(-survDC(1))*dexp(-(yscalar**2.d0)/(2.d0*sigmae))*prod_cag&
                *dexp(-uiiui(1)/2.d0)*1/(dsqrt(det)*2.d0*pi)
        end if

        return
    
        end function func2pred_bivGH3
    
    !=========================
    ! Prediction  : 2 effets        aleatoires
    !=========================
        double precision function  func1pred_bivGH2(frail,frail2)
        ! calcul de l integrant (numerateur de la fonction de prediction)
        use comon,only:nb1,etaydc,sigmae,s_cag_id,s_cag,&
        ut,utt,link,npp !nva3,vey
    
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1
        use prediction
        use optim
        implicit none    
    
        !double precision::XbetapredDCi
        !double precision,dimension(2)::survDC
        double precision,intent(in)::frail,frail2
        !double precision,dimension(netadc):: xea,xea2,ui
        double precision :: eps,finddet,det,prod_cag,alnorm
        integer :: j,k
        double precision::yscalar
        double precision,dimension(nmescur)::mu1
        integer :: jj,ier
        double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1):: Xea2
        double precision,dimension(nb1):: uii, Xea22
        double precision,dimension(1)::uiiui
        double precision,dimension(nb1,nb1)::mat
        double precision,external::survdcCM_pred
        double precision :: resultdc,abserr,resabs,resasc
        logical :: upper
        double precision,parameter::pi=3.141592653589793d0    
        
        upper = .false.      
        Xea2=0.d0
        Xea2(1,1) = frail
        Xea2(2,1) = frail2
        Xea22(1) = frail
        Xea22(2) = frail2    
        mat = matmul(ut,utt)    
        jj=0
        ! jjj = 0

        do j=1,nb1
            do k=j,nb1
                jj=j+k*(k-1)/2
                !jjj = jjj +1
                matv(jj)=mat(j,k)
                ! bb2vv(jjj)=bb2(j,k)    
            end do
        end do
        ier = 0
        eps = 1.d-10 

        call dsinvj(matv,nb1,eps,ier)    
        mat=0.d0

        do j=1,nb1
            do k=1,nb1
                if (k.ge.j) then
                    mat(j,k)=matv(j+k*(k-1)/2)
                else
                    mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
        end do    
    
        uii = matmul(Xea22,mat)
        det = finddet(matmul(ut,utt),nb1)
        
        uiiui=matmul(uii,Xea2)

        if(nmescur.gt.0) then
            mu1(1:nmescur) = mu(1:nmescur,1) +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if   
    
        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(1) = resultdc    
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(2),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(2) = resultdc    
        end if    
        yscalar= 0.d0
        prod_cag = 1.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
                   prod_cag =prod_cag*(1.d0-alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper))
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
        if(link.eq.1)  then
            func1pred_bivGH2 = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2))) &
                - survDC(2)**(exp(XbetapredDCi+etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2)))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*prod_cag *dexp(-uiiui(1)/2.d0)&
                *1/(dsqrt(det)*2.d0*pi)
        else                
            func1pred_bivGH2 = (dexp(-survDC(1))-dexp(-survDC(2))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*prod_cag *dexp(-uiiui(1)/2.d0)&
                *1/(dsqrt(det)*2.d0*pi)
        end if  

        return
    
        end function func1pred_bivGH2
    
        double precision function  func2pred_bivGH2(frail,frail2)
        ! calcul de l integrant (numerateur de la fonction de prediction)
        use comon,only:nb1,etaydc,sigmae,s_cag_id,s_cag,&
                  ut,utt,link,npp !nva3,vey
    
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1
        use prediction
        use optim
        implicit none  
    
        double precision,intent(in)::frail,frail2
        !double precision,dimension(netadc):: xea,xea2,ui
        integer :: j,k
        double precision::yscalar,prod_cag,eps,det,finddet,alnorm
        double precision,dimension(nmescur)::mu1
        integer :: jj,ier
        double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22
        double precision,dimension(1)::uiiui
        double precision,dimension(nb1,nb1)::mat
        double precision,external::survdcCM_pred
        double precision :: resultdc,abserr,resabs,resasc
        logical :: upper
        double precision,parameter::pi=3.141592653589793d0

        upper = .false.    
        Xea2=0.d0
        Xea2(1,1) = frail
        Xea2(2,1) = frail2
        Xea22(1) = frail
        Xea22(2) = frail2    
        mat = matmul(ut,utt)

        jj=0
        ! jjj = 0
        do j=1,nb1
            do k=j,nb1
                jj=j+k*(k-1)/2
                !jjj = jjj +1
                matv(jj)=mat(j,k)
                !bb2vv(jjj)=bb2(j,k)   
            end do
        end do
        ier = 0
        eps = 1.d-10    
        call dsinvj(matv,nb1,eps,ier)    
        mat=0.d0
        do j=1,nb1
            do k=1,nb1
                if (k.ge.j) then
                    mat(j,k)=matv(j+k*(k-1)/2)
                else
                    mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
        end do    
        uii = matmul(Xea22,mat)
        det = finddet(matmul(ut,utt),nb1)
        uiiui=matmul(uii,Xea2)   

        if(nmescur.gt.0) then
            mu1(1:nmescur) = mu(1:nmescur,1) +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if
         if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,Xea22)
            survDC(1) = resultdc    
        end if 
    
        yscalar= 0.d0
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
    
        yscalar = dsqrt(yscalar)
    
        if(link.eq.1) then
            func2pred_bivGH2 = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2))) ) &
                * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*prod_cag*dexp(-uiiui(1)/2.d0)&
                *1/(dsqrt(det)*2.d0*pi)
        else
            func2pred_bivGH2 = dexp(-survDC(1))*dexp(-(yscalar**2.d0)/(2.d0*sigmae))*prod_cag&
                *dexp(-uiiui(1)/2.d0)*1/(dsqrt(det)*2.d0*pi)
        end if

        return
    
        end function func2pred_bivGH2   
    
    !=========================
    ! Prediction  : 1 effet aleatoire
    !=========================
        double precision function  func1pred_bivGH1(frail1)
        ! calcul de l integrant (numerateur de la fonction de prediction)
        use comon,only:etaydc,sigmae,s_cag_id,&
                s_cag,ut,link,npp !nva3,vey,utt,etaydc2,netadc
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1!,it_cur,x2cur,z1cur
        use prediction
        use optim
        implicit none    
        
		! ==modification SCL (frail --> frail1) 10/04/2020 pour correction Rank mismatch. 
		! en effet, la subroutine integrationdc() qui appelle cette variable attend un 
		! vecteur et pas un scalaire. du coup, je renomme la variable pour creer le vecteur 
		! par la suite. ceci m'evite de modivier la definition de la subroutine ===
         double precision,intent(in)::frail1 
        ! double precision,dimension(netadc):: xea,xea2,ui
        double precision ::prod_cag,alnorm
        integer :: k
        double precision::yscalar
        double precision,dimension(:),allocatable::mu1
        double precision,dimension(:),allocatable::frail ! /* ajout scl 10/04/2020 pour correction Rank mismatch */
        logical :: upper
        double precision,external::survdcCM_pred
        double precision :: resultdc,abserr,resabs,resasc
        double precision,parameter::pi=3.141592653589793d0
        
		! =========== ajout scl 10/04/2020 pour correction Rank mismatch ========
		allocate(frail(1))
		frail(1) = frail1
		! =========== Fin ajout scl ================
		
        upper = .false.
    
        if(nmescur.gt.0) then
            allocate(mu1(nmescur))
        else
            allocate(mu1(1))
        end if
    
        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,frail)
            survDC(1) = resultdc    
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(2),resultdc,abserr,resabs,resasc,1,b1,npp,frail)
            survDC(2) = resultdc    
        end if  

        if(nmescur.gt.0) then
            mu1(1:nmescur) = mu(1:nmescur,1) +frail(1)*Z2(1:nmescur,1) !scl
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if
    
        yscalar= 0.d0
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

        yscalar = dsqrt(yscalar) 
        if(link.eq.1)  then    
            func1pred_bivGH1 = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail(1))) & 
                - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail(1)))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*prod_cag &
                *dexp( - (frail(1)**2.d0)/(2.d0*ut(1,1)**2))&
                *1/dsqrt(ut(1,1)*2.d0*pi)
        else
            func1pred_bivGH1 = (dexp(-survDC(1))-dexp(- survDC(2))) &
                * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*prod_cag&
                *dexp( - (frail(1)**2.d0)/(2.d0*ut(1,1)**2))&
                *1/dsqrt(ut(1,1)*2.d0*pi)
        end if
        
        deallocate(mu1, frail)

        return
    
        end function func1pred_bivGH1
    
    
        double precision function  func2pred_bivGH1(frail)
        ! calcul de l integrant (numerateur de la fonction de prediction)
        use comon,only:etaydc,sigmae,s_cag_id,s_cag,&
               ut,link,npp,nb1!,nva3,vey,etaydc2,netadc,utt
        use donnees_indiv,only:nmescur,mu,z2,ycurrent,b1!,x2cur,z1cur,it_cur
        use prediction
        use optim
        
        implicit none    
        
        double precision,dimension(nb1),intent(in)::frail
        integer :: k
        double precision::yscalar,prod_cag,alnorm
        double precision,dimension(:),allocatable::mu1
        logical :: upper
        double precision,external::survdcCM_pred
        double precision :: resultdc,abserr,resabs,resasc
        double precision,parameter::pi=3.141592653589793d0
    
        upper = .false.
    
        if(nmescur.gt.0) then
            allocate(mu1(nmescur))
        else
            allocate(mu1(1))
        end if
    
        if(link.eq.2) then
            call integrationdc(survdcCM_pred,0.d0,predtime_cm(1),resultdc,abserr,resabs,resasc,1,b1,npp,frail)
            survDC(1) = resultdc    
        end if
    
        if(nmescur.gt.0) then
            mu1(1:nmescur) = mu(1:nmescur,1) +frail(1)*Z2(1:nmescur,1)
        else
            mu1(1:nmescur)  = mu(1:nmescur,1)
        end if
    
        yscalar= 0.d0
        prod_cag = 1.d0
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then    
                     prod_cag = prod_cag*(1.d0-alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper))            !dlog(alnorm((mu1(k)-s_cag)/sqrt(sigmae),upper)) !
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
        if(link.eq.1) then    
            func2pred_bivGH1 = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail(1))) ) &
                * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*prod_cag&
                *dexp( - (frail(1)**2.d0)/(2.d0*ut(1,1)**2))&
                *1/dsqrt(ut(1,1)*2.d0*pi)
        else    
            func2pred_bivGH1 = dexp(-survDC(1)-(yscalar**2.d0)/(2.d0*sigmae))*prod_cag&
                *dexp( - (frail(1)**2.d0)/(2.d0*ut(1,1)**2))&
                *1/dsqrt(ut(1,1)*2.d0*pi)
        end if    
    
        deallocate(mu1)
        return
    
        end function func2pred_bivGH1
    
            !========================================
            !========================================
    
        double precision function survdcCM_pred(tps,i,bh,np,frail)
    
        use tailles
        use comon
        use betatttps
        use donnees_indiv
        use prediction,only: XbetapredDCi
    
        integer::j,i,np,k,n,numpa
        double precision::tps
        double precision,dimension(-2:np)::the2
        double precision::bbb,su
        double precision,dimension(1)::current_m
        double precision,dimension(np)::bh
        double precision,dimension(nb1)::frail    
    
        k=0
        j=0
        su=0.d0
        bbb=0.d0
        numpa = i    
        X2cur(1,1) = 1.d0
        X2cur(1,2) = tps
        
        if(nva3.gt.2) then
            do k=3,nva3
                X2cur(1,k) = dble(vey(it_cur+1,k))    
            end do
        end if
    
        Z1cur(1,1) = 1.d0
        if(nb1.eq.2)  Z1cur(1,2) = tps    
        current_m = 1.d0
        if(nea.gt.1) then
            current_m =dot_product(X2cur(1,1:nva3),bh((np-nva3+1):np))&
                    +dot_product(Z1cur(1,1:nb1),frail(1:nb1))
        else
            current_m= dot_product(X2cur(1,1:nva3),bh((np-nva3+1):np))+Z1cur(1,1:nb1)*frail(1:nb1)
        end if
    
        select case(typeof)
            case(0) ! calcul du risque splines    
                if(netar+netadc.ge.1) then
                    n = (np-nva-effet-indic_ALPHA-1-nb1 - netadc - netar)/(effet+1) !nst        !to znaczy ze dzielimy lliczbe wezlow na 2
                else
                    n = (np-nva-effet-indic_ALPHA-nb1 - netadc - netar)/(effet+1)
                endif    
                do k=1,n
                    if(typeJoint.eq.2) the2(k-3)=(bh(k))**2.d0
                    if(typeJoint.eq.3) the2(k-3)=(bh(k+n))**2.d0
                end do
                
                call susps(tps,the2,nz1+2,su,bbb,zi)    
               
            case(2) ! calcul du risque weibull
                if(typeJoint.eq.2) then
                    betaD = bh(1)**2
                    etaD = bh(2)**2
                else
                    betaD = bh(3)**2
                    etaD = bh(4)**2
                end if    
                if (tps.eq.0.d0) tps = 1d-12 ! utile car log(0) => -Inf    
                bbb = (betaD*dexp((betaD-1.d0)*dlog(tps))/(etaD**betaD)) !((tps/etaD)**betaD)!    
        end select    
    
        survdcCM_pred =bbb*dexp(XbetapredDCi)*dexp(etaydc(1)*current_m(1))!+cdc(i)*etaydc1*current_mean(1)
    
            return
    
        end function survdcCM_pred