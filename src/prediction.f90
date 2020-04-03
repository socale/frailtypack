

! ============================================== prediction Joint

    subroutine predict(np,b,nz,nbintervR,nbintervDC,nva1,nva2,nst,typeof0,typevent,zi,HIHOut,time,timedc, &
    ntimeAll,npred0,predTime,window,predtimerec,nrec0,vaxpred0,vaxdcpred0, &
    predAll1,predAll2,predAll3,predAll1R,predAlllow1,predAllhigh1, &
    predAlllow2,predAllhigh2,predAlllow3,predAllhigh3,predAlllow1R,predAllhigh1R,&
    icproba,nsample,intcens,trunctime,lowertime,uppertime,movingwindow,timeAll,modeltype,indic_alpha)

    implicit none

    integer::i,ii,iii,j,npred0,nrec0,nsample,typeof
    integer,intent(in)::np,nz,nbintervR,nbintervDC,nva1,nva2,nst,&
            typeof0, typevent,ntimeAll,icproba,intcens,movingwindow,modeltype,indic_alpha
    double precision,dimension(np),intent(in)::b
    double precision,dimension(nz+6),intent(in)::zi
    double precision,dimension(np,np),intent(in)::HIHOut
    double precision,dimension(nbintervR+1),intent(in)::time
    double precision,dimension(nbintervDC+1),intent(in)::timedc
    double precision,dimension(npred0,nva1),intent(in)::vaxpred0
    double precision,dimension(npred0,nva2),intent(in)::vaxdcpred0
    double precision,dimension(1,npred0)::XbetapredR,XbetapredDC,XbetapredRalea,XbetapredDCalea
    integer,dimension(npred0)::nreci
    double precision::predTime,window,predTime2,scR,shR,scDC,shDC, &
    scRalea,shRalea,scDCalea,shDCalea,alea,thesigalea,alphaalea,alpha,thesig
    double precision::ss11,ss12,ss21,ss22,ss31,ss32,ss11R,ss12R
    double precision,dimension(npred0)::predProba1,predProba2,predProba3,predProbaR1
    double precision,dimension(npred0,ntimeAll),intent(out)::predAll1,predAll2,predAll3,predAll1R
    double precision,dimension(npred0,ntimeAll),intent(out)::predAlllow1,predAllhigh1,predAlllow2, &
                predAllhigh2,predAlllow3,predAllhigh3, &
                predAlllow1R,predAllhigh1R
    double precision,dimension(nz+2)::theR,theDC,theRalea,theDCalea
    double precision,dimension(2)::surv,survDC,survDCalea,lam
    double precision,dimension(npred0)::survLT,survLTalea,survU,survL,survUalea,survLalea
    double precision,dimension(npred0,nrec0+2)::survR,hazR,survRalea,hazRalea
    double precision,dimension(nrec0+2)::survRi,hazRi,survRialea,hazRialea
    double precision,dimension(ntimeAll)::timeAll
    double precision,dimension(nsample,np)::balea
    double precision,dimension(nsample,npred0)::predProbaalea1,predProbaalea2,&
                                                predProbaalea3,predProbaalea1R
    double precision,dimension(1,nva1)::coefBetaalea
    double precision,dimension(1,nva2)::coefBetadcalea
    double precision,dimension(npred0,nrec0)::predtimerec
    double precision,dimension(npred0,nrec0+2)::predtimerec2
    double precision,dimension(npred0)::trunctime,lowertime,uppertime,lowertime2,uppertime2
    double precision,dimension(1,nva1)::coefBeta
    double precision,dimension(1,nva2)::coefBetadc
    
    !Myriam
    ! thesig = theta when joint gamma modelling; = sigma when joint log-normal modelling
    ! typevent = 2 if prediction for terminal event ; =3 if prediction for recurrent event ; =1 if prediction for the both
    ! modeltype =0 if joint gamma modelling; =1 if joint log-normal modelling

    typeof = typeof0
    coefBeta(1,:) = b((np-nva1-nva2+1):(np-nva2))
    coefBetadc(1,:) = b((np-nva2+1):np)

    XbetapredR = matmul(coefBeta,transpose(vaxpred0))
    XbetapredDC = matmul(coefBetadc,transpose(vaxdcpred0))

    ! Determine when to calculate Survival function (gap time)
    ! and the number of recurrences in the prediction (for each pred i)
    
    

!     timeAll(1) = predTime + window
!     do i=2,ntimeAll
!         timeAll(i) = timeAll(i-1) + window
!     end do

    if (icproba.eq.1) then ! generation des parametres
        do j=1,nsample
            do i=1,np
                call rnorm(b(i),sqrt(HIHOut(i,i)),alea)
                balea(j,i) = alea
            end do
        end do
    end if
    
    do iii=1,ntimeAll
        do i=1,npred0
            if (movingwindow.eq.1) then 
                predtimerec2(i,1) = predTime
            else 
                predtimerec2(i,1) = timeAll(iii) - window
            endif

            !predtimerec2(i,2) = predtimerec(i,1)
            nreci(i) = 0
            do ii=2,nrec0+1
                if (predtimerec(i,ii-1).le.predtimerec2(i,1)) then ! check if relapse happened before prediction time
                    predtimerec2(i,ii) = predtimerec(i,ii-1)
                else ! otherwise 0
                    predtimerec2(i,ii) = 0
                endif
                if ((ii.gt.2).and.(predtimerec2(i,ii-1).gt.0).and.(predtimerec2(i,ii).eq.0)) then
                    nreci(i) = ii-2
                endif
            end do
            if (predtimerec2(i,nrec0+1).gt.0) nreci(i) = nrec0
            predtimerec2(i,nrec0+2) = timeAll(iii)

            if (intcens.eq.1) then
                if (uppertime(i).gt.predtimerec2(i,1)) then
                    uppertime2(i) = 0.d0
                    lowertime2(i) = 0.d0
                else
                    uppertime2(i) = uppertime(i)
                    lowertime2(i) = lowertime(i)
                endif
            endif
            
        end do
        
        ! Calcul des risques de base
        ! A chaque fois, calculé pour : 
        ! DC au temps de base (predtimerec2(1,1)) et à l'horizon (predtimerec2(1,nrec0+2))
        ! Recurrence au temps de base et pour chaque temps de rechute entré (predtimerec2(i,ii))
        ! pour chaque prediction demandée

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
                if (intcens.eq.1) then 
                    call survival_frailty(trunctime(i),theR,theDC,nz+2,zi,surv,lam,nst)
                    survLT(i) = surv(1)
                    if (trunctime(i).eq.0.d0) survLT(i) = 1.d0
                    call survival_frailty(lowertime2(i),theR,theDC,nz+2,zi,surv,lam,nst)
                    survL(i) = surv(1)
                    call survival_frailty(uppertime2(i),theR,theDC,nz+2,zi,surv,lam,nst) !!
                    survU(i) = surv(1)
                endif

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

            case(1)

            predTime2 = predtimerec2(1,1)
            call survivalj_cpm(predTime2,b(1:(nbintervR+nbintervDC)),nbintervR,nbintervDC &
            ,time,timedc,surv)
            survR(:,1) = surv(1)
            survDC(1) = surv(2)
            do i=1,npred0
                if (intcens.eq.1) then
                    call survivalj_cpm(trunctime(i),b(1:(nbintervR+nbintervDC)),nbintervR,nbintervDC & !!
                    ,time,timedc,surv)
                    survLT(i) = surv(1)
                    if (trunctime(i).eq.0.d0) survLT(i) = 1.d0
                    call survivalj_cpm(lowertime2(i),b(1:(nbintervR+nbintervDC)),nbintervR,nbintervDC & !!
                    ,time,timedc,surv)
                    survL(i) = surv(1)
                    call survivalj_cpm(uppertime2(i),b(1:(nbintervR+nbintervDC)),nbintervR,nbintervDC & !!
                    ,time,timedc,surv)
                    survU(i) = surv(1)
                endif

                do ii=1,nrec0
                    predTime2 = predtimerec2(i,ii+1)
                    call survivalj_cpm(predTime2,b(1:(nbintervR+nbintervDC)),&
                    nbintervR,nbintervDC,time,timedc,surv)
                    survR(i,ii+1) = surv(1)
                end do
            end do
            predTime2 = predtimerec2(1,nrec0+2)
            call survivalj_cpm(predTime2,b(1:(nbintervR+nbintervDC)),nbintervR,nbintervDC &
            ,time,timedc,surv)
            survR(:,nrec0+2) = surv(1)
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
                if (intcens.eq.1) then
                    survLT(i) = exp(-(trunctime(i)/scR)**shR) !!
                    if (trunctime(i).eq.0.d0) survLT(i) = 1.d0
                    survL(i) = exp(-(lowertime2(i)/scR)**shR) !!
                    survU(i) = exp(-(uppertime2(i)/scR)**shR) !!
                endif

                do ii=1,nrec0+1
                    survR(i,ii+1) = exp(-(predtimerec2(i,ii+1)/scR)**shR)
                    hazR(i,ii+1) = (shR/scR)*((predtimerec2(i,ii+1)/scR)**(shR-1))
                end do
            end do
            survDC(2) = exp(-(predtimerec2(1,nrec0+2)/scDC)**shDC)
        end select

        thesig = b(np-nva1-nva2-1)*b(np-nva1-nva2-1)
        if(indic_alpha.eq.1) then 
            alpha = b(np-nva1-nva2)
        else 
            alpha = 1.d0
        end if

        do i=1,npred0
            if (intcens.eq.0) then
                survRi = survR(i,:)
                hazRi = hazR(i,:)
                if (modeltype == 0) then
                    call gaulagJpred(ss11,ss12,ss21,ss22,ss31,ss32,ss11R,ss12R,&
                                 thesig,alpha,typevent,XbetapredR(1,i),XbetapredDC(1,i),&
                                 survRi,hazRi, survDC,nrec0,nreci(i))
                else                     
                    call gauherJpred(ss11,ss12,ss21,ss22,ss31,ss32,ss11R,ss12R,&
                                 thesig,alpha,typevent,XbetapredR(1,i),XbetapredDC(1,i),&
                                 survRi,hazRi, survDC,nrec0,nreci(i))
                endif

                if ((typevent.eq.1).or.(typevent.eq.2)) then
                    predProba1(i) = ss11/ss12
                    predProba2(i) = ss21/ss22
                    predProba3(i) = ss31/ss32
                endif
                if ((typevent.eq.1).or.(typevent.eq.3)) then
                    predProbaR1(i) = ss11R/ss12R
                endif
                
            else
                if (modeltype == 0) then
                    call gaulagJpredic(ss21,ss22,thesig,alpha,XbetapredR(1,i),XbetapredDC(1,i),hazRi, &
                        survDC,nrec0,survL(i),survU(i),survLT(i))
                else 
                    call gauherJpredic(ss21,ss22,thesig,alpha,XbetapredR(1,i),XbetapredDC(1,i),hazRi, &
                        survDC,nrec0,survL(i),survU(i),survLT(i))
                endif
                
                predProba1(i) = 0.d0
                predProba2(i) = ss21/ss22
                predProba3(i) = 0.d0
            endif
               
        end do

        if ((typevent.eq.1).or.(typevent.eq.2)) then
            predAll1(:,iii) = predProba1
            predAll2(:,iii) = predProba2
            predAll3(:,iii) = predProba3
        endif
        if ((typevent.eq.1).or.(typevent.eq.3)) then
            predAll1R(:,iii) = predProbaR1            
        endif

        !=============================================
        ! Variabilite des proba predites
        ! Creation d'un vecteur balea, qui correspond au vecteur b où chaque parametre
        ! est tiré au sort selon sa loi
!        seProba1(:)=0.d0; seProba2(:)=0.d0; seProba3(:)=0.d0;seProba4(:)=0.d0;
!        lowProba1(:)=0.d0; lowProba2(:)=0.d0; lowProba3(:)=0.d0;lowProba4(:)=0.d0;
!        highProba1(:)=0.d0; highProba2(:)=0.d0; highProba3(:)=0.d0;highProba4(:)=0.d0;
!        predProbaalea1(:,:)=0.d0;predProbaalea2(:,:)=0.d0;
!        predProbaalea3(:,:)=0.d0;predProbaalea4(:,:)=0.d0;

        if (icproba.eq.1) then ! calcul de l'intervalle de confiance seulement si demande

            do j=1,nsample
                ss11 = 0.d0
                ss12 = 0.d0
                ss21 = 0.d0
                ss22 = 0.d0
                ss31 = 0.d0
                ss32 = 0.d0
                ss11R = 0.d0
                ss12R = 0.d0
                XbetapredRalea = 0.d0
                XbetapredDCalea = 0.d0
                survRalea = 0.d0
                hazRalea = 0.d0
                survDCalea = 0.d0
                survLTalea = 0.d0

                coefBetaalea(1,:) = balea(j,(np-nva1-nva2+1):(np-nva2))
                coefBetadcalea(1,:) = balea(j,(np-nva2+1):np)

                XbetapredRalea = matmul(coefBetaalea,transpose(vaxpred0))
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
                        if (intcens.eq.1) then 
                            call survival_frailty(trunctime(i),theRalea,theDCalea,nz+2,zi,surv,lam,nst)
                            survLTalea(i) = surv(1)
                            if (trunctime(i).eq.0.d0) survLTalea(i) = 1.d0
                            call survival_frailty(lowertime(i),theRalea,theDCalea,nz+2,zi,surv,lam,nst)
                            survLalea(i) = surv(1)
                            call survival_frailty(uppertime(i),theRalea,theDCalea,nz+2,zi,surv,lam,nst) !!
                            survUalea(i) = surv(1)
                        endif

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

                    case(1)

                    predTime2 = predtimerec2(1,1)
                    call survivalj_cpm(predTime2,balea(j,1:(nbintervR+nbintervDC)),nbintervR,nbintervDC &
                    ,time,timedc,surv)
                    survRalea(:,1) = surv(1)
                    survDCalea(1) = surv(2)
                    do i=1,npred0
                        if (intcens.eq.1) then
                            call survivalj_cpm(trunctime(i),balea(j,1:(nbintervR+nbintervDC)),nbintervR,nbintervDC & !!
                            ,time,timedc,surv)
                            survLTalea(i) = surv(1)
                            if (trunctime(i).eq.0.d0) survLTalea(i) = 1.d0
                            call survivalj_cpm(lowertime(i),balea(j,1:(nbintervR+nbintervDC)),nbintervR,nbintervDC & !!
                            ,time,timedc,surv)
                            survLalea(i) = surv(1)
                            call survivalj_cpm(uppertime(i),balea(j,1:(nbintervR+nbintervDC)),nbintervR,nbintervDC & !!
                            ,time,timedc,surv)
                            survUalea(i) = surv(1)
                        endif

                        do ii=1,nrec0
                            predTime2 = predtimerec2(i,ii+1)
                            call survivalj_cpm(predTime2,balea(j,1:(nbintervR+nbintervDC)),&
                            nbintervR,nbintervDC,time,timedc,surv)
                            survRalea(i,ii+1) = surv(1)
                        end do
                    end do
                    predTime2 = predtimerec2(1,nrec0+2)
                    call survivalj_cpm(predTime2,balea(j,1:(nbintervR+nbintervDC)),nbintervR,nbintervDC &
                    ,time,timedc,surv)
                    survRalea(:,nrec0+2) = surv(1)
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
                        if (intcens.eq.1) then
                            survLTalea(i) = exp(-(trunctime(i)/scRalea)**shRalea) !!
                            if (trunctime(i).eq.0.d0) survLTalea(i) = 1.d0
                            survLalea(i) = exp(-(lowertime(i)/scRalea)**shRalea) !!
                            survUalea(i) = exp(-(uppertime(i)/scRalea)**shRalea) !!
                        endif

                        do ii=1,nrec0+1
                            survRalea(i,ii+1) = exp(-(predtimerec2(i,ii+1)/scRalea)**shRalea)
                            hazRalea(i,ii+1) = (shRalea/scRalea)*((predtimerec2(i,ii+1)/scRalea)**(shRalea-1))
                        end do
                    end do
                    survDCalea(2) = exp(-(predtimerec2(1,nrec0+2)/scDCalea)**shDCalea)
                end select

                if(indic_alpha.eq.1) then 
                    alphaalea = balea(j,np-nva1-nva2)
                else 
                    alphaalea = 1.d0
                endif
                thesigalea = balea(j,np-nva1-nva2-1)*balea(j,np-nva1-nva2-1)
                do i=1,npred0
                    if (intcens.eq.0) then
                        survRialea = survRalea(i,:)
                        hazRialea = hazRalea(i,:)
                        if (modeltype == 0) then
                            call gaulagJpred(ss11,ss12,ss21,ss22,ss31,ss32, &
                                        ss11R,ss12R,thesigalea,alphaalea,typevent,XbetapredRalea(1,i), &
                                        XbetapredDCalea(1,i),survRialea,hazRialea,&
                                        survDCalea,nrec0,nreci(i))
                        else 
                            call gauherJpred(ss11,ss12,ss21,ss22,ss31,ss32, &
                                        ss11R,ss12R,thesigalea,alphaalea,typevent,XbetapredRalea(1,i), &
                                        XbetapredDCalea(1,i),survRialea,hazRialea,&
                                        survDCalea,nrec0,nreci(i))
                        endif
                        
                        if ((typevent.eq.1).or.(typevent.eq.2)) then
                            predProbaalea1(j,i) = ss11/ss12
                            predProbaalea2(j,i) = ss21/ss22
                            predProbaalea3(j,i) = ss31/ss32
                        endif
                        if ((typevent.eq.1).or.(typevent.eq.3)) then
                            predProbaalea1R(j,i) = ss11R/ss12R
                        endif
                    else
                        if (modeltype == 0) then
                            call gaulagJpredic(ss21,ss22,thesigalea,alphaalea,XbetapredRalea(1,i),XbetapredDCalea(1,i), &
                                        hazRialea,survDCalea,nrec0,survLalea(i),survUalea(i),survLTalea(i))
                        else 
                            call gauherJpredic(ss21,ss22,thesigalea,alphaalea,XbetapredRalea(1,i),XbetapredDCalea(1,i), &
                                        hazRialea,survDCalea,nrec0,survLalea(i),survUalea(i),survLTalea(i))
                        endif
                        predProbaalea1(j,i) = 0.d0
                        predProbaalea2(j,i) = ss21/ss22
                        predProbaalea3(j,i) = 0.d0
                    endif
                end do
            end do

            ! utilisation de la fonction percentile2 de aaUseFunction
            if ((typevent.eq.1).or.(typevent.eq.2)) then
                do i=1,npred0
                    call percentile2(predProbaalea1(:,i),nsample,predAlllow1(i,iii),predAllhigh1(i,iii))
                    call percentile2(predProbaalea2(:,i),nsample,predAlllow2(i,iii),predAllhigh2(i,iii))
                    call percentile2(predProbaalea3(:,i),nsample,predAlllow3(i,iii),predAllhigh3(i,iii))
                end do
            endif
            if ((typevent.eq.1).or.(typevent.eq.3)) then
                do i=1,npred0
                    call percentile2(predProbaalea1R(:,i),nsample,predAlllow1R(i,iii),predAllhigh1R(i,iii))!Myriam : prediction P1 pour recurrence
                end do
            endif
            
        endif ! calcul de l'intervalle de confiance seulement si demande

    end do
    end subroutine predict

!=========================
! Prediction 1 : exactement j recurrences
!=========================
    double precision function func1pred1(frail,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
    ! calcul de l integrant (numerateur de la fonction de prediction)
    implicit none

    integer::nrec0,recj
    double precision,intent(in)::frail
    double precision::XbetapredRi,XbetapredDCi
    double precision,dimension(2)::survDC
    double precision,dimension(nrec0+2)::survRi
    double precision::ptheta,palpha,logGammaJ

    func1pred1 = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi)) &
    - survDC(2)**((frail**palpha) * exp(XbetapredDCi))) &
    * (frail**recj) &
    * (survRi(1)**(frail * exp(XbetapredRi))) &
    * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta)))) 
    
      
    return

    end function func1pred1

    double precision function func2pred1(frail,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
    ! calcul de l integrant (denominateur de la fonction de prediction)
    implicit none

    integer::nrec0,recj
    double precision,intent(in)::frail
    double precision::XbetapredRi,XbetapredDCi
    double precision,dimension(2)::survDC
    double precision,dimension(nrec0+2)::survRi
    double precision::ptheta,palpha,logGammaJ

    func2pred1 = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi))) &
    * (frail**recj) &
    * (survRi(1)**(frail * exp(XbetapredRi))) &
    * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / (ptheta**(1.d0/ptheta) *dexp( logGammaJ(1.d0/ptheta))) )
    return

    end function func2pred1

!=========================
! Prediction 2 : au moins j recurrences
!=========================
    double precision function func1pred2(frail,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)

    implicit none

    integer::nrec0,recj
    double precision,intent(in)::frail
    double precision::XbetapredRi,XbetapredDCi
    double precision,dimension(2)::survDC
    double precision,dimension(nrec0+2)::survRi
    double precision::ptheta,palpha,logGammaJ

    func1pred2 = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi)) &
    - survDC(2)**((frail**palpha) * exp(XbetapredDCi))) &
    * (frail**recj) &
    * ((survRi(recj+1))**(frail * exp(XbetapredRi))) &
    * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta)))) 

    return

    end function func1pred2

    double precision function func2pred2(frail,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)

    implicit none

    integer::nrec0,recj
    double precision,intent(in)::frail
    double precision::XbetapredRi,XbetapredDCi
    double precision,dimension(2)::survDC
    double precision,dimension(nrec0+2)::survRi
    double precision::ptheta,palpha,logGammaJ

    func2pred2 = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi))) &
    * (frail**recj) &
    * ((survRi(recj+1))**(frail * exp(XbetapredRi))) & !!
    * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta))))

    return

    end function func2pred2

!=========================
! Prediction 2 : au moins j recurrences pour la censure par intervalle !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=========================
    double precision function func1pred2ic(frail,ptheta,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!

    implicit none

    !integer::nrec0
    double precision,intent(in)::frail
    double precision::XbetapredRi,XbetapredDCi
    double precision,dimension(2)::survDC
    double precision::survLT,survL,survU !!
    double precision::ptheta,palpha,logGammaJ

    if ((survL.eq.1.d0).or.(survU.eq.1.d0)) then
        func1pred2ic = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi)) &
        - survDC(2)**((frail**palpha) * exp(XbetapredDCi))) &
        !* (frail**recj) &
        !* ((survL**(frail * exp(XbetapredRi))-survU**(frail * exp(XbetapredRi)))
        /(survLT**(frail * exp(XbetapredRi))) & !!
        * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta))))
    else
        func1pred2ic = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi)) &
        - survDC(2)**((frail**palpha) * exp(XbetapredDCi))) &
        !* (frail**recj) &
        * ((survL**(frail * exp(XbetapredRi))-survU**(frail * exp(XbetapredRi)))/survLT**(frail * exp(XbetapredRi))) & !!
        * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta))))
    endif

    return

    end function func1pred2ic

    double precision function func2pred2ic(frail,ptheta,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!

    implicit none

    !integer::nrec0
    double precision,intent(in)::frail
    double precision::XbetapredRi,XbetapredDCi
    double precision,dimension(2)::survDC
    double precision::survLT,survL,survU!!
    double precision::ptheta,palpha,logGammaJ

    if ((survL.eq.1.d0).or.(survU.eq.1.d0)) then
        func2pred2ic = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi))) &
        !* (frail**recj) &
        !* ((survL**(frail * exp(XbetapredRi))-survU**(frail * exp(XbetapredRi)))
        /(survLT**(frail * exp(XbetapredRi))) & !!
        * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta))))
    else
        func2pred2ic = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi))) &
        !* (frail**recj) &
        * ((survL**(frail * exp(XbetapredRi))-survU**(frail * exp(XbetapredRi)))/survLT**(frail * exp(XbetapredRi))) & !!
        * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta))))
    endif

    return

    end function func2pred2ic

!=========================
! Prediction 3 : ignorant les recurrences precedentes
!=========================
    double precision function func1pred3(frail,ptheta,palpha,XbetapredDCi,survDC)

    implicit none

    double precision,intent(in)::frail
    double precision::XbetapredDCi
    double precision,dimension(2)::survDC
    double precision::ptheta,palpha,logGammaJ

    func1pred3 = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi)) &
    - survDC(2)**((frail**palpha) * exp(XbetapredDCi))) &
    * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta))) )

    return

    end function func1pred3

    double precision function func2pred3(frail,ptheta,palpha,XbetapredDCi,survDC)

    implicit none

    double precision,intent(in)::frail
    double precision::XbetapredDCi
    double precision,dimension(2)::survDC
    double precision::ptheta,palpha,logGammaJ

    func2pred3 = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi))) &
    * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta))))

    return

    end function func2pred3

! Myriam - 07/2016 - !!!
!=========================
! Prediction 1 pour les evenements recurrents, tenant compte de exactement j recurrences
!=========================
      double precision function func1pred1_rec(frail,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
    ! calcul de l integrant (numerateur de la fonction de prediction)
      implicit none

      integer::nrec0,recj
      double precision,intent(in)::frail
      double precision::XbetapredRi,XbetapredDCi
      double precision,dimension(2)::survDC
      double precision,dimension(nrec0+2)::survRi
      double precision::ptheta,palpha,logGammaJ

      func1pred1_rec = ((survRi(1)**(frail* exp(XbetapredRi)) &
      - survRi(nrec0+2)**(frail* exp(XbetapredRi))) &
      * survDC(1)**((frail**palpha) * exp(XbetapredDCi)) &
      * (frail**recj) &
      * (survRi(recj+1)**(frail * exp(XbetapredRi))) &
      * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta)) / &
      (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta))))
      return

      end function func1pred1_rec

      double precision function func2pred1_rec(frail,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
    ! calcul de l integrant (denominateur de la fonction de prediction)
      implicit none

      integer::nrec0,recj
      double precision,intent(in)::frail
      double precision::XbetapredRi,XbetapredDCi
      double precision,dimension(2)::survDC
      double precision,dimension(nrec0+2)::survRi
      double precision::ptheta,palpha,logGammaJ

      func2pred1_rec = (survRi(1)**(frail* exp(XbetapredRi)) &
      * (survDC(1)**((frail**palpha) * exp(XbetapredDCi))) &
      * (frail**recj) &
      * (survRi(recj+1)**(frail * exp(XbetapredRi))) &
      * (frail**(1.d0/ptheta -1.d0) * exp(-frail/ptheta))/ &
      (ptheta**(1.d0/ptheta) * dexp(logGammaJ(1.d0/ptheta)))) 

      return
    
      end function func2pred1_rec

!=========================
! Calcul des intégrales pour un modèle gamma
!=========================
    subroutine gaulagJpred(ss11,ss12,ss21,ss22,ss31,ss32,ss11R,ss12R,&
        ptheta,palpha,typevent,XbetapredRi,XbetapredDCi,survRi,hazRi,survDC,nrec0,recj) !! 

!    use tailles
    use comon,only:typeof
    use donnees,only:w,x,w1,x1
    
    implicit none

    integer,intent(in)::nrec0,recj
    double precision,intent(out)::ss11,ss12,ss21,ss22,ss31,ss32,ss11R,ss12R
    double precision::XbetapredRi,XbetapredDCi,var1
    double precision::auxfunca11,auxfunca12,auxfunca21,auxfunca22, &
                      auxfunca31,auxfunca32,auxfunca11_rec,auxfunca12_rec
    double precision,external :: func1pred1,func2pred1,func1pred2, & 
                                 func2pred2,func1pred3,func2pred3, &
                                 func2pred1_rec, func1pred1_rec
    double precision,dimension(2)::survDC
    !double precision::survLT !!
    double precision,dimension(nrec0)::survRi,hazRi,hazRi2
    double precision::ptheta,palpha
    integer:: j,typevent

! gauss laguerre
! func1 est l integrant, ss le resultat de l integrale sur 0 ,  +infty

    hazRi2 = hazRi
    ss11=0.d0
    ss12=0.d0
    ss21=0.d0
    ss22=0.d0
    ss31=0.d0
    ss32=0.d0
    ss11R=0.d0
    ss12R=0.d0

    if (typeof == 0)then  
! Will be twice the average value of the function,since the ten
! weights (five numbers above each used twice) sum to 2.
        if ((typevent == 1).or.(typevent == 2)) then
            do j=1,20
                var1 = x(j)
                auxfunca11 = func1pred1(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                
                ss11 = ss11 + w(j)*(auxfunca11)
                auxfunca12 = func2pred1(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss12 = ss12 + w(j)*(auxfunca12)
                auxfunca21 = func1pred2(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss21 = ss21 + w(j)*(auxfunca21)
                auxfunca22 = func2pred2(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss22 = ss22 + w(j)*(auxfunca22)
                auxfunca31 = func1pred3(var1,ptheta,palpha,XbetapredDCi,survDC)
                ss31 = ss31 + w(j)*(auxfunca31)
                auxfunca32 = func2pred3(var1,ptheta,palpha,XbetapredDCi,survDC)
                ss32 = ss32 + w(j)*(auxfunca32)
            end do
        endif
        if ((typevent == 1).or.(typevent == 3)) then 
            do j=1,20
                var1 = x(j)
                auxfunca11_rec = func1pred1_rec(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss11R = ss11R + w(j)*(auxfunca11_rec)
                auxfunca12_rec = func2pred1_rec(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss12R = ss12R + w(j)*(auxfunca12_rec)
            end do
        endif
    else
        if ((typevent == 1).or.(typevent == 2)) then
            do j=1,32
                var1 = x1(j)
                auxfunca11 = func1pred1(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss11 = ss11 + w1(j)*(auxfunca11)
                auxfunca12 = func2pred1(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss12 = ss12 + w1(j)*(auxfunca12)
                auxfunca21 = func1pred2(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss21 = ss21 + w1(j)*(auxfunca21)
                auxfunca22 = func2pred2(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss22 = ss22 + w1(j)*(auxfunca22)
                auxfunca31 = func1pred3(var1,ptheta,palpha,XbetapredDCi,survDC)
                ss31 = ss31 + w1(j)*(auxfunca31)
                auxfunca32 = func2pred3(var1,ptheta,palpha,XbetapredDCi,survDC)
                ss32 = ss32 + w1(j)*(auxfunca32)
            end do
        endif
        if ((typevent == 1).or.(typevent == 3)) then  
            do j=1,32
                var1 = x1(j)
                !Myriam : prediction pour recidive
                auxfunca11_rec = func1pred1_rec(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss11R = ss11R + w1(j)*(auxfunca11_rec)
                auxfunca12_rec = func2pred1_rec(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                ss12R = ss12R + w1(j)*(auxfunca12_rec)
            end do
        endif
    end if
    end subroutine gaulagJpred

    subroutine gaulagJpredic(ss21,ss22,ptheta,palpha,XbetapredRi,XbetapredDCi,hazRi,survDC,nrec0,survL,survU,survLT) !! 

!    use tailles
    use comon,only:typeof
    use donnees,only:w,x,w1,x1

    implicit none

    integer,intent(in)::nrec0
    double precision,intent(out)::ss21,ss22
    double precision::XbetapredRi,XbetapredDCi,var1
    double precision::auxfunca21,auxfunca22
    double precision,external ::func1pred2ic,func2pred2ic
    double precision,dimension(2)::survDC
    double precision::survLT,survL,survU !!
    double precision,dimension(nrec0)::hazRi,hazRi2
    double precision::ptheta,palpha
    integer:: j

! gauss laguerre
! func1 est l integrant, ss le resultat de l integrale sur 0 ,  +infty

    hazRi2 = hazRi
    ss21=0.d0
    ss22=0.d0

    if (typeof == 0)then  
! Will be twice the average value of the function,since the ten
! weights (five numbers above each used twice) sum to 2.
        do j=1,20
            var1 = x(j)
            auxfunca21 = func1pred2ic(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!
            ss21 = ss21 + w(j)*(auxfunca21)
            auxfunca22 = func2pred2ic(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!
            ss22 = ss22 + w(j)*(auxfunca22)
        end do
    else
        do j=1,32
            var1 = x1(j)
            auxfunca21 = func1pred2ic(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!
            ss21 = ss21 + w1(j)*(auxfunca21)
            auxfunca22 = func2pred2ic(var1,ptheta,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!
            ss22 = ss22 + w1(j)*(auxfunca22)
        end do
    end if
    end subroutine gaulagJpredic
    
!=============================================================================================================
!=============================================================================================================
!========================================Modele log normal====================================================    
    
!=========================
! Prediction 1 : exactement j recurrences (modele log normal)
!=========================
      double precision function func1pred1LogN(frail,psig2,palpha,&
                XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
    ! calcul de l integrant (numerateur de la fonction de prediction)
      implicit none

      integer::nrec0,recj
      double precision,intent(in)::frail
      double precision::XbetapredRi,XbetapredDCi
      double precision,dimension(2)::survDC
      double precision,dimension(nrec0+2)::survRi
      double precision::psig2,palpha
      double precision,parameter::pi=3.141592653589793d0

      func1pred1LogN = ((survDC(1)**((dexp(frail*palpha))& 
      * exp(XbetapredDCi)) - survDC(2)**((dexp(frail*palpha))& 
      * exp(XbetapredDCi))) * (dexp(frail)**recj) &
      * (survRi(1)**(dexp(frail) * exp(XbetapredRi))) &
      * dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))

      return

      end function func1pred1LogN

      double precision function func2pred1LogN(frail,psig2,palpha,&
                 XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
    ! calcul de l integrant (denominateur de la fonction de prediction)
      implicit none

      integer::nrec0,recj
      double precision,intent(in)::frail
      double precision::XbetapredRi,XbetapredDCi
      double precision,dimension(2)::survDC
      double precision,dimension(nrec0+2)::survRi
      double precision::psig2,palpha
      double precision,parameter::pi=3.141592653589793d0

      func2pred1LogN = ((survDC(1)**((dexp(frail*palpha)) &
      * exp(XbetapredDCi))) * (dexp(frail)**recj) &
      * (survRi(1)**(dexp(frail) * exp(XbetapredRi))) &
      *dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))

      return

      end function func2pred1LogN

!=========================
! Prediction 2 : au moins j recurrences (modele log normal)
!=========================
    double precision function func1pred2LogN(frail,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)

    implicit none

    integer::nrec0,recj
    double precision,intent(in)::frail
    double precision::XbetapredRi,XbetapredDCi
    double precision,dimension(2)::survDC
    double precision,dimension(nrec0+2)::survRi
    double precision::psig2,palpha
    double precision,parameter::pi=3.141592653589793d0

    func1pred2LogN = ((survDC(1)**((dexp(frail*palpha)) * exp(XbetapredDCi)) &
    - survDC(2)**((dexp(frail*palpha)) * exp(XbetapredDCi))) &
    * (dexp(frail)**recj) &
    * ((survRi(recj+1))**(dexp(frail) * exp(XbetapredRi))) &
    *dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))

    return

    end function func1pred2LogN

    double precision function func2pred2LogN(frail,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)

    implicit none

    integer::nrec0,recj
    double precision,intent(in)::frail
    double precision::XbetapredRi,XbetapredDCi
    double precision,dimension(2)::survDC
    double precision,dimension(nrec0+2)::survRi
    double precision::psig2,palpha
    double precision,parameter::pi=3.141592653589793d0

    func2pred2LogN = ((survDC(1)**((dexp(frail*palpha)) * exp(XbetapredDCi))) &
    * (dexp(frail)**recj) &
    * ((survRi(recj+1))**(dexp(frail) * exp(XbetapredRi))) & !!
     *dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))

    return
    
    end function func2pred2LogN

!=========================
! Prediction 2 : au moins j recurrences (modele log normal) pour la censure par intervalle !!!!!!!!
!=========================
    double precision function func1pred2LogNicLogN(frail,psig2,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!

    implicit none

    !integer::nrec0
    double precision,intent(in)::frail
    double precision::XbetapredRi,XbetapredDCi
    double precision,dimension(2)::survDC
    double precision::survLT,survL,survU !!
    double precision::psig2,palpha
    double precision,parameter::pi=3.141592653589793d0

    if ((survL.eq.1.d0).or.(survU.eq.1.d0)) then
        func1pred2LogNicLogN = ((survDC(1)**((dexp(frail*palpha)) * exp(XbetapredDCi)) &
        - survDC(2)**((dexp(frail*palpha)) * exp(XbetapredDCi))) &
        !* (frail**recj) &
        !* ((survL**(frail * exp(XbetapredRi))-survU**(frail * exp(XbetapredRi)))
        /(survLT**(dexp(frail) * exp(XbetapredRi))) & !!
        *dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))
    else
        func1pred2LogNicLogN = ((survDC(1)**((frail**palpha) * exp(XbetapredDCi)) &
        - survDC(2)**((frail**palpha) * exp(XbetapredDCi))) &
        !* (frail**recj) &
        * ((survL**(dexp(frail) * exp(XbetapredRi))-survU**(dexp(frail) &
        * exp(XbetapredRi)))/survLT**(dexp(frail) * exp(XbetapredRi))) & !!
         *dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))
    endif

    return

    end function func1pred2LogNicLogN

    double precision function func2pred2LogNicLogN(frail,psig2,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!

    implicit none

    !integer::nrec0
    double precision,intent(in)::frail
    double precision::XbetapredRi,XbetapredDCi
    double precision,dimension(2)::survDC
    double precision::survLT,survL,survU!!
    double precision::psig2,palpha
    double precision,parameter::pi=3.141592653589793d0

    if ((survL.eq.1.d0).or.(survU.eq.1.d0)) then
        func2pred2LogNicLogN = ((survDC(1)**((dexp(frail*palpha)) * exp(XbetapredDCi))) &
        !* (frail**recj) &
        !* ((survL**(frail * exp(XbetapredRi))-survU**(frail * exp(XbetapredRi)))
        /(survLT**(dexp(frail) * exp(XbetapredRi))) & !!
        *dexp(-(dexp(frail)**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))
    else
        func2pred2LogNicLogN = ((survDC(1)**((dexp(frail*palpha)) * exp(XbetapredDCi))) &
        !* (frail**recj) &
        * ((survL**(dexp(frail) * exp(XbetapredRi))-survU**(dexp(frail)&
        * exp(XbetapredRi)))/survLT**(dexp(frail) * exp(XbetapredRi))) & !!
        *dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))
    endif
    
    return
    
    end function func2pred2LogNicLogN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=========================
! Prediction 3 : ignorant les recurrences precedentes (modele log normal)
!=========================
    double precision function func1pred3LogN(frail,psig2,palpha,XbetapredDCi,survDC)
    
    implicit none
    
    double precision,intent(in)::frail
    double precision::XbetapredDCi
    double precision,dimension(2)::survDC
    double precision::psig2,palpha
    double precision,parameter::pi=3.141592653589793d0
    
    func1pred3LogN = ((survDC(1)**((dexp(frail*palpha)) * exp(XbetapredDCi)) &
    - survDC(2)**((dexp(frail*palpha)) * exp(XbetapredDCi))) &
    *dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))
    
    return
    
    end function func1pred3LogN

    double precision function func2pred3LogN(frail,psig2,palpha,XbetapredDCi,survDC)
    
    implicit none
    
    double precision,intent(in)::frail
    double precision::XbetapredDCi
    double precision,dimension(2)::survDC
    double precision::psig2,palpha
    double precision,parameter::pi=3.141592653589793d0
    
    func2pred3LogN = ((survDC(1)**((dexp(frail*palpha)) * exp(XbetapredDCi))) &
   *dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))
    
    return
    
    end function func2pred3LogN

! Myriam - 08/2016 - !!!
!=========================
! Prediction 1 pour les evenements recurrents, tenant compte de exactement j recurrences (modele log normal)
!=========================
      double precision function func1pred1LogN_rec(frail,psig2,palpha,&
        XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
    ! calcul de l integrant (numerateur de la fonction de prediction)
      implicit none
    
      integer::nrec0,recj
      double precision,intent(in)::frail
      double precision::XbetapredRi,XbetapredDCi
      double precision,dimension(2)::survDC
      double precision,dimension(nrec0+2)::survRi
      double precision::psig2,palpha
      double precision,parameter::pi=3.141592653589793d0

      func1pred1LogN_rec = ((survRi(1)**(dexp(frail)* exp(XbetapredRi)) &
      - survRi(nrec0+2)**(dexp(frail)* exp(XbetapredRi))) &
      * survDC(1)**((dexp(frail*palpha)) * exp(XbetapredDCi)) &
      * (dexp(frail)**recj) &
      * (survRi(recj+1)**(dexp(frail) * exp(XbetapredRi))) &
      * dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))

      return
    
      end function func1pred1LogN_rec

      double precision function func2pred1LogN_rec(frail,psig2,palpha,&
            XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
    ! calcul de l integrant (denominateur de la fonction de prediction)
      implicit none

      integer::nrec0,recj
      double precision,intent(in)::frail
      double precision::XbetapredRi,XbetapredDCi
      double precision,dimension(2)::survDC
      double precision,dimension(nrec0+2)::survRi
      double precision::psig2,palpha
      double precision,parameter::pi=3.141592653589793d0

      func2pred1LogN_rec = (survRi(1)**(dexp(frail)* exp(XbetapredRi)) &
      * (survDC(1)**((dexp(frail*palpha)) * exp(XbetapredDCi))) &
      * (dexp(frail)**recj) &
      * (survRi(recj+1)**(dexp(frail) * exp(XbetapredRi)))&
      * dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2)))

      return

      end function func2pred1LogN_rec
    
!=================================================
! Calcul des intégrales pour un modèle log-normal
!=================================================    
    subroutine gauherJpred(ss11,ss12,ss21,ss22,ss31,ss32,ss11R,ss12R,&
        psig2,palpha,typevent,XbetapredRi,XbetapredDCi,survRi,hazRi,survDC,nrec0,recj) !!
    
!    use tailles
        use comon,only:typeof
        use donnees,only:w2,x2,w3,x3
        
        implicit none
        
        integer,intent(in)::nrec0,recj
        double precision,intent(out)::ss11,ss12,ss21,ss22,ss31,ss32,ss11R,ss12R
        double precision::XbetapredRi,XbetapredDCi,var1
        double precision::auxfunca11,auxfunca12,auxfunca21,auxfunca22,&
              auxfunca31,auxfunca32,auxfunca11_rec,auxfunca12_rec
        double precision,external :: func1pred1LogN,func2pred1LogN,&
              func1pred2LogN,func2pred2LogN,func1pred3LogN,&
              func2pred3LogN,func1pred1LogN_rec,func2pred1LogN_rec
        double precision,dimension(2)::survDC
        !double precision::survLT !!
        double precision,dimension(nrec0)::survRi,hazRi,hazRi2
        double precision::psig2,palpha
        integer:: j, typevent

! gauss laguerre
! func1 est l integrant, ss le resultat de l integrale sur 0 ,  +infty

        hazRi2 = hazRi
        ss11=0.d0
        ss12=0.d0
        ss21=0.d0
        ss22=0.d0
        ss31=0.d0
        ss32=0.d0
        ss11R=0.d0
        ss12R=0.d0

        if (typeof == 0)then  
! Will be twice the average value of the function,since the ten
! weights (five numbers above each used twice) sum to 2.
            if ((typevent == 1).or.(typevent == 2)) then
                do j=1,20
                    var1 = x2(j)
                    auxfunca11 = func1pred1LogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss11 = ss11 + w2(j)*(auxfunca11)
                    auxfunca12 = func2pred1LogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss12 = ss12 + w2(j)*(auxfunca12)
                    auxfunca21 = func1pred2LogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss21 = ss21 + w2(j)*(auxfunca21)
                    auxfunca22 = func2pred2LogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss22 = ss22 + w2(j)*(auxfunca22)
                    auxfunca31 = func1pred3LogN(var1,psig2,palpha,XbetapredDCi,survDC)
                    ss31 = ss31 + w2(j)*(auxfunca31)
                    auxfunca32 = func2pred3LogN(var1,psig2,palpha,XbetapredDCi,survDC)
                    ss32 = ss32 + w2(j)*(auxfunca32)
                end do
            endif
            if ((typevent == 1).or.(typevent == 3)) then
                do j=1,20
                    var1 = x2(j)
                    !Myriam : prediction pour recidive
                    auxfunca11_rec = func1pred1LogN_rec(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss11R = ss11R + w2(j)*(auxfunca11_rec)
                    auxfunca12_rec = func2pred1LogN_rec(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss12R = ss12R + w2(j)*(auxfunca12_rec)
                end do
            endif
        else
            if ((typevent == 1).or.(typevent == 2)) then
                do j=1,32
                    var1 = x3(j)
                    auxfunca11 = func1pred1LogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss11 = ss11 + w3(j)*(auxfunca11)
                    auxfunca12 = func2pred1LogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss12 = ss12 + w3(j)*(auxfunca12)
                    auxfunca21 = func1pred2LogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss21 = ss21 + w3(j)*(auxfunca21)
                    auxfunca22 = func2pred2LogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss22 = ss22 + w3(j)*(auxfunca22)
                    auxfunca31 = func1pred3LogN(var1,psig2,palpha,XbetapredDCi,survDC)
                    ss31 = ss31 + w3(j)*(auxfunca31)
                    auxfunca32 = func2pred3LogN(var1,psig2,palpha,XbetapredDCi,survDC)
                    ss32 = ss32 + w3(j)*(auxfunca32)
                end do
            endif
            if ((typevent == 1).or.(typevent == 3)) then  
                do j=1,32
                    !Myriam : prediction pour recidive
                    auxfunca11_rec = func1pred1LogN_rec(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss11R = ss11R + w3(j)*(auxfunca11_rec)
                    auxfunca12_rec = func2pred1LogN_rec(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survRi,survDC,nrec0,recj)
                    ss12R = ss12R + w3(j)*(auxfunca12_rec)            
                end do
            end if
        endif
    end subroutine gauherJpred
    
    subroutine gauherJpredic(ss21,ss22,psig2,palpha,XbetapredRi,XbetapredDCi,hazRi,survDC,nrec0,survL,survU,survLT) !! 
    
!    use tailles
    use comon,only:typeof
    use donnees,only:w2,x2,w3,x3
    
    implicit none

    integer,intent(in)::nrec0
    double precision,intent(out)::ss21,ss22
    double precision::XbetapredRi,XbetapredDCi,var1
    double precision::auxfunca21,auxfunca22
    double precision,external ::func1pred2LogNicLogN,func2pred2LogNicLogN
    double precision,dimension(2)::survDC
    double precision::survLT,survL,survU !!
    double precision,dimension(nrec0)::hazRi,hazRi2
    double precision::psig2,palpha
    integer:: j

! gauss laguerre
! func1 est l integrant, ss le resultat de l integrale sur 0 ,  +infty

    hazRi2 = hazRi
    ss21=0.d0
    ss22=0.d0
    
    if (typeof == 0)then  
! Will be twice the average value of the function,since the ten
! weights (five numbers above each used twice) sum to 2.
        do j=1,20
            var1 = x2(j)
            auxfunca21 = func1pred2LogNicLogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!
            ss21 = ss21 + w2(j)*(auxfunca21)
            auxfunca22 = func2pred2LogNicLogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!
            ss22 = ss22 + w2(j)*(auxfunca22)
        end do
    else
        do j=1,32
            var1 = x3(j)
            auxfunca21 = func1pred2LogNicLogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!
            ss21 = ss21 + w3(j)*(auxfunca21)
            auxfunca22 = func2pred2LogNicLogN(var1,psig2,palpha,XbetapredRi,XbetapredDCi,survDC,survL,survU,survLT) !!
            ss22 = ss22 + w3(j)*(auxfunca22)
        end do
    end if
    end subroutine gauherJpredic
    
!======================================================================
!======================================================================

 !Pour tirer au sort aléatoirement dans une loi normale de moyenne m et d'écart-type s
    subroutine rnorm(m,s,res)
    
    double precision,intent(in)::m,s
    double precision,intent(out)::res
    double precision::alea1,alea2,UNIRAN
    double precision,parameter::pi=3.1415926536
    
    alea1 = UNIRAN()
    alea2 = UNIRAN()

    res = m + (s * sqrt(-2.d0*log(alea1)) * cos(2.d0*pi* alea2))
    
    end subroutine rnorm
