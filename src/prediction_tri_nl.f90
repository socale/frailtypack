
! ============================================== pred iction Joint

    
        subroutine predicttrinl(np,b,nz,nva10,nva20,nva30,nva40,nb_re0,which_random0,box_cox0,&
            nzyr,nzyd,link0,nst,typeof0,zi0,HIHOut, &
        ntimeAll,npred0,predTime,window,predtimerec,nrec0,nrecy0,yy0,matzy0,vaxpred0,&
        vaxdcpred0,vaxypred0,groupey,uniGroupey,nsujety, &
            nsujet,predAll1,predAlllow1,predAllhigh1, &
        icproba,nsample,movingwindow,timeAll,s_cag_id0,s_cag0)
    !model - type of a model 0- recur/survie, 1- longi/survie, 2-longi/recur/survie
            use donnees_indiv,only:nmescur,mu,ycurrent,b1,X2cur,Z1cur,nmescurr !it_cur,z2
            use comon,only:etaydc,sigmae,netadc,s_cag_id,s_cag,ut,utt,nva,link,npp,&
            nea,vey,nb1,netar,indic_Alpha,nva1,nva2, nva3,nva4,effet,zi,nz1,typeof,nb_re,alpha,&
            etayr,typeJoint,it,K_G0,K_D0,lambda,y0,mat,det,ziy,which_random,box_cox1,box_cox_par
            use lois_normales
            use prediction
        !     use ParametresPourParallelisation
            use optim
        implicit none
    
    
        integer::i,ii,iii,j,k,kk,jj
        integer,intent(in)::np,nz,nva10,nva20,nva30,nva40,nb_re0,nzyd,nzyr,nst,typeof0,ntimeAll,&
                                    icproba,movingwindow,nsujety,s_cag_id0,nsujet,link0,which_random0
        integer::npred0,nrec0,nsample,nrecy0
            double precision,dimension(nsujety),intent(in)::yy0
            double precision,dimension(npred0,nsujety)::yy_matrice
            double precision,dimension(2)::box_cox0
            integer,dimension(nsujety),intent(in)::groupey
        double precision,dimension(np),intent(in)::b
        double precision,dimension(nz+6),intent(in)::zi0
        double precision,dimension(np,np),intent(in)::HIHOut
            double precision,dimension(nsujet,nva10),intent(in)::vaxpred0
        double precision,dimension(npred0,nva20),intent(in)::vaxdcpred0
            double precision,dimension(nsujety,nva30+nva40),intent(in)::vaxypred0
        double precision,dimension(1,npred0)::XbetapredDC,XbetapredDCalea
            double precision,dimension(1,nsujety) :: XbetapredY,XbetapredYalea
            double precision,dimension(1,nsujet) :: XbetapredR,XbetapredRalea
            double precision,dimension(nsujety,2) :: matzy0
        integer,dimension(npred0)::nreci,nreci_all,nrecyi,uniGroupey
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
        double precision,dimension(1,nva30+nva40)::coefBetaYalea
        double precision,dimension(1,nva20)::coefBetadcalea
            double precision,dimension(1,nva10)::coefBetaRalea
        double precision,dimension(npred0,nrec0)::predtimerec
        double precision,dimension(npred0,nrec0+2)::predtimerec2
        double precision,dimension(1,nva20)::coefBetadc
            double precision,dimension(1,nva30+nva40)::coefBetay
            double precision,dimension(1,nva10)::coefBetaR
    
            integer :: ndim, restar,nf2,ier,jjj
    double precision::finddet,eps
       double precision,dimension((nb_re0+1)*(nb_re0+1+1)/2)::matv
       
       
       !---------------------------------------------------------------------------

            typejoint = 3
    s_cag_id = s_cag_id0
    
    s_cag = s_cag0
    
            link = link0
    netadc = nzyd
    netar = nzyr

    nea = nb_re0+1
    nb1 = nb_re0
      
        nb_re = nb_re0!+ (nb_re0*(nb_re0-1))/2.d0
    typeof = typeof0
    indic_Alpha = 1
    nva1 = nva10
    nva2 = nva20
    nva3 = nva30
    nva4 = nva40
    
     effet = 1
    box_cox1 = int(box_cox0(1))
        if(box_cox1.eq.1)then 
        box_cox_par = box_cox0(2)
        else 
        box_cox_par = int(box_cox0(2))
        end if
        
    which_random = which_random0
 
        allocate(zi(-2:nz+3),b1(np))
    b1(1:np) = b(1:np)
    
    zi(-2:nz+3) = zi0(1:nz+6)
    
    
    nz1= nz
    npp=np
    
    nva = nva1+nva2+nva3+nva4
            coefBetaR(1,:) = b((np-nva+1):(np-nva4-nva3-nva2))
    XbetapredR = matmul(coefBetaR,transpose(vaxpred0))
    
        coefBetaDC(1,:) = b((np-nva4-nva3-nva2+1):(np-nva4-nva3))
    
        XbetapredDC = matmul(coefBetadc,transpose(vaxdcpred0))
            coefBetaY(1,:) = b((np-nva4-nva3+1):np)
    
            XbetapredY = matmul(coefBetay,transpose(vaxypred0))
            !  XbetapredDC = matmul(coefBetadc,transpose(vaxdcpred0))
    

            allocate(survRi(nrec0+2),hazRi(nrec0+2))
    
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
    k=1
    
            ndim =  netar+1
    
    
            restar = 0
            nf2 = 1
    
    
       
                    allocate(Ut(nea,nea),Utt(nea,nea),mat(nea,nea),etaydc(netadc),etayr(netar))
    
    
                    ut = 0.d0
                    utt = 0.d0
            allocate(ycurrent(nrecy0),mu(nrecy0,2),ziy(nsujety,2))
            
            ziy = matzy0
    
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
    
    
                            do while ((groupey(k).eq.uniGroupey(i)).and.(k.lt.nsujety))
                            !       write(*,*)npred0,i,k,vaxypred0(k,2),predtimerec2(i,1)
                                    if(ziy(k,1).le.predtimerec2(i,1)) then ! check if measurement happened before prediction time
                                            nrecyi(i) = nrecyi(i)+1
                                    end if
                                    nreci_all(i) = nreci_all(i) + 1
                                    if(k.eq.(nsujety-1)) then
                                            if(ziy(k+1,1).le.predtimerec2(i,1)) then ! check if measurement happened before prediction time
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
        
                    ! les y jusqu'? predtimerec2
                    yy_matrice = 0.d0
                    it = 1
                    do kk=1,npred0
    
                            if(nrecyi(kk).gt.0) then
                            yy_matrice(kk,1:nrecyi(kk)) = yy0(it:(it+nrecyi(kk)-1))
    
                            end if
                            it = it+nreci_all(kk)
    
                    end do
    
    !       write(*,*)iii,'predtimred2'
    !       do i=1,npred0
    !               write(*,*)predtimerec2(i,1:nrec0+2)
    !       end do
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
    
        K_G0 =  b(np-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)
        K_D0 =  b(np-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)
        lambda =  b(np-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)
        y0 = b(np-nva-nb_re-1-netadc - netar-effet-indic_alpha)
        
    
    
    
    
                    alpha = b(np-nva-nb_re-1-netadc - netar)
    
                    sigmae = b(np-nva-nb_re)*b(np-nva-nb_re)
    
            etaydc(1:netadc) = b(np-nva-nb_re-netadc:np-nva-nb_re-1)
      
            etayr(1:netar) =b(np-nva-nb_re-netadc - netar:np-nva-nb_re-netadc -1)
           
    
            Ut = 0.d0
            Utt = 0.d0
    
                                do jj=1,nb1
                do k=1,jj
                
                    if(jj.eq.k) then
                   Ut(jj,k)=sqrt(b(np-nva-nb_re+k)**2.d0)
                  Utt(k,jj)=sqrt(b(np-nva-nb_re+k)**2.d0)
                 
                end if
                end do
            end do
    
    
                Ut(nea,nea) =  b(np-nva-nb_re-1-netadc - netar-1)
                Utt(nea,nea) =  b(np-nva-nb_re-1-netadc - netar-1)
    
    
    
     mat = matmul(ut,utt)
    

        jj=0
    ! jjj = 0
        do j=1,nea
        do k=j,nea
        jj=j+k*(k-1)/2
    !    jjj = jjj +1
    
        matv(jj)=mat(j,k)
        !  bb2vv(jjj)=bb2(j,k)
    
        end do
        end do
        ier = 0
        eps = 1.d-10
    
 
        call dsinvj(matv,nea,eps,ier)
    
        mat=0.d0
        do j=1,nea
                do k=1,nea
                            if (k.ge.j) then
                mat(j,k)=matv(j+k*(k-1)/2)
                else
                mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
                end do
                
            
    det = finddet(matmul(ut,utt),nea)
  
            !theta = b(np-nva1-nva2-1)*b(np-nva1-nva2-1)
            !alpha = b(np-nva1-nva2)
            it = 1
            do i=1,npred0
                            ycurrent  =0.d0
                            mu = 0.d0
                         !   z2 = 0.d0
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
    
                                    mu(1:nmescur,1) = matmul(vaxypred0(it:(it+nmescur-1),1:(nva3)),b((np-nva4-nva3+1):(np-nva4)))
                                    mu(1:nmescur,2) = matmul(vaxypred0(it:(it+nmescur-1),(nva3+1):(nva3+nva4)),b((np-nva4+1):np))
  
        
                             
                                     
    
                                    end if
              if(nb1.eq.1) call gauherPred_tri2_nl(ss11,1)
           
            if(nb1.eq.2) call gauherPred_tri3_nl(ss11,1)
            if(nb1.eq.3)call gauherPred_tri4_nl(ss11,1)
            if(nb1.eq.4)call gauherPred_tri5_nl(ss11,1)
        
             if(nb1.eq.1) call gauherPred_tri2_nl(ss12,2)
           
            if(nb1.eq.2)   call gauherPred_tri3_nl(ss12,2)
        if(nb1.eq.3)      call gauherPred_tri4_nl(ss12,2)
            if(nb1.eq.4)call gauherPred_tri5_nl(ss12,2)
        !    write(*,*)nb1,ss11,ss12
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
    
            K_G0 =  balea(j,np-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)
        K_D0 =  balea(j,np-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)
        lambda =  balea(j,np-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)
        y0 = balea(j,np-nva-nb_re-1-netadc - netar-effet-indic_alpha)
        
    
    
                                    alpha = balea(j,np-nva-nb_re-1-netadc - netar)
    
                    sigmae =balea(j,np-nva-nb_re)*balea(j,np-nva-nb_re)
    
    
            etaydc(1:netadc) = balea(j,np-nva-nb_re-netadc:np-nva-nb_re-1)
           
            etayr(1:netar) =balea(j,np-nva-nb_re-netadc - netar:np-nva-nb_re-netadc -1)
     
    
            Ut = 0.d0
            Utt = 0.d0
    
                               do jj=1,nb1
                do k=1,jj
                
                    if(jj.eq.k) then
                   Ut(jj,k)=sqrt(balea(j,np-nva-nb_re+k)**2.d0)
                  Utt(k,jj)=sqrt(balea(j,np-nva-nb_re+k)**2.d0)
                 
                end if
                end do
            end do
    
    
                Ut(nea,nea) =  balea(j,np-nva-nb_re-1-netadc - netar-1)
                Utt(nea,nea) =  balea(j,np-nva-nb_re-1-netadc - netar-1)
    
    
     mat = matmul(ut,utt)
    

        jj=0
    ! jjj = 0
        do jjj=1,nea
        do k=j,nea
        jj=jjj+k*(k-1)/2
    !    jjj = jjj +1
    
        matv(jj)=mat(jjj,k)
        !  bb2vv(jjj)=bb2(j,k)
    
        end do
        end do
        ier = 0
        eps = 1.d-10
    
 
        call dsinvj(matv,nea,eps,ier)
    
        mat=0.d0
        do jjj=1,nea
                do k=1,nea
                            if (k.ge.jjj) then
                mat(jjj,k)=matv(jjj+k*(k-1)/2)
                else
                mat(jjj,k)=matv(k+jjj*(jjj-1)/2)
                end if
            end do
                end do
    det = finddet(matmul(ut,utt),nea)
  
    
    
    
            !theta = b(np-nva1-nva2-1)*b(np-nva1-nva2-1)
            !alpha = b(np-nva1-nva2)
            it = 1
            do i=1,npred0
                      !      z2 = 0.d0
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
                    mu(1:nmescur,1) = matmul(vaxypred0(it:(it+nmescur-1),1:(nva3)),balea(j,(np-nva4-nva3+1):(np-nva4)))
                    mu(1:nmescur,2) = matmul(vaxypred0(it:(it+nmescur-1),(nva3+1):(nva3+nva4)),balea(j,(np-nva4+1):np))
                    end if
  
              
    
                 if(nb1.eq.2) call gauherPred_tri2_nl(ss11,1)
           
            if(nb1.eq.2) call gauherPred_tri3_nl(ss11,1)
            if(nb1.eq.3)call gauherPred_tri4_nl(ss11,1)
            if(nb1.eq.4)call gauherPred_tri5_nl(ss11,1)
        
             if(nb1.eq.1) call gauherPred_tri2_nl(ss12,2)
           
            if(nb1.eq.2)   call gauherPred_tri3_nl(ss12,2)
        if(nb1.eq.3)      call gauherPred_tri4_nl(ss12,2)
            if(nb1.eq.4)call gauherPred_tri5_nl(ss12,2)

    
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
        !write(*,*)"here"
            !stop
            deallocate(mu,ycurrent,ziy,survRi,hazRi,ut,utt,mat)
            deallocate(vey,X2cur,Z1cur,b1,zi,etaydc,etayr)
            !   deallocate(ut)
        end subroutine predicttrinl
    
    
    
                               !***********************************
            !********* Gauss-Hermit pour la dimension 5 - mod?le trviarie b_G0,b_D0, b_lambda, b_y0, v
            !*************************************
    
            SUBROUTINE gauherPred_tri5_nl(ss,choix)
    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only : typeof !auxig,netadc
            use donnees_indiv,only : frailpol4 !numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca
        external::gauherPred_tri4_nl
        integer::j
    
        ss=0.d0
        if (typeof.eq.0) then
            do j=1,20
    
                frailpol4 = x2(j)*sqrt(2.d0)
                            call gauherPred_tri4_nl(auxfunca,choix)
                    ss = ss+w2(j)*(auxfunca)*sqrt(2.d0)
    
            end do
        
        else
    
                            do j=1,32
                            frailpol4= x3(j)*sqrt(2.d0)
                            call gauherPred_tri4_nl(auxfunca,choix)
    
                            ss = ss+w3(j)*(auxfunca)*sqrt(2.d0)
                            end do
    
    
        endif
    
        return
    
        END SUBROUTINE gauherPred_tri5_nl
    
    
    
                                !***********************************
            !********* Gauss-Hermit pour la dimension 4 - mod?le trviarie b_10,b_11, b_12, v
            !*************************************
    
            SUBROUTINE gauherPred_tri4_nl(ss,choix)
    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only:typeof !auxig,netadc
            use donnees_indiv,only : frailpol3 !numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca
        external::gauherPred_tri3_nl
        integer::j
    
        ss=0.d0
        if (typeof.eq.0) then
            do j=1,20
    
                frailpol3 = x2(j)*sqrt(2.d0)
                            call gauherPred_tri3_nl(auxfunca,choix)
                    ss = ss+w2(j)*(auxfunca)*sqrt(2.d0)
    
            end do
            
        else
    
                            do j=1,32
                            frailpol3= x3(j)*sqrt(2.d0)
                            call gauherPred_tri3_nl(auxfunca,choix)
    
                            ss = ss+w3(j)*(auxfunca)*sqrt(2.d0)
                            end do
    
    
        endif
    
        return
    
        END SUBROUTINE gauherPred_tri4_nl
    
    
                                    !***********************************
            !********* Gauss-Hermit pour la dimension 3 - mod?le trviarie b_10,b_11, v
            !*************************************
    
            SUBROUTINE gauherPred_tri3_nl(ss,choix)
    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only:typeof !auxig,netadc
            use donnees_indiv,only : frailpol2 !numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca
        external::gauherPred_tri2_nl
        integer::j
    
        ss=0.d0
        if (typeof.eq.0) then
            do j=1,20
    
                frailpol2 = x2(j)*sqrt(2.d0)
                            call gauherPred_tri2_nl(auxfunca,choix)
                    ss = ss+w2(j)*(auxfunca)*sqrt(2.d0)
    
            end do
                
        else
    
                            do j=1,32
                            frailpol2 = x3(j)*sqrt(2.d0)
                            call gauherPred_tri2_nl(auxfunca,choix)
    
                            ss = ss+w3(j)*(auxfunca)*sqrt(2.d0)
                            end do
    
    
        endif
    
        return
    
        END SUBROUTINE gauherPred_tri3_nl
    
    !===============================================
            SUBROUTINE gauherPred_tri2_nl(ss,choix)
    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only:typeof !auxig
            use donnees_indiv,only : frailpol !numpat
        Implicit none
    
        double precision,intent(out)::ss
            integer,intent(in)::choix
        double precision::auxfunca
        external::gauherPred_tri_nl
        integer::j
    
    
        ss=0.d0
        if (typeof.eq.0) then
            do j=1,20
                frailpol = x2(j)*sqrt(2.d0)
                                    call gauherPred_tri_nl(auxfunca,choix)
                                    ss = ss+w2(j)*(auxfunca)*sqrt(2.d0)
            end do
            
        else
                            do j=1,32
                                    frailpol = x3(j)*sqrt(2.d0)
                                    call gauherPred_tri_nl(auxfunca,choix)
                                    ss = ss+w3(j)*(auxfunca)*sqrt(2.d0)
                    end do
    
        endif
    
        return
    
        END SUBROUTINE gauherPred_tri2_nl
    
    
            SUBROUTINE gauherPred_tri_nl(ss,choix)
    
        use tailles
        use donnees,only:x2,w2,x3,w3
        use comon,only:typeof,which_random !auxig,nb1
            use donnees_indiv,only : frailpol,frailpol2,frailpol3,frailpol4
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca,func1pred1GHtri_nl,func2pred1GHtri_nl,&
            func1pred2GHtri_nl,func2pred2GHtri_nl,&
            func1pred3GHtri_nl,func2pred3GHtri_nl,&
            func1pred4GHtri_nl,func2pred4GHtri_nl,&
            func1pred5GHtri_nl,func2pred5GHtri_nl,&
            func1pred6GHtri_nl,func2pred6GHtri_nl,&
            func1pred7GHtri_nl,func2pred7GHtri_nl,&
            func1pred8GHtri_nl,func2pred8GHtri_nl,&
            func1pred9GHtri_nl,func2pred9GHtri_nl,&
            func1pred10GHtri_nl,func2pred10GHtri_nl,&
            func1pred11GHtri_nl,func2pred11GHtri_nl,&
            func1pred12GHtri_nl,func2pred12GHtri_nl,&
            func1pred13GHtri_nl,func2pred13GHtri_nl,&
            func1pred14GHtri_nl,func2pred14GHtri_nl,&
            func1pred15GHtri_nl,func2pred15GHtri_nl
        
        integer::j

        ss=0.d0
            auxfunca = 0.d0
        if (typeof.eq.0) then
            do j=1,20
                 select case(which_random)
                    case(1)
                        if(choix.eq.1) then
                            auxfunca=func1pred1GHtri_nl(frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred1GHtri_nl(frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(2)
                        if(choix.eq.1) then
                            auxfunca=func1pred2GHtri_nl(frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred2GHtri_nl(frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(3)
                        if(choix.eq.1) then
                            auxfunca=func1pred3GHtri_nl(frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred3GHtri_nl(frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(4)
                        if(choix.eq.1) then
                            auxfunca=func1pred4GHtri_nl(frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred4GHtri_nl(frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(5)
                        if(choix.eq.1) then
                            auxfunca=func1pred5GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred5GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(6)
                        if(choix.eq.1) then
                            auxfunca=func1pred6GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred6GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(7)
                        if(choix.eq.1) then
                            auxfunca=func1pred7GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred7GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(8)
                        if(choix.eq.1) then
                            auxfunca=func1pred8GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred8GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(9)
                        if(choix.eq.1) then
                            auxfunca=func1pred9GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred9GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(10)
                        if(choix.eq.1) then
                            auxfunca=func1pred10GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred10GHtri_nl(frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(11)
                        if(choix.eq.1) then
                            auxfunca=func1pred11GHtri_nl(frailpol3,frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred11GHtri_nl(frailpol3,frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(12)
                        if(choix.eq.1) then
                            auxfunca=func1pred12GHtri_nl(frailpol3,frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred12GHtri_nl(frailpol3,frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(13)
                        if(choix.eq.1) then
                            auxfunca=func1pred13GHtri_nl(frailpol3,frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred13GHtri_nl(frailpol3,frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(14)
                        if(choix.eq.1) then
                            auxfunca=func1pred14GHtri_nl(frailpol3,frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred14GHtri_nl(frailpol3,frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                    case(15)
                        if(choix.eq.1) then
                            auxfunca=func1pred15GHtri_nl(frailpol4,frailpol3,frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred15GHtri_nl(frailpol4,frailpol3,frailpol2,frailpol,x2(j)*sqrt(2.d0))
                        endif
                end select    
            
               ss = ss+w2(j)*(auxfunca)*sqrt(2.d0)
            end do
    
        else
            do j=1,32
                        select case(which_random)
                    case(1)
                        if(choix.eq.1) then
                            auxfunca=func1pred1GHtri_nl(frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred1GHtri_nl(frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(2)
                        if(choix.eq.1) then
                            auxfunca=func1pred2GHtri_nl(frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred2GHtri_nl(frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(3)
                        if(choix.eq.1) then
                            auxfunca=func1pred3GHtri_nl(frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred3GHtri_nl(frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(4)
                        if(choix.eq.1) then
                            auxfunca=func1pred4GHtri_nl(frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred4GHtri_nl(frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(5)
                        if(choix.eq.1) then
                            auxfunca=func1pred5GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred5GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(6)
                        if(choix.eq.1) then
                            auxfunca=func1pred6GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred6GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(7)
                        if(choix.eq.1) then
                            auxfunca=func1pred7GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred7GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(8)
                        if(choix.eq.1) then
                            auxfunca=func1pred8GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred8GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(9)
                        if(choix.eq.1) then
                            auxfunca=func1pred9GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred9GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(10)
                        if(choix.eq.1) then
                            auxfunca=func1pred10GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred10GHtri_nl(frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(11)
                        if(choix.eq.1) then
                            auxfunca=func1pred11GHtri_nl(frailpol3,frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred11GHtri_nl(frailpol3,frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(12)
                        if(choix.eq.1) then
                            auxfunca=func1pred12GHtri_nl(frailpol3,frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred12GHtri_nl(frailpol3,frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(13)
                        if(choix.eq.1) then
                            auxfunca=func1pred13GHtri_nl(frailpol3,frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred13GHtri_nl(frailpol3,frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(14)
                        if(choix.eq.1) then
                            auxfunca=func1pred14GHtri_nl(frailpol3,frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred14GHtri_nl(frailpol3,frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                    case(15)
                        if(choix.eq.1) then
                            auxfunca=func1pred15GHtri_nl(frailpol4,frailpol3,frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        else if(choix.eq.2) then
                            auxfunca = func2pred15GHtri_nl(frailpol4,frailpol3,frailpol2,frailpol,x3(j)*sqrt(2.d0))
                        endif
                end select    
                ss = ss+w3(j)*(auxfunca)*sqrt(2.d0)
                end do
            endif
    
    
    
        return
    
        END SUBROUTINE gauherPred_tri_nl
    
    
    
    
    
    
     !=========================
    ! Prediction  : 1+1 effet aleatoire (case 1)
    !=========================
    double precision function func1pred1GHtri_nl(frail,frail2)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc$
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
             Xea22(1) = frail
            Xea22(2) = frail2
         
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
      !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred1GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        frail2*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+frail2*alpha))) &
                                        * exp((frail2+etayr(1)*frail))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail2+etayr(1)*frail))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred1GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred1GHtri_nl(frail,frail2)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
              Xea22(1) = frail
            Xea22(2) = frail2
          
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
          !   if(link.eq.1) then
            func2pred1GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail&
                                    +frail2*alpha)) ) &
                                        * exp((frail2+etayr(1)*frail))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail2+etayr(1)*frail))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred1GHtri_nl
    
    
    
    
    
     !=========================
    ! Prediction  : 1+1 effet aleatoire (case 2)
    !=========================
    double precision function func1pred2GHtri_nl(frail,frail2)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
             Xea22(1) = frail
            Xea22(2) = frail2
          
   
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)

       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
      !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred2GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail&
                                        +frail2*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail&
                                        +frail2*alpha))) &
                                        * exp((frail2+etayr(1)*frail))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail2+etayr(1)*frail))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
                                        
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred2GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred2GHtri_nl(frail,frail2)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
             Xea22(1) = frail
            Xea22(2) = frail2
           
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred2GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail&
                                        +frail2*alpha)) ) &
                                        * exp((frail2+etayr(1)*frail))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail2+etayr(1)*frail))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred2GHtri_nl
    
    
    
    
    
    
     !=========================
    ! Prediction  : 1+1 effet aleatoire (case 3)
    !=========================
    double precision function func1pred3GHtri_nl(frail,frail2)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
             Xea22(1) = frail
            Xea22(2) = frail2
          
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(1)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
     !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred3GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        frail2*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        frail2*alpha))) &
                                        * exp((frail2+etayr(1)*frail))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail2+etayr(1)*frail))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred3GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred3GHtri_nl(frail,frail2)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea22(1) = frail
            Xea22(2) = frail2
          
          uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(1)-lambda+mu(1:nmescur,2))*&
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred3GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        frail2*alpha)) ) &
                                        * exp((frail2+etayr(1)*frail))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail2+etayr(1)*frail))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred3GHtri_nl
    
    
    
    
     !=========================
    ! Prediction  : 1+1 effet aleatoire (case 4)
    !=========================
    double precision function func1pred4GHtri_nl(frail,frail2)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
             Xea22(1) = frail
            Xea22(2) = frail2
          
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda-Xea22(1)+mu(1:nmescur,2))*&
        (dexp(-dexp(lambda+Xea22(1))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
        !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred4GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+frail2*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+frail2*alpha))) &
                                        * exp((frail2+etayr(1)*frail))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail2+etayr(1)*frail))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred4GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred4GHtri_nl(frail,frail2)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea22(1) = frail
            Xea22(2) = frail2
               
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda-Xea22(1)+mu(1:nmescur,2))*&
        (dexp(-dexp(lambda+Xea22(1))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
        !   if(link.eq.1) then
            func2pred4GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        frail2*alpha)) ) &
                                        * exp((frail2+etayr(1)*frail))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail2+etayr(1)*frail))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred4GHtri_nl
    
    
    
    
     !=========================
    ! Prediction  : 2+1 effet aleatoire (case 5)
    !=========================
    double precision function func1pred5GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
  
   !if(frail.le.-7.62.and.frail2.eq.-7.62.and.frail3.eq.-7.62) write(*,*)mat
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
      !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred5GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha))) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail+&
                                        etayr(2)*frail2))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

            return
    
        end function func1pred5GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred5GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
                Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred5GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) ) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail&
                                        +etayr(2)*frail2))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred5GHtri_nl
    
    
    
    
    
    
     !=========================
    ! Prediction  : 2+1 effet aleatoire (case 6)
    !=========================
    double precision function func1pred6GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
              Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(2)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
     !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred6GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha))) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail+&
                                        etayr(2)*frail2))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred6GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred6GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(2)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred6GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) ) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail&
                                        +etayr(2)*frail2))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred6GHtri_nl
    
    
    
    
    
    
     !=========================
    ! Prediction  : 2+1 effet aleatoire (case 7)
    !=========================
    double precision function func1pred7GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
    
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda-Xea22(2)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(2))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
         !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred7GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha))) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail+&
                                        etayr(2)*frail2))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred7GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred7GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda-Xea22(2)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(2))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred7GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) ) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail&
                                        +etayr(2)*frail2))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred7GHtri_nl
    
    
    
    
    
     !=========================
    ! Prediction  : 2+1 effet aleatoire (case 8)
    !=========================
    double precision function func1pred8GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
    
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(2)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
       !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred8GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha))) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail+&
                                        etayr(2)*frail2))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred8GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred8GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(2)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
        !   if(link.eq.1) then
            func2pred8GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) ) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail&
                                        +etayr(2)*frail2))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred8GHtri_nl
    
    
    
    
     !=========================
    ! Prediction  : 2+1 effet aleatoire (case 9)
    !=========================
    double precision function func1pred9GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
    
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda-Xea22(2)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(2))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
       !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred9GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha))) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail+&
                                        etayr(2)*frail2))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred9GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred9GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda-Xea22(2)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(2))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred9GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) ) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail&
                                        +etayr(2)*frail2))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred9GHtri_nl
    
    
    
    
    
     !=========================
    ! Prediction  : 2+1 effet aleatoire (case 10)
    !=========================
    double precision function func1pred10GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
    
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(1)-lambda-Xea22(2)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(2))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
       !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred10GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha))) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail+&
                                        etayr(2)*frail2))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred10GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred10GHtri_nl(frail,frail2,frail3)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(1)-lambda-Xea22(2)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(2))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred10GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+frail3*alpha)) ) &
                                        * exp((frail3+etayr(1)*frail+etayr(2)*frail2))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail3+etayr(1)*frail&
                                        +etayr(2)*frail2))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred10GHtri_nl
    
    
    
    
     !=========================
    ! Prediction  : 3+1 effet aleatoire (case 11)
    !=========================
    double precision function func1pred11GHtri_nl(frail,frail2,frail3,frail4)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3,frail4
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea2(4,1) = frail4
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            Xea22(4) = frail4
    
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(3)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
     !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred11GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha))) &
                                        * exp((frail4+etayr(1)*frail+etayr(2)*frail2+etayr(3)*frail3))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail4+etayr(1)*frail+&
                                        etayr(2)*frail2+etayr(3)*frail3))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred11GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred11GHtri_nl(frail,frail2,frail3,frail4)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3,frail4
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
             Xea2(4,1) = frail4
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            Xea22(4) = frail4
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(3)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred11GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha)) ) &
                                        * exp((frail4+etayr(1)*frail+etayr(2)*frail2+etayr(3)*frail3))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail4+etayr(1)*frail&
                                        +etayr(2)*frail2+etayr(3)*frail3))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred11GHtri_nl
    
    
    
    
    
    
     !=========================
    ! Prediction  : 3+1 effet aleatoire (case 12)
    !=========================
    double precision function func1pred12GHtri_nl(frail,frail2,frail3,frail4)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3,frail4
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea2(4,1) = frail4
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            Xea22(4) = frail4
    
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda-Xea22(3)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(3))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
      !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred12GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha))) &
                                        * exp((frail4+etayr(1)*frail+etayr(2)*frail2+etayr(3)*frail3))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail4+etayr(1)*frail+&
                                        etayr(2)*frail2+etayr(3)*frail3))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred12GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred12GHtri_nl(frail,frail2,frail3,frail4)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3,frail4
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
             Xea2(4,1) = frail4
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            Xea22(4) = frail4
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0-lambda-Xea22(3)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(3))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred12GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha)) ) &
                                        * exp((frail4+etayr(1)*frail+etayr(2)*frail2+etayr(3)*frail3))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail4+etayr(1)*frail&
                                        +etayr(2)*frail2+etayr(3)*frail3))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred12GHtri_nl
    
    
    
    
     !=========================
    ! Prediction  : 3+1 effet aleatoire (case 13)
    !=========================
    double precision function func1pred13GHtri_nl(frail,frail2,frail3,frail4)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3,frail4
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea2(4,1) = frail4
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            Xea22(4) = frail4
    
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(2)-lambda-Xea22(3)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(3))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
      !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred13GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha))) &
                                        * exp((frail4+etayr(1)*frail+etayr(2)*frail2+etayr(3)*frail3))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail4+etayr(1)*frail+&
                                        etayr(2)*frail2+etayr(3)*frail3))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred13GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred13GHtri_nl(frail,frail2,frail3,frail4)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3,frail4
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
             Xea2(4,1) = frail4
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            Xea22(4) = frail4
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(2)-lambda-Xea22(3)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(3))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred13GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha)) ) &
                                        * exp((frail4+etayr(1)*frail+etayr(2)*frail2+etayr(3)*frail3))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail4+etayr(1)*frail&
                                        +etayr(2)*frail2+etayr(3)*frail3))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred13GHtri_nl
    
    
    
    
    
     !=========================
    ! Prediction  : 3+1 effet aleatoire (case 14)
    !=========================
    double precision function func1pred14GHtri_nl(frail,frail2,frail3,frail4)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3,frail4
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea2(4,1) = frail4
              Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            Xea22(4) = frail4
        
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(2)-lambda-Xea22(3)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(3))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
      !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred14GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha))) &
                                        * exp((frail4+etayr(1)*frail+etayr(2)*frail2+etayr(3)*frail3))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail4+etayr(1)*frail+&
                                        etayr(2)*frail2+etayr(3)*frail3))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred14GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred14GHtri_nl(frail,frail2,frail3,frail4)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey,
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3,frail4
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
             Xea2(4,1) = frail4
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            Xea22(4) = frail4
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(2)-lambda-Xea22(3)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(3))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
          !   if(link.eq.1) then
            func2pred14GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+frail4*alpha)) ) &
                                        * exp((frail4+etayr(1)*frail+etayr(2)*frail2+etayr(3)*frail3))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail4+etayr(1)*frail&
                                        +etayr(2)*frail2+etayr(3)*frail3))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred14GHtri_nl
    
    
    
    
    
    
    
     !=========================
    ! Prediction  : 4+1 effet aleatoire (case 15)
    !=========================
    double precision function func1pred15GHtri_nl(frail,frail2,frail3,frail4,frail5)
        ! calcul de l integrant (numerateur de la fonction de prediction)
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& ! !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3,frail4,frail5
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
    

            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea2(4,1) = frail4
            Xea2(5,1) = frail5
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            Xea22(4) = frail4
            Xea22(5) = frail5
    
    
    
                    uii = matmul(Xea22,mat)
    
    
                    uiiui=matmul(uii,Xea2)
    
       if(nmescur.gt.0) then
    
         mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(3)-lambda-Xea22(4)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(4))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
        
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag =  prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
                    else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
     
            yscalar = dsqrt(yscalar)
     !     if(link.eq.1) then        only for random effects link, no current level yet
            func1pred15GHtri_nl = (survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+etaydc(4)*frail4+frail5*alpha)) &
                                        - survDC(2)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+etaydc(4)*frail4+frail5*alpha))) &
                                        * exp((frail5+etayr(1)*frail+etayr(2)*frail2+etayr(3)*frail3&
                                        +etayr(4)*frail4))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail5+etayr(1)*frail+&
                                        etayr(2)*frail2+etayr(3)*frail3+etayr(4)*frail4))) &
                                        * dexp(-(yscalar**2.d0)/(2.d0*sigmae))*exp(prod_cag)&
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
       !     else if(link.eq.2) then
       ! to do
      !       end if

    
            return
    
        end function func1pred15GHtri_nl
    
    
            !*************************************************
        double precision  function func2pred15GHtri_nl(frail,frail2,frail3,frail4,frail5)
        ! calcul de l integrant (denominateur de la fonction de prediction)
            use optim
                    use comon,only:etaydc,sigmae,etayr,& !netadc,netar
                    s_cag,s_cag_id,alpha,box_cox1,box_cox_par,& ! !ut,utt,link,npp
                    nb1,mat,det,K_G0,K_D0,lambda,y0,ziy,it,nea !nva3,nva4,vey
            use donnees_indiv,only:nmescur,mu,ycurrent,nmescurr !b1,x2cur,z1cur,z2
    
            use prediction
        implicit none
    
    
    double precision,intent(in)::frail,frail2,frail3,frail4,frail5
            double precision,dimension(nmescur,1)::mu1
            double precision :: yscalar,prod_cag !eps,finddet
            integer :: k !j,jj,ier
            double precision,dimension((nb1+1),1)::  Xea2
            double precision,dimension((nb1+1)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,external::survdcCM_pred,survRCM_pred
        double precision,parameter::pi=3.141592653589793d0
            logical :: upper
            double precision :: alnorm
    !double precision :: resultdc,resultR,abserr,resabs,resasc
    
            upper = .false.
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
             Xea2(4,1) = frail4
             Xea2(5,1) = frail5
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
            Xea22(4) = frail4
            Xea22(5) = frail5
            
       
    
                    uii = matmul(Xea22,mat)
          
    
                    uiiui=matmul(uii,Xea2)
    
    
    
        if(nmescur.gt.0) then
    
                 mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy(it:(it+nmescur-1),1)+&
            ziy(it:(it+nmescur-1),2)*dexp(K_D0+Xea22(3)-lambda-Xea22(4)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(4))*ziy(it:(it+nmescur-1),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        end if
    
            yscalar = 0.d0
                            prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
    
    
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
        
    
            yscalar = dsqrt(yscalar)
    
         !   if(link.eq.1) then
            func2pred15GHtri_nl = ((survDC(1)**(exp(XbetapredDCi+etaydc(1)*frail+&
                                        etaydc(2)*frail2+etaydc(3)*frail3+etaydc(4)*frail4+frail5*alpha)) ) &
                                        * exp((frail5+etayr(1)*frail+etayr(2)*frail2+etayr(3)*frail3&
                                        +etaydc(4)*frail4))**nmescurr &
                                        * (survRi(1)**( exp(XbetapredRi+frail5+etayr(1)*frail&
                                        +etayr(2)*frail2+etayr(3)*frail3+etayr(4)*frail4))) &
                                        * exp(-(yscalar**2.d0)/(2.d0*sigmae)))*exp(prod_cag) &
                                        *dexp(-uiiui(1)/2.d0)/(det*2.d0*pi)**(nea/2.d0)
         !   else if(link.eq.2) then
         
         !   end if
    
    
            return
    
        end function func2pred15GHtri_nl
    
    
    
    
    
    
    
    
    