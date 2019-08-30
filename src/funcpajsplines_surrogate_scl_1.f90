


!========================          FUNCPAJ_SPLINES         ====================
    double precision function funcpajsplines_surrogate_1(b,np,id,thi,jd,thj,k0)
    
    use tailles
    use comon
    use residusM
    use var_surrogate ! fichier Aparameters_scl.f90
    use fonction_A_integrer ! integrant (fichier Integrant_scl.f90)
    use GaussHermi_mult ! pour la fonction d'integration (fichier Integrale_mult_scl.f90)
    use monteCarlosMult_Gaus ! pour integration par monte carlo (fichier GuassHermi_mult_scl.f90)
        
    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :
    
    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(2),intent(in)::k0
    double precision,intent(in)::thi,thj
    
    integer::n,i,j,k,vj,ig,choix,l,vcdiag,nsujet_trial
    integer,dimension(ngmax)::cpt
    double precision::pe1,pe2,inv,som1,som2,res,vet,vet2,h1 !som,inc
    double precision,dimension(3):: resultatInt
    
    double precision,dimension(-2:npmax):: the1,the2
    double precision,dimension(np)::bh
    double precision,dimension(ngmax)::res2,res1dc,res2dc &
    ,res3dc,integrale1,integrale2
!AD: for death,change dimension 
    double precision,dimension(ndatemax)::dut1
    double precision,dimension(ndatemaxdc)::dut2
!AD:end
    double precision,dimension(0:ndatemax)::ut1
    double precision,dimension(0:ndatemaxdc)::ut2
    !double precision,dimension(:),allocatable::frail
    double precision::logGammaJ,pourgam !c3,c4,int
    double precision,dimension(ntrials)::integrale3

!    !print*,'debut funcpa'
    choix=0
    ig=0
    k=0
    vj=0
    n=0
    j=0
    ut1=0.d0
    ut2=0.d0
    dut2=0.d0
    dut1=0.d0        
    do i=1,np
        bh(i)=b(i)
    end do 

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj    
    n = (np-nva-nparamfrail)/nst
    !n = (np-nva-effet-indic_ALPHA)/nst
    !!print*,"nparamfrail=",nparamfrail,"effet=",effet,"indic_ALPHA=",indic_ALPHA
    !!print*,"(np-nva-nparamfrail)/nst",(np-nva-nparamfrail)/nst,"n=",n

    do i=1,n
        the1(i-3)=(bh(i))*(bh(i))
        j = n+i 
        if (nst.eq.2) then
            the2(i-3)=(bh(j))*(bh(j))
        endif
    end do
    
    
    if(effet.eq.1) then
        if(logNormal==1)then
            theta2 = bh(np-nva-nparamfrail+indice_eta+indice_theta)**2.d0
            !alpha  = bh(np-nva-nparamfrail+indice_eta+indice_theta+indic_alpha)
            !sigma2 = bh(np-nva-nparamfrail+indice_eta+indice_theta+indic_alpha+indice_sigma)**2
            sig2=theta2 ! je fais appel a sig2 car c'est la variable utilisee dans la suite des procedures pour le joint classique
        else
            theta = bh(np-nva-nparamfrail+indice_eta+indice_theta)**2.d0
            sig2=theta ! je fais appel a sig2 car c'est la variable utilisee dans la suite des procedures pour le joint classique
        endif
        
        if(indice_eta==0)then
            eta=1.d0! on fixe eta a 1
        else
            eta = bh(np-nva-nparamfrail+indice_eta)
        endif
        !!print*,"theta2=",theta2
    endif
!----------  calcul de ut1(ti) et ut2(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

!AD:modify
    dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
    dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
    
    ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))
    ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
    
    ut1(0) = 0.d0
    ut2(0) = 0.d0
     
!//// NEW AMADOU vvv :
!--- strate1
    som1 = 0.d0
    vj = 0
    do i=2,ndate-1
        do k = 2,n-2
            if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                som1 = som1 + the1(j-4)
                vj  = j
                endif   
            endif
        end do 
    
        ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i)) &
        +(the1(j-1)*im1(i))+(the1(j)*im(i))
        
        dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i)) &
        +(the1(j-1)*mm1(i))+(the1(j)*mm(i))

    end do
                   
!--- strate2
    vj = 0
    som2 = 0.d0

    do i=2,ndatedc-1
        do k = 2,n-2
            if (((datedc(i)).ge.(zi(k-1))).and.(datedc(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                som2 = som2 + the2(j-4)
                vj  = j
                endif   
            endif
        end do

        if(nst.eq.2)then
            ut2(i) = som2 +(the2(j-3)*im3dc(i))+(the2(j-2)*im2dc(i)) &
            +(the2(j-1)*im1dc(i))+(the2(j)*imdc(i))
            dut2(i) = (the2(j-3)*mm3dc(i))+(the2(j-2)*mm2dc(i)) &
            +(the2(j-1)*mm1dc(i))+(the2(j)*mmdc(i))
        endif
            
    end do

!-------------fin strate2  
    i = n-2
    h1 = (zi(i)-zi(i-1))
    
    ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
    dut1(ndate) = (4.d0*the1(i-1)/h1)

        
    ut2(ndatedc)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)!am the1(i-4)
    dut2(ndatedc) = (4.d0*the2(i-1)/h1)
         
!    !print*,ndatemaxdc
!    !print*,'ut2',ut2
    
    
!//// fin NEW AMADOU-vvv
!AD:end
!-------------------------------------------------------
!---------- calcul de la vraisemblance ------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------

    do k=1,ng
        res1(k) = 0.d0
        res2(k) = 0.d0
        res3(k) = 0.d0
        res1dc(k) = 0.d0
        res2dc(k) = 0.d0
        res3dc(k) = 0.d0
        cpt(k) = 0
        integrale1(k) = 0.d0
        integrale2(k) = 0.d0
        !integrale3_(k) = 0.d0
        aux1(k)=0.d0
        aux2(k)=0.d0
    end do
    
    integrale3 = 0.d0
    const_res1=0.d0
    const_aux1=0.d0
    const_res4=0.d0
    const_res5=0.d0
    res2s=0.d0
    res2_dcs=0.d0

!*******************************************         
!-----avec un effet aleatoire dans le modele
!*********************************************

    inv = 1.d0/sigma2

!ccccccccccccccccccccccccccccccccccccccccc
!     pour les donnees recurrentes
!ccccccccccccccccccccccccccccccccccccccccc

    !!print*,"g=",g

    do i=1,nsujet 
        cpt(g(i))=cpt(g(i))+1  
        if(nva1.gt.0)then
            vet = 0.d0   
            do j=1,nva1
                vet =vet + bh(np-nva+j)*dble(ve(i,j))
!                  !write(*,*)'*** funcpaj_splines vet',vet,ve(i,j),i,j
            end do
        else
            vet=1.d0
        endif
            
        if((c(i).eq.1))then
            res2s(pourtrial(i)) = res2s(pourtrial(i))+dlog(dut1(nt1(i)))+vet
        endif  
        if ((res2s(pourtrial(i)).ne.res2s(pourtrial(i))).or.(abs(res2s(pourtrial(i))).ge. 1.d30)) then
            funcpajsplines_surrogate_1=-1.d9
            goto 123
        end if    
        
        ! scl pour calcul integrale
        !const_res1(pourtrial(i))=const_res1(pourtrial(i))+ut1(nt1(i))*dexp(vet)
        !if ((const_res1(pourtrial(i)).ne.const_res1(pourtrial(i))).or.(abs(const_res1(pourtrial(i))).ge. 1.d30)) then
        !    funcpajsplines=-1.d9
        !    goto 123
        !end if
        
        const_res4(g(i)) = const_res4(g(i))+ut1(nt1(i))* dexp(vet)
        if ((const_res4(g(i)).ne.const_res4(g(i))).or.(abs(const_res4(g(i))).ge. 1.d30)) then !scl si risque cumule >10^30
            funcpajsplines_surrogate_1=-1.d9
            goto 123
        end if
    end do
!    !print*,'const_res1',const_res1
!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces 
!ccccccccccccccccccccccccccccccccccccccccc 

    do k=1,ng  
        if(nva2.gt.0)then
            vet2 = 0.d0   
            do j=1,nva2
                vet2 =vet2 + bh(np-nva2+j)*dble(vedc(k,j))
            end do
            !vet2 = dexp(vet2)!scl
        else
            vet2=1.d0
        endif
        !!print*,"cdc=",cdc
        !stop
        if(cdc(k).eq.1)then
            !!print*,"k=",k,"size(nt1dc)",size(nt1dc),"nt1dc(k)=",nt1dc(k),"dut2(nt1dc(k))=",dut2(nt1dc(k))
            !!print*,"dut2(1:10)",dut2(1:10)
            res2_dcs(pourtrial(k)) =res2_dcs(pourtrial(k))+dlog(dut2(nt1dc(k)))+vet2
            !!print*,"res2_dcs(pourtrial(k))=",res2_dcs(pourtrial(k))
            !stop
            if ((res2_dcs(pourtrial(k)).ne.res2_dcs(pourtrial(k))).or.(abs(res2_dcs(pourtrial(k))).ge. 1.d30)) then
                funcpajsplines_surrogate_1=-1.d9
                goto 123
            end if    
        endif 
      
        ! scl pour calcul integrale
        !const_aux1(pourtrial(k))=const_aux1(pourtrial(k))+ut2(nt1dc(k))*dexp(vet2)
        !if ((const_aux1(pourtrial(k)).ne.const_aux1(pourtrial(k))).or.(abs(const_aux1(pourtrial(k))).ge. 1.d30)) then
        !    funcpajsplines=-1.d9
        !    goto 123
        !end if
        const_res5(g(k)) = const_res5(g(k))+ut2(nt1dc(k))* dexp(vet2) 
!        !print*,"je suis la bien funcpajsplines_surr ligne 260, k=",k,"const_res5(k)=",const_res5(k)
        if ((const_res5(g(k)).ne.const_res5(g(k))).or.(abs(const_res5(g(k))).ge. 1.d30)) then
            funcpajsplines_surrogate_1=-1.d9
            goto 123
        end if
    end do
    !!print*,"res2_dcs(k)=",sum(res2_dcs(:))
    !!print*,"const_res5=",const_res5
    !!print*,"const_res4=",const_res4
    !stop
!!print*,"size(const_res5)",size(const_res5)
!**************INTEGRALES ****************************
!!print*,"ng=",ng,"ntrials=",ntrials

    !================================================================================
    !==========distribution lognormale des effects aleatoires==============================
    !================================================================================
    
    if (logNormal==1) then 
        select case(methodInt)
        case(0) ! estimation par monte carlo
            posind_i=1
            !call cpu_time(c3)
            do ig=1,ntrials 
                !auxig=ig pas utiliser
                
                allocate(mu(nsujeti(ig))) !initialisation du vecteur des moyennes des effects alatoires
                allocate(vc(nsujeti(ig),nsujeti(ig))) !initialisation de la matrice de variance covariance pour le MC
                mu=0.d0
                vc=0.d0
                do l=1,nsujeti(ig)
                    vc(l,l)=theta2    
                end do        
                vcdiag=0 ! la matrice n'est plus diagonale
                resultatInt=0.d0
                !calcul de l'integrale par monte carlo pour l'integrale multiple et quadrature adaptative ou pas pour l'integrale su vsi et vti
                if(vectorisation.eq.1) then ! on vectorise, reduction du temps de calcul
                    call monteCarlosMult(integrant_indiv_1MCA,mu,vc,nsim,vcdiag,ig,resultatInt)
                else ! on ne vectorise pas
                    resultatInt=monteCarlos_ind(integrant_indiv_1MC,mu,vc,nsim,vcdiag)
                endif
                !!print*,"funcpa: resultatInt(ig)=",resultatInt(ig),"ig=",ig
                !stop
                
                posind_i=posind_i+nsujeti(ig)
                if(resultatInt(1).eq.0.d0) then
                    integrale3(ig)=0.1d-300
                    !print*,"integrale nulle et affectation de la valeur 0.1d-300, trials:",ig
                else
                    integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                end if
                !!print*,"funcpajsplines_surr ligne 312 nsujeti(ig)=",nsujeti(ig),"resultatInt=",integrale3(ig)  
                deallocate(mu,vc)
            end do
!            call cpu_time(c4)
!            !print*,"t1=",c4-c3
        !!print*,"cas Monte carlo"
        !!print*,integrale3
        case(1)! quadature classique (non-adaptative) ou pseudo-adaptative selon le contenu de la variable adaptative
!            !print*,"funcpajsplines_surrogate_1_scl.f90 lige 307: quadrature pas encore implémenté"
            posind_i=1
            l=1
            res = 0.d0
            do ig=1,ntrials 
                nsujet_trial=nsujeti(ig)
                !allocate(frail(nsujet_trial))
                !allocate(vc(nsujeti(ig),nsujeti(ig))) !initialisation de la matrice de variance covariance
                !allocate(vcinv(nsujeti(ig),nsujeti(ig)))
                !vc=0.d0
                !do l=1,nsujeti(ig)
                !    vc(l,l)=theta2    
                !    vcinv(l,l)=1/theta2 !inverse de la matrice des covariances
                !end do    
                !initialisation du vecteur frail avec xx1(1)
                !frail=xx1(1)
                !k = size(frail,dim=1)
                !!print*,"k=",k
                !!print*,"frail=",frail        
                resultatInt=0.d0
                if(vectorisation.eq.0) then ! on ne vectorise pas
                    resultatInt=gauss_HermMult(integrant_indiv_1,gauherJ_scl,npoint,nsujet_trial)
                else ! on vectorise, reduction du temps de calcul
                    resultatInt=gauss_HermMultA(integrant_indiv_1A,npoint,nsujet_trial)
                endif
                !!print*,"funcpa resultatInt(1)=",-(1/2.d0)*(nsujeti(ig)*dlog(2.d0*pi)+dlog(theta2**nsujeti(ig)))+resultatInt(1)
                !stop
                !resultatInt=gauss_HermMult(integrant_indiv_1,npoint,nsujet_trial)
                !!print*,resultatInt(1),npoint,posind_i,nsujeti(ig)
                posind_i=posind_i+nsujeti(ig)
                !!print*,"funcpajsplines_surr ligne 285: resultatInt=",resultatInt(1),"ni=",nsujeti(ig)  
                if(resultatInt(1).eq.0.d0) then
                    integrale3(ig)=0.1d-300
                    !print*,"integrale nulle et affectation de la valeur 0.1d-300, trials:",ig
                else
                    integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                endif 
                !deallocate(vc,frail,vcinv)
            end do
        !determinant=theta2**nsujeti(ig)
        !!print*,"cas quadrature"
        !!print*,-(1/2.d0)*(nsujeti(ig)*dlog(2.d0*pi)+dlog(determinant))+dlog(integrale3)
        !stop
        end select

    !!print*,"res=",res
    
    !************* FIN INTEGRALES **************************
                   
        res=0.d0
        select case(methodInt)
        case(0) ! estimation par monte carlo
        do k=1,ntrials!ng  
            if(cpt(k).gt.0)then
                if(sigma2.gt.(1.d-8)) then      
                    res= res + res2s(k) &
                    + res2_dcs(k)&
                    !+ integrale3(k)
                    + dlog(integrale3(k))
                else
                    res= res + res2s(k) &
                    + res2_dcs(k)  &
                    !+ integrale3(k)
                    + dlog(integrale3(k))
                endif
                
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpajsplines_surrogate_1=-1.d9
                    goto 123
                end if    
            endif
!            !print*,'k',k
            !!print*,'res',res,'dlog(integrale3(k))', dlog(integrale3(k))     
        end do
        case(1) !quadrature non adaptative
        !!print*,"sum(res2s)=",SUM(res2s(1:ntrials)),"size(res2s)=",size(res2s)
        !!print*,"sum(res2_dcs)=",SUM(res2_dcs(1:ntrials)),"size(res2_dcs)=",size(res2_dcs)
        do k=1,ntrials!ng  
            if(cpt(k).gt.0)then
                if(sigma2.gt.(1.d-8)) then  
                    determinant=theta2**nsujeti(k) ! produit des elements de la diagonale
                        res= res + res2s(k) &
                    + res2_dcs(k)&
                    -(1/2.d0)*(nsujeti(k)*dlog(2.d0*pi)+dlog(determinant))&
                    !+ integrale3(k)
                    + dlog(integrale3(k))
                else
                    determinant=theta2**nsujeti(k) ! produit des elements de la diagonale
                    res= res + res2s(k) &
                    + res2_dcs(k)&
                    -(1/2.d0)*(nsujeti(k)*dlog(2.d0*pi)+dlog(determinant))&
                    !+ integrale3(k)
                    + dlog(integrale3(k))
                endif
                
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpajsplines_surrogate_1=-1.d9
                    goto 123
                end if    
            endif
!            !print*,'k',k
            !!print*,'res',res,'dlog(integrale3(k))', dlog(integrale3(k))     
        end do
        end select
    
    !================================================================================
    !==========distribution gamma des effets aleatoires==============================
    !================================================================================
    
    else
        ! calcul integral
        posind_i=1
        l=1
        res = 0.d0
        !!print*,ntrials
        !stop
        do ig=1,ntrials 
            nsujet_trial=nsujeti(ig)
            resultatInt=0.d0
            if(vectorisation.eq.0) then ! on ne vectorise pas
                resultatInt=gauss_HermMult(integrant_indiv_1,gaulagJ_scl,npoint,nsujet_trial)
            else ! on vectorise, reduction du temps de calcul
                resultatInt=gauss_HermMultA(integrant_indiv_1A,npoint,nsujet_trial)
            endif
            posind_i=posind_i+nsujeti(ig)
            if(resultatInt(1).eq.0.d0) then
                integrale3(ig)=0.1d-300
                !print*,"integrale nulle et affectation de la valeur 0.1d-300, trials:",ig
            else
                integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
            endif 
        end do
    
    !************* FIN INTEGRALES **************************   
    
        res=0.d0
        pourgam=0.d0
        do k=1,ntrials!ng  
            if(cpt(k).gt.0)then
                if(theta2.gt.(1.d-8)) then  
                    pourgam=nsujeti(k)*(logGammaJ(1./theta)+(1./theta)*dlog(theta)) 
                    res= res + res2s(k)+ res2_dcs(k)-pourgam &
                    + integrale3(k)
                else
                    pourgam=nsujeti(k)*(logGammaJ(1./theta)+(1./theta)*dlog(theta)) 
                    res= res + res2s(k)+ res2_dcs(k)-pourgam &
                    + integrale3(k)
                endif
                
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpajsplines_surrogate_1=-1.d9
                    goto 123
                end if    
            endif
            !!print*,nsujeti(k),(logGammaJ(1./theta)+(1./theta)*dlog(theta))
            !!print*,'res',res,'dlog(integrale3(k))', dlog(integrale3(k))     
        end do

    endif
!!print*,res2s,res2_dcs,sum(res2s)+sum(res2_dcs),-pourgam,integrale3(k)
!!print*,"res =",res
!stop
    
    !================================================================================
    !=========================calcul de la penalisation==============================
    !================================================================================

    pe1 = 0.d0
    pe2 = 0.d0
    
    do i=1,n-3
        pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
        *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
        the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
        m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
        the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
        m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
        *the1(i)*m1m(i))
        if(nst.eq.1)then
            pe2=0.d0
        else
            pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
            *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
            the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
            m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
            the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
            m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
            *the2(i)*m1m(i))
        endif
    end do
    
    pe = k0(1)*pe1 + k0(2)*pe2 
    resnonpen = res
    res = res - pe
    !!print*,"funcpajsplines ligne 434, vraisemblance penalisee res=",res
    !cpteu=cpteu+1
    !if(cpteu.eq.1000) then
    !    stop
    !endif

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpajsplines_surrogate_1=-1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123

    else
        funcpajsplines_surrogate_1 = res 
        ! section encore a definir en fonction de la suite
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=cdc(k)
        end do
    end if
!Ad:
123     continue

    return
    
    end function funcpajsplines_surrogate_1



