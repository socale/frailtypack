    
!========================    FUNCPARES RESIDUS MATRINGALE DENSITE A POSTERIORI       ====================

!!!!
!!!! Calcul Residus shared log-normal
!!!!
    double precision function funcpasres(uu,np,id,thi,jd,thj)

    use comon
    use residusM

    implicit none

    integer,intent(in)::id,jd,np
    double precision,dimension(np)::bh
    double precision,dimension(np),intent(in)::uu
    double precision,intent(in)::thi,thj
    double precision::frail

    bh=uu
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    frail = bh(1)*bh(1)

    funcpasres = dexp(nig(indg)*frail - &
    dexp(frail)*cumulhaz(indg) - (frail**2.d0)/(2.d0*sig2))
    
    return

    end function funcpasres


!!!!
!!!! Calcul Residus joint
!!!!
! la differenciation dans le calcul entre joint cluster
! et classique se fait dans le funcpa

    double precision function funcpajres(uu,np,id,thi,jd,thj)

    use comon
        use residusM
    
    implicit none

    integer,intent(in)::id,jd,np
    double precision,intent(in)::thi,thj
    double precision,dimension(np)::uu,bh
    double precision::frail1,res
    double precision,parameter::pi=3.141592653589793d0

    bh=uu
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj
    
    frail1=bh(1)*bh(1)

!    funcpajres = frail1**(Ndc(indg) + Nrec(indg) + 1.d0/theta - 1.d0 + &
!    alpha * (Nrec(indg) + Ndc(indg))) * dexp(-frail1*(1.d0/theta + &
!    Rrec(indg))) * dexp(-(frail1**alpha)*Rdc(indg))

    res = frail1**(Nrec(indg) + 1.d0/theta - 1.d0 + &
    alpha * Ndc(indg)) * dexp(-frail1*(1.d0/theta + &
    Rrec(indg))) * dexp(-(frail1**alpha)*Rdc(indg))

!    res = dexp(-(frail1**alpha)*Rdc(indg)) * frail1**Nrec(indg) &
!    * frail1**(1.d0/theta - 1.d0) * frail1**(alpha * Ndc(indg)) &
!    * dexp(-frail1*(1.d0/theta + Rrec(indg)))

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then ! limite Ã  e+300, ce n'est pas une log-vraisemblance
        funcpajres=-1.d9
        goto 222
    end if

    funcpajres = res

222    continue

    return
    
    end function funcpajres

!!!!
!!!! Calcul Residus joint log-normal
!!!!

    double precision function funcpajres_log(uu,np,id,thi,jd,thj)

    use comon
        use residusM
    
    implicit none

    integer,intent(in)::id,jd,np
    double precision,intent(in)::thi,thj
    double precision,dimension(np)::uu,bh
    double precision::frail,res
    double precision,parameter::pi=3.141592653589793d0

    bh=uu
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj
    
    frail=bh(1)*bh(1)

    res = dexp((Nrec(indg)+alpha*Ndc(indg))*frail- &
    dexp(frail)*Rrec(indg) - dexp(alpha*frail)*Rdc(indg)- &
    (frail**2.d0)/(2.d0*sig2))

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpajres_log=-1.d9
        goto 222
    end if

    funcpajres_log = res

222    continue

    return
    
    end function funcpajres_log

!!!!
!!!! Calcul Residus nested
!!!!

    double precision function funcpanres(uu,np,id,thi,jd,thj)
    
    use comon,only:alpha,eta
        use residusM
    use commun

    
    implicit none

    integer,intent(in)::id,jd,np
    double precision,intent(in)::thi,thj    
    double precision,dimension(np),intent(in)::uu
    integer::j
    double precision,dimension(np)::bh
    double precision::frail1,prod1,prod2,prod3,res
    double precision,dimension(np-1)::frail2

    bh=uu
    
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj    
    
    frail1=bh(1)*bh(1)
    
    do j=1,n_ssgbygrp(indg)
        frail2(j)=bh(j+1)*bh(j+1)
    end do

    prod1 = 1.d0
    prod2 = 1.d0
    prod3 = 1.d0
    
    do j=1,n_ssgbygrp(indg)
        prod1 = prod1 * (frail2(j)**mij(indg,j)) * dexp(-frail1 * frail2(j) * cumulhaz1(indg,j))
        prod2 = prod2 * frail2(j)**((1.d0/eta) - 1.d0) * dexp(-frail2(j)/eta)
        prod3 = prod3 * dexp(-frail1 * frail2(j) * cumulhaz0(indg,j))
    end do

    res = frail1**(mid(indg)+1.d0/alpha - 1.d0) * prod1 * prod3 * dexp(-frail1/alpha) * prod2
    
    if ((res.ne.res).or.(abs(res).ge. 1.d300)) then
        funcpanres=-1.d9
        goto 333
    end if

    funcpanres = res

333    continue

    return
    
    end function funcpanres


        
!!!!
!!!! Calcul Residus joint nested (AK 12/12/2016)
!!!!
    
    double precision function funcpajres_fam(uu,np,id,thi,jd,thj)
    
    use comon,only:alpha,eta,xi,theta,cdc,fsize
    use residusM
    use commun
    
    implicit none

    integer,intent(in)::id,jd,np
    double precision,intent(in)::thi,thj    
    double precision,dimension(np),intent(in)::uu
    integer::j
    double precision,dimension(np)::bh
    double precision::frail1,prod1,prod2,prod3,res,prod4,prod5
    double precision,dimension(np-1)::frail2

    bh=uu
   
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj    
    
    frail1=bh(1)*bh(1)
    
    do j=1,fsize(indg)
        frail2(j)=bh(j+1)*bh(j+1)
    end do

    prod1 = 1.d0
    prod2 = 1.d0
    prod3 = 1.d0
    prod4 = 1.d0
    prod5 = 1.d0        

    
    do j=1,fsize(indg)
        prod1 = prod1 * (frail2(j)**Nrec_ind(sum(fsize(1:(indg-1)))+j)) * dexp(-frail1**xi * frail2(j) * cumulhaz1(indg,j))
        prod2 = prod2 *  dexp(-frail2(j)/theta)
        prod3 = prod3 * dexp(-frail1**xi * frail2(j) * cumulhaz0(indg,j))
        prod4 = prod4 * dexp(-frail1 * frail2(j)**alpha * cumulhazdc(indg,j))
        prod5 = prod5 * frail2(j)**(Nrec_ind(sum(fsize(1:(indg-1)))+j)+alpha*cdc(sum(fsize(1:(indg-1)))+j))
    end do

    res = prod1 * prod3 * dexp(-frail1/eta) * prod2 * &! frail1**(1.d0/eta - 1.d0) *
            prod4 * prod5 * frail1**(Nrec_fam(indg)*xi+Ndc_fam(indg))
            
    if ((res.ne.res).or.(abs(res).ge. 1.d300)) then
        funcpajres_fam=-1.d9
        goto 333
    end if
   funcpajres_fam = res

333    continue
 
    return    
    end function funcpajres_fam  
        
    
!!!!
!!!! Calcul Residus additive
!!!!

    double precision function funcpaares(uu,np,id,thi,jd,thj)
    
    use comon
        use residusM
    !use additiv,only:mid,Xbeta
    use additiv,only:ve2,ut1,ut2    

    implicit none

    integer,intent(in)::id,jd,np
    double precision,intent(in)::thi,thj
    double precision,dimension(np),intent(in)::uu
    integer::k,ip
    double precision,dimension(np)::bh
    double precision::frail1,frail2,som1,som2
    double precision,parameter::pi=3.141592653589793d0
    double precision::result
    double precision,dimension(1,2)::apres
    double precision,dimension(2,1)::avant    
    double precision,dimension(1,1)::res    
    double precision::vet


    bh=uu
!    kapa=k0(1)
        
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj    
    
    frail1=bh(1)
    frail2=bh(2)
    
    apres(1,1)=frail1
    apres(1,2)=frail2    
    
    avant(1,1)=frail1
    avant(2,1)=frail2
    
    res = (-1.d0/2) * matmul(apres,matmul(invsigma,avant))    
    result = res(1,1)

    som1 = 0.d0
    do k=1,nsujet
        if(g(k) == indg)then
            if(c(k).eq.1)then
                som1 = som1 + frail1 + frail2 * ve2(k,1) + dlog(som_Xbeta(indg))
            end if
        end if
    end do
    
    som2 = 0.d0

    do k=1,nsujet
        if(nva.gt.0 .and.g(k).eq.indg)then
            vet = 0.d0 
            do ip=1,nva
                vet = vet + b_temp(np-nva+ip)*ve(k,ip)
            end do
            vet = dexp(vet)
        else
            vet=1.d0
        endif
        if(typeof==0) then
            if(g(k) == indg)then    
                if(stra(k).eq.1)then
                    som2 = som2  - 1.d0 * ut1(nt1(k)) * dexp(frail1 + frail2 * ve2(k,1) + dlog(vet))
                end if
                if(stra(k).eq.2)then
                    som2 = som2 - 1.d0 * ut2(nt1(k)) * dexp(frail1 + frail2 * ve2(k,1) + dlog(vet))
                end if
            end if
        else
            if(g(k) == indg)then    
                som2 = som2  - 1.d0 * cumulhaz(g(k)) * dexp(frail1 + frail2 * ve2(k,1) + dlog(vet))
            end if

        end if
    end do

    funcpaares = 1.d0/(2.d0 * pi * dsqrt(detSigma))* dexp(som1 + som2 + result)
    
    return
    
    end function funcpaares    




!========================    FUNCPARES MULTIVE RESIDUS MATRINGALE DENSITE A POSTERIORI       ====================


    double precision function funcpamultires(uu,np,id,thi,jd,thj)
    
    use comonmultiv,only:eta,theta,alpha,alpha1,alpha2
    use residusMmultiv
    
    IMPLICIT NONE
    
    integer,intent(in)::id,jd,np
    double precision,dimension(np)::bh
    double precision,dimension(np),intent(in)::uu
!    double precision,dimension(3),intent(in)::k0
    double precision,intent(in)::thi,thj
    double precision::frail1,frail2
    double precision,parameter::pi=3.141592653589793d0
!    double precision,dimension(3)::kappa_tmp

!    kappa_tmp = k0
    bh=uu

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj    

    
    frail1=bh(1)
    frail2=bh(2) 


!---------- calcul de la penalisation -------------------

    funcpamultires= frail1*(Ndc(indg)*alpha1+Nrec(indg)) &
    +frail2*(Ndc(indg)*alpha2+Nrec2(indg)) &
    -dexp(frail1)*Rrec(indg)-dexp(frail2)*Rrec2(indg) &
    -dexp(frail1*alpha1+frail2*alpha2)*Rdc(indg) &
    +(2.d0*((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0) &
    *frail1*frail2/sqrt(theta*eta)  &
    -(frail1**2.d0)/theta -(frail2**2.d0)/eta) &
    /(2.d0*(1.d0-((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0)**2.d0)) 

    return
    
    end function funcpamultires


    
     


!!!!
!!!! =============================== Calcul Residus joint bivariate - longitudinales et deces ==============================
    
    
    
        double precision function funcpajres_biv(uu,np,id,thi,jd,thj)
    
            use optim
        use comon
            use donnees_indiv,only:z1cur,x2cur,current_mean,b1
            use residusM
    
        implicit none
    
        integer,intent(in)::id,jd,np
        double precision,intent(in)::thi,thj
        double precision,dimension(np)::uu,bh
        double precision::res,yscalar,prod_cag,eps
            double precision::finddet,alnorm
            double precision,dimension(nb1)::b_vec,uii
            double precision,dimension(nb1*(nb1+1)/2)::matv
            double precision,dimension(nb1,1)::b_vecT
            double precision,dimension(nb1,nb1)::mat_B
            double precision,dimension(1)::uiiui
            integer::jj,j,k,ier
            logical :: upper
        double precision,parameter::pi=3.141592653589793d0
            double precision :: resultdc
    double precision,external::survdcCM
        double precision :: abserr,resabs,resasc
            upper=.false.
        bh=uu
    
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
    
        b_vec(1:nb1) = bh(1:nb1)
        b_vecT(1:nb1,1) = bh(1:nb1)
          
          matv = 0.d0
            mu1_res = 0.d0
    
            if(nmesy(indg).gt.0) then
            !if(nb1.eq.1)mu1_res(1:nmesy(indg)) = XbetaY_res(1,indg) +Zet(1:nmesy(indg),1:netadc)*b_vec
            mu1_res(1:nmesy(indg)) = XbetaY_res(1,it_res:(it_res+nmesy(indg)-1)) &
                    +MATMUL(Zet(it_res:(it_res+nmesy(indg)-1),1:nb1),b_vec)
            end if
    
            !********* Left-censoring ***********
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 0,nmesy(indg)-1
                            if(yy(k+it_res).le.s_cag) then
                                    prod_cag = prod_cag*(1.d0-alnorm((mu1_res(k+1)-s_cag)/sqrt(sigmae),upper))
                            else
                                    yscalar = yscalar + (yy(k+it_res)-mu1_res(k+1))**2
                            end if
                    end do
            else
                    do k=0,nmesy(indg)-1
                            yscalar = yscalar + (yy(k+it_res)-mu1_res(k+1))**2.d0
                    end do
            end if
    
            yscalar = dsqrt(yscalar)
    
    
            mat_B = matmul(ut,utt)
            det = finddet(matmul(ut,utt),nb1)
    
        if(nb1.ge.2) then
                    jj=0
                    do j=1,nb1
                            do k=j,nb1
                                    jj=j+k*(k-1)/2
                                    matv(jj)=mat_B(j,k)
                            end do
                    end do
                    ier = 0
                    eps = 1.d-10
                    call dsinvj(matv,nb1,eps,ier)
                    mat_b=0.d0
                    do j=1,nb1
                            do k=1,nb1
                                    if (k.ge.j) then
                                            mat_b(j,k)=matv(j+k*(k-1)/2)
                                    else
                                            mat_b(j,k)=matv(k+j*(j-1)/2)
                                    end if
                            end do
                    end do
            else
    
                    matv(1) = 1.d0/mat_B(1,1)
    
                    Mat_B(1,1) = matv(1)
            end if
    
            if(link.eq.2) then
    
            call integrationdc(survdcCM,t0dc(indg),t1dc(indg),resultdc,abserr,resabs,resasc,indg,b1,npp,b_vec)
    
            X2cur = 0.d0
    
                    X2cur(1,1) = 1.d0
            X2cur(1,2) =t1dc(indg)
            if((nva3-2).gt.0) then
                do k=3,nva3
                        X2cur(1,k) = dble(vey(it_res,k))
                    end do
            end if
    
            Z1cur(1,1) = 1.d0
                    if(nb1.eq.2)  Z1cur(1,2) = t1dc(indg)
            current_mean = 0.d0
    
            current_mean(1) =dot_product(X2cur(1,1:nva3),b1((npp-nva3+1):npp))&
                                    +dot_product(Z1cur(1,1:nb1),b_vec(1:nb1))
            end if
    
            uii = matmul(b_vec,mat_b)
            uiiui=matmul(uii,b_vecT)
    
            if(prod_cag.lt.0.1d-321)prod_cag= 0.1d-321
    
            res = 0.d0
            if(link.eq.1) then
                res = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+&
                Ndc(indg)*dot_product(etaydc,b_vec) - &
                Rdc(indg)*dexp(dot_product(etaydc,b_vec))-uiiui(1)/2.d0
          
            else
                res = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                                    +Ndc(indg)*(etaydc(1)*current_mean(1))  &
                                    -uiiui(1)/2.d0 - resultdc!      -       Rdc(indg)*dexp(etaydc1*current_mean(1))
            
            end if
    
        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
            funcpajres_biv=-1.d9
            goto 222
        end if
            funcpajres_biv = res
    
    222    continue
    
        return
    
        end function funcpajres_biv
    
    
    
    
    !!!!
    !!!! =============================== Calcul Residus joint bivariate - longitudinales et deces ==========
    
    
    
        double precision function funcpajres_tri(uu,np,id,thi,jd,thj)
    
            use optim
        use comon
            use donnees_indiv,only:z1cur,x2cur,current_mean,b1
            use residusM
    
        implicit none
    
        integer,intent(in)::id,jd,np
        double precision,intent(in)::thi,thj
        double precision,dimension(np)::uu,bh
        double precision::res,yscalar,prod_cag,eps !frail1,frail2,frail3
            double precision::finddet,alnorm
            double precision,dimension(nea)::b_vec,uii
            double precision,dimension(nea*(nea+1)/2)::matv
            double precision,dimension(nea,1)::b_vecT
            double precision,dimension(nea,nea)::mat_B
            double precision,dimension(1)::uiiui
            integer::jj,j,k,ier,ii
            logical :: upper
        double precision,parameter::pi=3.141592653589793d0
            double precision :: resultdc,resultr
            double precision,external::survdcCM,survRCM
        double precision :: abserr,resabs,resasc,ress
            double precision,dimension(1)::current_meanR
            upper=.false.
            
        bh=uu
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
    
        b_vec(1:nea) = bh(1:nea)
        b_vecT(1:nea,1) = bh(1:nea)
         
            matv = 0.d0
            mu1_res = 0.d0
    
            if(nmesy(indg).gt.0) then
            if(nb1.eq.2)mu1_res(1:nmesy(indg)) = XbetaY_res(1,indg) +MATMUL(Zet(1:nmesy(indg),1:nb1),b_vec(1:nb1))
            end if
    
            !********* Left-censoring ***********
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 0,nmesy(indg)-1
                            if(yy(k+it_res).le.s_cag) then
                                    prod_cag = prod_cag*(1.d0-alnorm((mu1_res(k+1)-s_cag)/sqrt(sigmae),upper))
                                    !0.5d0*(1.d0-erf((mu1_res(k+1)-s_cag)/(sigmae*dsqrt(2.d0))))        
                            else
                                    yscalar = yscalar + (yy(k+it_res)-mu1_res(k+1))**2
                            end if
                    end do
            else
                    do k=0,nmesy(indg)-1
                            yscalar = yscalar + (yy(k+it_res)-mu1_res(k+1))**2
                    end do
            end if
    
            yscalar = dsqrt(yscalar)
    
    
            mat_B = matmul(ut,utt)
            det = finddet(matmul(ut,utt),nea)
    
    
                    jj=0
                    do j=1,nea
                            do k=j,nea
                                    jj=j+k*(k-1)/2
                                    matv(jj)=mat_B(j,k)
                            end do
                    end do
                    ier = 0
                    eps = 1.d-10
                    call dsinvj(matv,nea,eps,ier)
                    mat_b=0.d0
                    do j=1,nea
                            do k=1,nea
                                    if (k.ge.j) then
                                            mat_b(j,k)=matv(j+k*(k-1)/2)
                                    else
                                            mat_b(j,k)=matv(k+j*(j-1)/2)
                                    end if
                            end do
                    end do
    
        ress = 0.d0
        current_mean = 0.d0
        current_meanR = 0.d0
            if(link.eq.2) then
                    
                
                    do ii=it_res_rec,it_res_rec+nmesrec(indg)-1
                            resultR = 0.d0
    
                            call integrationdc(survRCM,t0(ii),t1(ii),resultR,abserr,resabs,resasc,ii,b1,npp,b_vec)
                            ress = ress + resultR   !c'est deja res1-res3
    
    
                            if((c(ii).eq.1))then
                                    X2cur(1,1) = 1.d0
                                    X2cur(1,2) = t1(ii)
                                    if(nva3.gt.0) then
                                            do k=3,nva3
                                                    X2cur(1,k) = dble(vey(it_res,k))
                                            end do
                                    end if
    
                                    Z1cur(1,1) = 1.d0
                                    if(nb1.eq.2)Z1cur(1,2) = t1(ii)
                                    current_meanR = current_meanR + MATMUL(X2cur,b1((npp-nva3+1):npp))&
                                            +MATMUL(Z1cur,b_vec(1:nb1))
                            end if
                    end do
    
    
            call integrationdc(survdcCM,t0dc(indg),t1dc(indg),resultdc,abserr,resabs,resasc,indg,b1,npp,b_vec)
    
            X2cur = 0.d0
    
                    X2cur(1,1) = 1.d0
            X2cur(1,2) =t1dc(indg)
            if((nva3-2).gt.0) then
                do k=3,nva3
                        X2cur(1,k) = dble(vey(it_res,k))
                    end do
            end if
    
            Z1cur(1,1) = 1.d0
                    if(nb1.eq.2)  Z1cur(1,2) = t1dc(indg)
            
    
            current_mean(1) =dot_product(X2cur(1,1:nva3),b1((npp-nva3+1):npp))&
                                    +dot_product(Z1cur(1,1:nb1),b_vec(1:nb1))
            end if
    
    
            uii = matmul(b_vec,mat_b)
            uiiui=matmul(uii,b_vecT)
    
            if(prod_cag.lt.0.1d-321)prod_cag= 0.1d-321
            
            res = 0.d0
            if(link.eq.1) then
                 res = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+&
                       Ndc(indg)*(dot_product(etaydc,b_vec(1:nb1))+b_vec(nb1+1)*alpha) - &
                       Rdc(indg)*dexp(dot_product(etaydc,b_vec(1:nb1))+b_vec(1+nb1)*alpha)-uiiui(1)/2.d0 &
                        +Nrec(indg)*(b_vec(nb1+1)+dot_product(etayr,b_vec(1:nb1)))-&
                       Rrec(indg)*dexp(b_vec(nb1+1)+dot_product(etayr,b_vec(1:nb1)))
         
            else
                  
                res = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+&
                      Ndc(indg)*(etaydc(1)*current_mean(1)+b_vec(nb1+1)*alpha) - &
                      Rdc(indg)*dexp(etaydc(1)*current_mean(1)+b_vec(nb1+1)*alpha)-uiiui(1)/2.d0 &
                      +Nrec(indg)*(b_vec(nb1+1)+etayr(1)*current_meanR(1))&
                      -ress*dexp(b_vec(nb1+1)+etayr(1)*current_meanR(1))
              
            end if
    
        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
            funcpajres_tri=-1.d9
            goto 222
        end if
    
        funcpajres_tri = res
    
    222    continue
    
        return
    
        end function funcpajres_tri
 
 
 
!!!!
!!!! =============================== Calcul Residus joint bivariate - longitudinales et deces ==============================
    
    
    
        double precision function funcpares_uni(uu,np,id,thi,jd,thj)
    
            use optim
        use comon
            use donnees_indiv,only:b1,mu !z1cur,x2cur,current_mean
            use residusM
 !   use ParametresPourParallelisation
    
        implicit none
    
        integer,intent(in)::id,jd,np
        double precision,intent(in)::thi,thj
        double precision,dimension(np)::uu,bh
        double precision::frail1,frail2,frail3,frail4,res,yscalar,prod_cag,eps
            double precision::finddet,alnorm
            double precision,dimension(nb1)::b_vec,uii
            double precision,dimension(nb1*(nb1+1)/2)::matv
            double precision,dimension(nb1,1)::b_vecT
            double precision,dimension(nb1,nb1)::mat_B
            double precision,dimension(1)::uiiui
            integer::jj,j,k,ier
            logical :: upper
        double precision,parameter::pi=3.141592653589793d0
        !double precision :: abserr,resabs,resasc
            upper=.false.
        bh=uu
    
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
    
        frail1=bh(1)!*bh(1)
        frail2 = 0.d0
        frail3 = 0.d0
        
        b_vec(1) = frail1
        b_vecT(1,1) = frail1
        
         if(nb1.ge.2) then
            frail2 = bh(2)!*bh(2)
            b_vec(2) = frail2
            b_vecT(2,1) = frail2
        end if
        
        if(nb1.ge.3) then
            frail3 = bh(3)!*bh(2)
            b_vec(3) = frail3
            b_vecT(3,1) = frail3
        end if
    
        if(nb1.ge.4) then
            frail4 = bh(4)!*bh(2)
            b_vec(4) = frail4
            b_vecT(4,1) = frail4
        end if
    
    Ut = 0.d0
            Utt = 0.d0
    
                    do j=1,nb1
                do k=1,j
                    if(j.eq.k) then
    
                       Ut(j,k)=b1(npp-nva-nb_re+k)!sqrt(bh(np-nva-nb_re+k)**2.d0)
                  Utt(k,j)=b1(npp-nva-nb_re+k)!sqrt(bh(np-nva-nb_re+k)**2.d0)
                  
            end if
                end do
            end do
  

  
            matv = 0.d0
            mu1_res = 0.d0
    
            if(nmesy(indg).gt.0) then
 

    !    mu1_res(1:nmesy(indg)) =dexp(b1(npp-nva-nb_re-1)+b_vec(1)) * &
    !    dexp(dexp( b1(npp-nva-nb_re-1-3)+b_vec(2)+mu(1:nmesy(indg),1))*&
    !    ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
    !    dexp(b1(npp-nva-nb_re-1-2)+b_vec(3)-&
    !    (b1(npp-nva-nb_re-1-1))+&!+b_vec(2)
    !    mu(1:nmesy(indg),2))*&!*ziy((it+1):(it+nmescur),3)+b_vec(3)
    !    (dexp(-dexp(b1(npp-nva-nb_re-1-1))*&!+b_vec(3)
    !    ziy(it_res:(it_res+nmesy(indg)-1),1)) -1))!)**0.4d0-1)/0.4d0
    

         select case(which_random)
        case(1)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+b_vec(1) +dexp( b1(npp-nva-nb_re-1-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)-&
        (b1(npp-nva-nb_re-1-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(2)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+dexp( b1(npp-nva-nb_re-1-3)+b_vec(1)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)-&
        (b1(npp-nva-nb_re-1-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(3)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+dexp( b1(npp-nva-nb_re-1-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)+b_vec(1)-&
        (b1(npp-nva-nb_re-1-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))*box_cox_par-box_cox1)/box_cox_par
            
        case(4)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+dexp( b1(npp-nva-nb_re-1-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)-&
        (b1(npp-nva-nb_re-1-1)+b_vec(1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1)+b_vec(1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(5)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+b_vec(1) +dexp( b1(npp-nva-nb_re-1-3)+b_vec(2)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)-&
        (b1(npp-nva-nb_re-1-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(6)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+b_vec(1) +dexp( b1(npp-nva-nb_re-1-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)+b_vec(2)-&
        (b1(npp-nva-nb_re-1-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(7)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+b_vec(1) +dexp( b1(npp-nva-nb_re-1-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)-&
        (b1(npp-nva-nb_re-1-1)+b_vec(2))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1)+b_vec(2))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(8)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1) +dexp( b1(npp-nva-nb_re-1-3)+b_vec(1)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)+b_vec(2)-&
        (b1(npp-nva-nb_re-1-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(9)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1) +dexp( b1(npp-nva-nb_re-1-3)+b_vec(1)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)-&
        (b1(npp-nva-nb_re-1-1)+b_vec(2))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1)+b_vec(2))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(10)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1) +dexp( b1(npp-nva-nb_re-1-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)+b_vec(1)-&
        (b1(npp-nva-nb_re-1-1)+b_vec(2))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1)+b_vec(2))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(11)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+b_vec(1) +dexp( b1(npp-nva-nb_re-1-3)+b_vec(2)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)+b_vec(3)-&
        (b1(npp-nva-nb_re-1-1)+b_vec(4))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1)+b_vec(4))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(12)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+b_vec(1) +dexp( b1(npp-nva-nb_re-1-3)+b_vec(2)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)-&
        (b1(npp-nva-nb_re-1-1)+b_vec(3))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1)+b_vec(3))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(13)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+b_vec(1) +dexp( b1(npp-nva-nb_re-1-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)+b_vec(2)-&
        (b1(npp-nva-nb_re-1-1)+b_vec(3))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1)+b_vec(3))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(14)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1) +dexp( b1(npp-nva-nb_re-1-3)+b_vec(1)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)+b_vec(2)-&
        (b1(npp-nva-nb_re-1-1)+b_vec(3))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1)+b_vec(3))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(15)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1)+b_vec(1) +dexp( b1(npp-nva-nb_re-1-3)+b_vec(2)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-2)+b_vec(3)-&
        (b1(npp-nva-nb_re-1-1)+b_vec(4))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-1)+b_vec(4))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        end select
        
           end if
    
    
!    write(*,*)'yy',yy(it_res:it_res+nmesy(indg)-1)
    
    
!if(indg.eq.400)write(*,*)b1,mu(2,1),b_vec(1)
    
            !********* Left-censoring ***********
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 0,nmesy(indg)-1
                            if(yy(k+it_res).le.s_cag) then
                            
                                    prod_cag = prod_cag*(1.d0-alnorm((mu1_res(k+1)-s_cag)/sqrt(sigmae),upper))
                            else
                                    yscalar = yscalar + (yy(k+it_res)-mu1_res(k+1))**2
                            end if
                    end do
            else
                    do k=0,nmesy(indg)-1
                            yscalar = yscalar + (yy(k+it_res)-mu1_res(k+1))**2.d0
                    end do
            end if
    
            yscalar = dsqrt(yscalar)
        !    write(*,*)mu1_res(1),yy(1)
    
            mat_B = matmul(ut,utt)
            det = finddet(matmul(ut,utt),nb1)
    
        if(nb1.ge.2) then
                    jj=0
                    do j=1,nb1
                            do k=j,nb1
                                    jj=j+k*(k-1)/2
                                    matv(jj)=mat_B(j,k)
                            end do
                    end do
                    ier = 0
                    eps = 1.d-10
                    call dsinvj(matv,nb1,eps,ier)
                    mat_b=0.d0
                    do j=1,nb1
                            do k=1,nb1
                                    if (k.ge.j) then
                                            mat_b(j,k)=matv(j+k*(k-1)/2)
                                    else
                                            mat_b(j,k)=matv(k+j*(j-1)/2)
                                    end if
                            end do
                    end do
            else
    
                    matv(1) = 1.d0/mat_B(1,1)
    
                    Mat_B(1,1) = matv(1)
            end if
    
          

            uii = matmul(b_vec,mat_b)
            uiiui=matmul(uii,b_vecT)
    
            if(prod_cag.lt.0.1d-321)prod_cag= 0.1d-321
    
            res = 0.d0
        
                    if(nb1.eq.1) then
                            res =  dlog(prod_cag) &
                                            -(yscalar**2.d0)/(2.d0*b1(npp-nva-nb_re)*b1(npp-nva-nb_re))&
                                            -(frail1**2.d0)/(2.d0*ut(1,1)**2) 
    
                    else if(nb1.ge.2) then
                                    res = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*b1(npp-nva-nb_re)*b1(npp-nva-nb_re))&
                                   -uiiui(1)/2.d0
                
                    end if
          

    !            if(CurrentProcessorID.eq.0.and.indg.eq.1) write(*,*)mu1_res(1:2),'y0',b1(npp-nva-nb_re-1),&
    !    'k_g0',b1(npp-nva-nb_re-1-3),'first cova',mu(1:2,1),&
    !    'temps',ziy(it_res:(it_res+2),1), 'dose',ziy(it_res:(it_res+2),2),&
    !    'k_d0',b1(npp-nva-nb_re-1-2),&
    !    'lambda',b1(npp-nva-nb_re-1-1),&!+b_vec(2)
    !    'second cova',mu(1:2,2)
        
        !  if(indg.eq.50)write(*,*)'vu', res,yy(it_res),mu1_res(1),CurrentProcessorID
        !  stop
!      write(*,*)'res',res,b_vec
          
!    if(indg.eq.1)write(*,*)yscalar,b_vec,yy(it_res:it_res+nmesy(indg)-1)
    !write(*,*)dlog(prod_cag) , -(yscalar**2.d0),(2.d0*b1(npp-nva-nb_re)*b1(npp-nva-nb_re)),&
     !                                       -(frail1**2.d0)/(2.d0*ut(1,1)**2) 

    
        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
           funcpares_uni=-1.d9
            goto 222
        end if
            funcpares_uni = res
    
    222    continue
    
        return
    
        end function funcpares_uni
    

    
        
    !!!!
    !!!! =============================== Calcul Residus joint trivariate nonlinear ==========
    
    
    
        double precision function funcpajres_triNL(uu,np,id,thi,jd,thj)
    
            use optim
        use comon
            use donnees_indiv,only:current_mean,b1,mu,sigmav !z1cur,x2cur
            use residusM
        implicit none
    
        integer,intent(in)::id,jd,np
        double precision,intent(in)::thi,thj
        double precision,dimension(np)::uu,bh
        double precision::res,yscalar,prod_cag,eps !frail1,frail2,frail3,frail4
            double precision::finddet,alnorm
            double precision,dimension(nea)::b_vec,uii
            double precision,dimension(nea*(nea+1)/2)::matv
            double precision,dimension(nea,1)::b_vecT
            double precision,dimension(nea,nea)::mat_B
            double precision,dimension(1)::uiiui
            integer::jj,j,k,ier !ii
            logical :: upper
        double precision,parameter::pi=3.141592653589793d0
            !double precision :: resultdc,resultr
            double precision,external::survdcCM,survRCM
        double precision :: ress !abserr,resabs,resasc
            double precision,dimension(1)::current_meanR
            upper=.false.
            
        bh=uu
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
    
    
        b_vec(1:nea) = bh(1:nea)
        b_vecT(1:nea,1) = bh(1:nea)
            
        Ut = 0.d0
           Utt = 0.d0
    
                   do j=1,nb1
               do k=1,j
                   if(j.eq.k) then
    
                Ut(j,k)=b1(npp-nva-nb_re+k)
                 Utt(k,j)=b1(npp-nva-nb_re+k)
                  
           end if
               end do
           end do
        
           Ut(nea,nea) = sqrt(sigmav**2.d0)
                Utt(nea,nea) = sqrt(sigmav**2.d0)
        
        etaydc(1:netadc) = b1(npp-nva-nb_re-netadc:npp-nva-nb_re-1)
                 etayr(1:netar) =b1(npp-nva-nb_re-netadc - netar:npp-nva-nb_re-netadc - 1)
        

            matv = 0.d0
            mu1_res = 0.d0
    
            if(nmesy(indg).gt.0) then

             select case(which_random)
        case(1)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+b_vec(1) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
            
        case(2)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+b_vec(1)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
    
        case(3)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)+b_vec(1)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))*box_cox_par-box_cox1)/box_cox_par
            
        case(4)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(5)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+b_vec(1) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+b_vec(2)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(6)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+b_vec(1) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)+b_vec(2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(7)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+b_vec(1) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(2))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(2))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(8)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+b_vec(1)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)+b_vec(2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(9)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+b_vec(1)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(2))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(2))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(10)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)+b_vec(1)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(2))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(2))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(11)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+b_vec(1) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+b_vec(2)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)+b_vec(3)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(4))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(4))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(12)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+b_vec(1) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+b_vec(2)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(3))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(3))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(13)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+b_vec(1) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)+b_vec(2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(3))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(3))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(14)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+b_vec(1)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)+b_vec(2)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(3))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(3))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        case(15)
                mu1_res(1:nmesy(indg)) =(( &
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha)+b_vec(1) +&
        dexp( b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)+b_vec(2)+mu(1:nmesy(indg),1))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1) + ziy(it_res:(it_res+nmesy(indg)-1),2)*&
        dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)+b_vec(3)-&
        (b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(4))+&
        mu(1:nmesy(indg),2))*&
        (dexp(-dexp(b1(npp-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)+b_vec(4))*&
        ziy(it_res:(it_res+nmesy(indg)-1),1)) -1)))**box_cox_par-box_cox1)/box_cox_par
            
        end select
        
        
        
        
        end if
    

    
            !********* Left-censoring ***********
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 0,nmesy(indg)-1
                            if(yy(k+it_res).le.s_cag) then
                                    prod_cag = prod_cag*(1.d0-alnorm((mu1_res(k+1)-s_cag)/sqrt(sigmae),upper))
                                           
                            else
                                    yscalar = yscalar + (yy(k+it_res)-mu1_res(k+1))**2
                            end if
                    end do
            else
                    do k=0,nmesy(indg)-1
                            yscalar = yscalar + (yy(k+it_res)-mu1_res(k+1))**2
                    end do
            end if
    
            yscalar = dsqrt(yscalar)
    
    
            mat_B = matmul(ut,utt)
            det = finddet(matmul(ut,utt),nea)
    
    
                    jj=0
                    do j=1,nea
                            do k=j,nea
                                    jj=j+k*(k-1)/2
                                    matv(jj)=mat_B(j,k)
                            end do
                    end do
                    ier = 0
                    eps = 1.d-10
                    call dsinvj(matv,nea,eps,ier)
                    mat_b=0.d0
                    do j=1,nea
                            do k=1,nea
                                    if (k.ge.j) then
                                            mat_b(j,k)=matv(j+k*(k-1)/2)
                                    else
                                            mat_b(j,k)=matv(k+j*(j-1)/2)
                                    end if
                            end do
                    end do
    
        ress = 0.d0
        current_mean = 0.d0
        current_meanR = 0.d0
           
    
    
            uii = matmul(b_vec,mat_b)
            uiiui=matmul(uii,b_vecT)
    
            if(prod_cag.lt.0.1d-321)prod_cag= 0.1d-321
      
            res = 0.d0
          

                            res = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+&
                                    Ndc(indg)*(dot_product(etaydc,b_vec(1:netadc))+b_vec(nea)*alpha) - &
                                    Rdc(indg)*dexp(dot_product(etaydc,b_vec(1:netadc))+b_vec(nea)*alpha)-uiiui(1)/2.d0 &
                                    +Nrec(indg)*(b_vec(nea)+dot_product(etayr,b_vec(1:netar)))-&
                                    Rrec(indg)*dexp(b_vec(nea)+dot_product(etayr,b_vec(1:netar)))
                        
        

    
    
        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
            funcpajres_triNL=-1.d9
            goto 222
        end if
    
        funcpajres_triNL = res
    
    222    continue
    
        return
    
        end function funcpajres_triNL
 