
! ! Score Test Commenges/Andersen for shared frailty model only
! !=============================================================================
!     double precision function scoretest(b,frailtypred)
! 
!     use residusM,only:cumulhaz
!     use tailles,only:npmax
!     use comon,only:nsujet,ng,nva,ve,t1,g,c
! 
!     implicit none
!     
!     double precision,dimension(npmax),intent(in)::b
!     double precision,dimension(ng),intent(in)::frailtypred
!     integer::i,j,N
!     double precision,dimension(ng)::Mij
!     double precision::vet,somme,integrale
! 
!     Mij = 0.d0
!     do i=1,nsujet
!         if(nva.gt.0)then
!             vet = 0.d0
!             do j=1,nva
!                 vet =vet + b(npmax-nva+j)*dble(ve(i,j))
!             end do
!             vet = dexp(vet)
!         else
!            vet=1.d0
!         endif
!         Mij(g(i)) = Mij(g(i)) + c(i) !- ut1(nt1(i))*vet
!     end do
! print*,Mij
!     somme = 0.d0
!     do i=1,ng
!         Mij(i) = Mij(i) - cumulhaz(i)
!         somme = somme + Mij(i)*Mij(i)
!     end do
! 
!     N = sum(c)
! 
!     call gaulagSC(integrale,b,frailtypred)
! print*,somme,N,integrale
!     scoretest = somme - N + integrale
! 
!     return
! 
!     end function scoretest
! 
! 
! ! gauss laguerre (integrale sur 0 , +infty)
!     subroutine gaulagSC(ss,b,frailtypred)
! 
!     use tailles,only:npmax
!     use comon,only:ng
!     use donnees,only:w,x
!     
!     implicit none
! 
!     double precision,intent(out)::ss
!     double precision,dimension(npmax),intent(in)::b
!     double precision,dimension(ng),intent(in)::frailtypred
! 
!     integer::j
!     double precision::auxfunca,func
!     external::func
! 
!     ss = 0.d0
!     do j=1,20
!         auxfunca = func(x(j),b,frailtypred)
!         ss = ss+w(j)*auxfunca
!     enddo
! 
!     return
! 
!     end subroutine gaulagSC
! 
! ! fonction integrale pour la quadrature
!     double precision function func(frail,b,frailtypred)
! 
!     use tailles,only:npmax
!     use comon,only:nsujet,ng,nva,ve,t1,g,typeof,nst,stra, &
!     nz1,nz2,zi,betacoef,nbintervR,ttt,etaR,etaD,betaR,betaD
!     
!     implicit none
! 
!     double precision,intent(in)::frail
!     double precision,dimension(npmax),intent(in)::b
!     double precision,dimension(ng),intent(in)::frailtypred
! 
!     integer::i,j,k,n,gg
!     double precision,dimension(ng)::p
!     double precision,dimension(nsujet)::dNij
!     double precision::vet,S0,psum,dN
!     double precision,dimension(-2:npmax)::the1,the2
!     double precision::bbb,su
! 
! ! calcul de S0 d'abord
!     S0 = 0.d0
!     do i=1,nsujet
!         if(nva.gt.0)then
!             vet = 0.d0
!             do j=1,nva
!                 vet =vet + b(npmax-nva+j)*dble(ve(i,j))
!             end do
!             vet = dexp(vet)
!         else
!            vet=1.d0
!         endif
!         if (t1(i).ge.frail) then
!             S0 = S0 + vet
!         endif
!     end do
! 
!     p = 0.d0
!     dNij = 0.d0
!     do i=1,nsujet
!         if(nva.gt.0)then
!             vet = 0.d0
!             do j=1,nva
!                 vet = vet + b(npmax-nva+j)*dble(ve(i,j))
!             end do
!             vet = dexp(vet)
!         else
!            vet=1.d0
!         endif
! 
!         select case(typeof)
!         
!             case(0) ! calcul du risque splines
! 
!             n = (npmax-nva-1)/nst
! 
!             do k=1,n
!                 the1(k-3)=b(k)*b(k)
!                 j = n+k
!                 if (nst.eq.2) then
!                     the2(k-3)=b(j)*b(j)
!                 endif
!             end do
! 
!             if (stra(i).eq.1) then
! !============== fonction de risque
!                 call susps(frail,the1,nz1,su,bbb,zi)
!     
! ! ! le risque du temps maximum n'est pas calculé dans la fonction susps
! !                 if (tps.eq.date(ndate)) then
! !                 bbb = 4.d0*the1(n-2-1)/(zi(n-2)-zi(n-2-1))
! !             endif
! !============ fonction de risque
!             endif
! 
!             if (stra(i).eq.2) then
!                 call susps(frail,the2,nz2,su,bbb,zi)
! !             if (tps.eq.date(ndate)) then
! !                 bbb = 4.d0*the2(n-2-1)/(zi(n-2)-zi(n-2-1))
! !             endif
!             endif
! 
!             case(1) ! calcul du risque piecewise
! 
!             betacoef = 0.d0
!             do k=1,nst*nbintervR
!                 betacoef(k)=b(k)**2
!             end do
! 
!             if (stra(i).eq.1) then
!                 do gg=1,nbintervR
!                     if((frail.ge.(ttt(gg-1))).and.(frail.lt.(ttt(gg))))then
!                         bbb = betacoef(gg)
!                     end if
!                 end do
! !                 if((tps.ge.(ttt(nbintervR))))then
! !                     bbb = betacoef(nbintervR)
! !                 end if
!             endif
! 
!             if (stra(i).eq.2) then
!                 do gg=1,nbintervR
!                     if((frail.ge.(ttt(gg-1))).and.(frail.lt.(ttt(gg))))then
!                         bbb = betacoef(gg+nbintervR)
!                     end if
!                 end do
! !             if((tps.ge.(ttt(nbintervR))))then
! !                 bbb = betacoef(nbintervR+nbintervR)
! !             end if
!             endif
! 
!             case(2) ! calcul du risque weibull
! 
!             if (nst.eq.1) then
!                 betaR = b(1)**2
!                 etaR = b(2)**2
!                 etaD = 0.d0
!                 betaD = 0.d0
!             else
!                 betaR = b(1)**2
!                 etaR = b(2)**2
!                 betaD = b(3)**2
!                 etaD = b(4)**2
!             end if
! 
! !             if (frail.eq.0.d0) frail = 1d-12 ! utile car log(0) => -Inf
! 
!             if (stra(i).eq.1) then
!                 ! ecriture en exp(log) pour virer l'exposant
!                 bbb = (betaR*dexp((betaR-1.d0)*dlog(frail))/(etaR**betaR))
!             endif
! 
!             if (stra(i).eq.2) then
!                 bbb = (betaD*dexp((betaD-1.d0)*dlog(frail))/(etaD**betaD))
!             endif
! 
!         end select
! 
!         if (t1(i).ge.frail) then
!             p(g(i)) = p(g(i)) + vet/S0
!             dNij = frailtypred(g(i))*bbb*vet
!         endif
!         dNij = dNij - bbb*vet
!     end do
!     
!     psum = 0.d0
!     do i=1,ng
!         psum = psum + p(i)*p(i)
!     end do
!     dN = sum(dNij)
! 
!     func = psum*dN
! 
!     return
!     
!     end function func

!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Shared
!=============================================================================
    subroutine ResidusMartingale(b,np,namesfuncres,Resmartingale,frailtypred,frailtyvar,frailtysd)

    use residusM
    use optimres
    use comon

    implicit none
    
    integer::np
    double precision,external::namesfuncres
    double precision,dimension(np),intent(in)::b
    double precision,dimension(ng),intent(out)::Resmartingale
    double precision,dimension(ng),intent(out)::frailtypred,frailtysd,frailtyvar

    
    vecuiRes=0.d0
    moyuiR=0.d0
    varuiR=0.d0
    cares=0.d0
    cbres=0.d0
    ddres=0.d0

! la prediction des effets aleatoires n'est pas la meme pour gamma ou log-normal
    if (logNormal.eq.0) then !gamma frailty
        do indg=1,ng
            post_esp(indg)=(nig(indg)+1/(b(np-nva)*b(np-nva)))/(cumulhaz(indg)+1/(b(np-nva)*b(np-nva)))
            post_SD(indg)=dsqrt((nig(indg)+1/(b(np-nva)*b(np-nva)))/((cumulhaz(indg)+1/(b(np-nva)*b(np-nva)))**2))

            Resmartingale(indg)=nig(indg)-(post_esp(indg))*cumulhaz(indg)

            frailtypred(indg) = post_esp(indg)
            frailtysd(indg) = post_SD(indg)
            frailtyvar(indg) = frailtysd(indg)**2
        end do
    else !log normal frailty
        do indg=1,ng
            vuu=0.9d0
            call marq98res(vuu,1,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)

            if (istopres.eq.1) then
                Resmartingale(indg)=nig(indg)-(dexp(vuu(1)*vuu(1)))*cumulhaz(indg)
                frailtypred(indg) = vuu(1)*vuu(1)
                frailtyvar(indg) = ((2.d0*vuu(1))**2)*vres(1)
                frailtysd(indg) = dsqrt(frailtyvar(indg))
            else
                ! non convergence ou erreur de calcul de la fonction a maximiser
                Resmartingale(indg) = 0.d0
                frailtypred(indg) = 0.d0
                frailtyvar(indg) = 0.d0
                frailtysd(indg) = 0.d0
            endif
        end do
    endif

    end subroutine ResidusMartingale


!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Joint
!=============================================================================

    subroutine ResidusMartingalej(b,np,namesfuncres,Resmartingale,Resmartingaledc,&
    frailtypred,frailtyvar)

    use residusM
    use optimres
    use comon

    implicit none
    
    integer::np
    double precision,external::namesfuncres
    double precision,dimension(np),intent(in)::b
    double precision,dimension(np)::bint
    double precision,dimension(ng),intent(out)::Resmartingale,Resmartingaledc
    double precision,dimension(ng),intent(out)::frailtypred,frailtyvar
    
    bint=b
    ResidusRec=0.d0
    Residusdc=0.d0
    vecuiRes=0.d0
    moyuiR=0.d0
    
    do indg=1,ng
        vuu=0.9d0
        call marq98res(vuu,1,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)
        if (istopres.eq.1) then 
            if (logNormal.eq.0) then !gamma frailty
                ResidusRec(indg)=Nrec(indg)-((vuu(1)*vuu(1)))*Rrec(indg)
                Residusdc(indg)=Ndc(indg)-((vuu(1)*vuu(1))**alpha)*Rdc(indg)
            else!log normal frailty
                ResidusRec(indg)=Nrec(indg)-(dexp(vuu(1)*vuu(1)))*Rrec(indg)
                Residusdc(indg)=Ndc(indg)-(dexp(vuu(1)*vuu(1)*alpha))*Rdc(indg)
            endif
            vecuiRes(indg) = vuu(1)*vuu(1)
            Resmartingale(indg) = ResidusRec(indg)
            Resmartingaledc(indg) = Residusdc(indg)
            frailtypred(indg) = vecuiRes(indg)
            frailtyvar(indg) = ((2.d0*vuu(1))**2)*vres(1)
        else
            ! non convergence ou erreur de calcul de la fonction a maximiser
            Resmartingale(indg) = 0.d0
            Resmartingaledc(indg) = 0.d0
            frailtypred(indg) = 0.d0
            frailtyvar(indg) = 0.d0
        endif
    end do    
    end subroutine ResidusMartingalej
        
!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES  Joint Nested (AK 12/12/2016)
!=============================================================================
    
    subroutine ResidusMartingalej_fam(namesfuncres,Resmartingale,Resmartingaledc,&
                    frailtypred,frailtyvar,frailtypredind,frailtyvarind)

    use residusM
    use optimres
    use comon
    use commun

    implicit none
    
    integer::i,j!,np
    !double precision,dimension(np),intent(in)::b
    double precision,external::namesfuncres
    double precision,dimension(ng),intent(out)::Resmartingale,Resmartingaledc
    double precision,dimension(ng),intent(out)::frailtypred,frailtyvar
    double precision,dimension(ng,ng),intent(out)::frailtypredind,frailtyvarind
    double precision,dimension(:),allocatable::vuuu
    double precision,dimension(:,:),allocatable::H_hess0
    integer::indiv

    cares=0.d0
    cbres=0.d0
    ddres=0.d0
    vecuiRes=0.d0
    moyuiR=0.d0
    indiv = 1 
    
    Resmartingale = Nrec_fam(1:nfam) !
    Resmartingaledc = Ndc_fam(1:nfam) !    

    do indg=1,nfam
        allocate(H_hess0(fsize(indg)+1,fsize(indg)+1))        
        allocate(vuuu(fsize(indg)+1),vres((fsize(indg)+1)*((fsize(indg)+1)+3)/2))
        vuuu=0.9d0

        call marq98res(vuuu,(fsize(indg)+1),nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)

        do i=1,fsize(indg)+1
            do j=i,fsize(indg)+1
                H_hess0(i,j)=vres((j-1)*j/2+i)
            end do
        end do
        do i=1,(fsize(indg)+1)
            do j=1,i-1
                H_hess0(i,j) = H_hess0(j,i)
            end do
        end do
    
        if (istopres.eq.1) then
            do i=1,fsize(indg)
                Resmartingale(indiv) =Nrec(indiv) -(((vuuu(1)*vuuu(1))**xi)*&
                        (vuuu(1+i)*vuuu(1+i)))*cumulhaz1(indg,i)
                Resmartingaledc(indiv) = Ndc(indiv) - ((vuuu(1)*vuuu(1))*&
                        ((vuuu(1+i)*vuuu(1+i))**alpha))*cumulhazdc(indg,i)
                  frailtypredind(indg,i) = vuuu(1+i)**2
                indiv = indiv + 1
            end do

            frailtypred(indg) = vuuu(1)**2
            frailtyvar(indg) = ((2.d0*vuuu(1))**2)*H_hess0(1,1)

            do i=1,fsize(indg)
                frailtyvarind(indg,i) = ((2.d0*vuuu(1+i))**2)*H_hess0(1+i,1+i)
            end do
        else
                 Resmartingale(indiv:(indiv+fsize(indg)-1)) = 0.d0
Resmartingaledc(indiv:(indiv+fsize(indg)-1))  = 0.d0

indiv = indiv + fsize(indg)
            frailtypredind(indg,:) = 0.d0
            frailtyvarind(indg,:) = 0.d0
            frailtyvar(indg) = 0.d0
            frailtypred(indg) = 0.d0
        end if
        deallocate(vuuu,vres,H_hess0)!,I_hess,H_hess)
    end do
  end subroutine ResidusMartingalej_fam


!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Nested
!=============================================================================
    
    subroutine ResidusMartingalen(namesfuncres,Resmartingale,frailtypred,maxng,frailtypredg,&
    frailtyvar,frailtyvarg,frailtysd,frailtysdg)

    use residusM
    use optimres
    !use comon,only:alpha,eta
    use commun
    use tailles,only:nssgbyg

    implicit none
    
    integer::i,j,maxng
    double precision,external::namesfuncres
    double precision,dimension(nssgbyg),intent(out)::Resmartingale
    double precision,dimension(ngexact),intent(out)::frailtypred,frailtysd,frailtyvar
    double precision,dimension(ngexact,maxng),intent(out)::frailtypredg,frailtysdg,frailtyvarg
    double precision,dimension(:),allocatable::vuuu
    double precision,dimension(:,:),allocatable::H_hess0
    integer:: indiv

    cares=0.d0
    cbres=0.d0
    ddres=0.d0
    indiv = 1
    Resmartingale = mid(1:ngexact) !mid

    do indg=1,ngexact

        allocate(H_hess0(n_ssgbygrp(indg)+1,n_ssgbygrp(indg)+1))
        
        allocate(vuuu(n_ssgbygrp(indg)+1),vres((n_ssgbygrp(indg)+1)*((n_ssgbygrp(indg)+1)+3)/2))

        vuuu=0.9d0

        call marq98res(vuuu,(n_ssgbygrp(indg)+1),nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)

        do i=1,n_ssgbygrp(indg)+1
            do j=i,n_ssgbygrp(indg)+1
                H_hess0(i,j)=vres((j-1)*j/2+i)
            end do
        end do
        do i=1,(n_ssgbygrp(indg)+1)
            do j=1,i-1
                H_hess0(i,j) = H_hess0(j,i)
            end do
        end do

        if (istopres.eq.1) then

            do i=1,n_ssgbygrp(indg)
            
                Resmartingale(indiv) =  mij(indg,i)- ((vuuu(1)*vuuu(1+i))**2)*cumulhaz1(indg,i)
                frailtypredg(indg,i) = vuuu(1+i)**2
                indiv  = indiv + 1
            end do

            frailtypred(indg) = vuuu(1)**2

            frailtysd(indg) = dsqrt(((2.d0*vuuu(1))**2)*H_hess0(1,1)) ! correction de la variance le 2 est dans le carre
            frailtyvar(indg) = ((2.d0*vuuu(1))**2)*H_hess0(1,1)

            do i=1,n_ssgbygrp(indg)
                frailtysdg(indg,i) = dsqrt(((2.d0*vuuu(1+i))**2)*H_hess0(1+i,1+i))
                frailtyvarg(indg,i) = ((2.d0*vuuu(1+i))**2)*H_hess0(1+i,1+i)
            end do

        else
            Resmartingale(indiv) = 0.d0
            frailtypredg(indg,:) = 0.d0
            frailtysdg(indg,:) = 0.d0
            frailtyvarg(indg,:) = 0.d0
            frailtysd(indg) = 0.d0
            frailtyvar(indg) = 0.d0
        end if

        deallocate(vuuu,vres,H_hess0)!,I_hess,H_hess)

    end do

    end subroutine ResidusMartingalen
    


!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Additive
!=============================================================================
    
    subroutine ResidusMartingalea(b,np,namesfuncres,Resmartingale,frailtypred,frailtyvar,frailtysd,&
    frailtypred2,frailtyvar2,frailtysd2,frailtycov)

    use parameters
    use residusM,only:indg,cumulhaz
    use optimres
    !use comon,only:alpha,eta,nst,nig
    use comon,only:nsujet,g,stra,nt1,nva,ve,typeof!,H_hess
    use additiv,only:ve2,ngexact,ut1,ut2,mid

    implicit none
    
    integer::np,k,ip,i,j,ier,istop,ni
    double precision::vet,ca,cb,dd,rl
    double precision,dimension(np),intent(in)::b
    double precision,external::namesfuncres    
    double precision,dimension(ngexact),intent(out)::Resmartingale
    double precision,dimension(ngexact),intent(out)::frailtypred,frailtysd,frailtyvar,frailtycov
    double precision,dimension(ngexact),intent(out)::frailtypred2,frailtysd2,frailtyvar2
    double precision,dimension(2,2)::H_hess0
    double precision,dimension(2)::vu
    double precision,dimension(2*(2+3)/2)::v

    vet=0.d0
    H_hess0=0.d0
    
    Resmartingale = mid

    do indg = 1,ngexact

        vu=0.0d0
        v=0.d0
        call marq98res(vu,2,ni,v,rl,ier,istop,ca,cb,dd,namesfuncres)

         do i=1,2
             do j=i,2
                 H_hess0(i,j)=v((j-1)*j/2+i)
             end do
         end do
         H_hess0(2,1) = H_hess0(1,2)
        
        

        do k=1,nsujet
            if(nva.gt.0 .and.g(k).eq.indg)then
                vet = 0.d0 
                do ip = 1,nva
                    vet = vet + b(np-nva +ip)*ve(k,ip)
                end do

                vet = dexp(vet)
            else
                vet=1.d0
            endif
            if(typeof==0) then
                if(g(k) == indg)then    
                    if(stra(k).eq.1)then
                        Resmartingale(indg) = Resmartingale(indg) - ut1(nt1(k)) * &
                            dexp(vu(1) + vu(2) * ve2(k,1) + dlog(vet))
                    end if
                    if(stra(k).eq.2)then
                        Resmartingale(indg) = Resmartingale(indg) - ut2(nt1(k)) * dexp(vu(1) &
                        + vu(2) * ve2(k,1) + dlog(vet))
                    end if
                end if
            else
                if(g(k) == indg)then    
                    Resmartingale(indg) = Resmartingale(indg) - cumulhaz(g(k)) * &
                    dexp(vu(1) + vu(2) * ve2(k,1) + dlog(vet))
                end if
            end if
        end do
    
        frailtypred(indg) = vu(1)
        frailtypred2(indg) = vu(2)

        if(istop==1) then
            frailtyvar(indg) = H_hess0(1,1)
            frailtysd(indg) = dsqrt(H_hess0(1,1))
            
            frailtyvar2(indg) = H_hess0(2,2)
            frailtysd2(indg) = dsqrt(H_hess0(2,2))
            
            frailtycov(indg) = H_hess0(1,2)
            
        else
            frailtysd(indg) = 0.d0
            frailtyvar(indg) = 0.d0 
            frailtysd2(indg) = 0.d0
            frailtyvar2(indg) = 0.d0 
            frailtycov(indg) = 0.d0
        end if
    end do

    end subroutine ResidusMartingalea
    

!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Joint multive
!=============================================================================
        
    subroutine Residus_Martingale_multive(b,np,names_func_res,Res_martingale,Res_martingaledc,Res_martingale2,&
    frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtyCorr)
    

    use residusMmultiv
!    use optim
    use optimres
    use comonmultiv

    implicit none
    
    integer::np
    double precision,external::names_func_res
    double precision,dimension(np),intent(in)::b
    double precision,dimension(np)::bint
    double precision,dimension(ng),intent(out)::Res_martingale,Res_martingaledc,Res_martingale2
    double precision,dimension(ng),intent(out)::frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtyCorr    
    double precision::ca,cb,dd,rl
    integer::ni,ier
    
    bint=b
    ResidusRec=0.d0
    Residusdc=0.d0
    ResidusRec2=0.d0!Residusmeta
    moyuiR=0.d0
  
    do indg=1,ng
        ca=0.d0
        cb=0.d0
        dd=0.d0
        ni=0
        vuu=0.1d0
                
        call marq98res(vuu,2,ni,vres,rl,ier,istopres,ca,cb,dd,names_func_res)
        
        ResidusRec(indg)=Nrec(indg)-dexp(vuu(1))**Rrec(indg)
        Residusdc(indg)=Ndc(indg)-dexp(vuu(1)*alpha1+vuu(2)*alpha2)*Rdc(indg)
        ResidusRec2(indg)=Nrec2(indg)-dexp(vuu(2))*Rrec2(indg)
        
        Res_martingale(indg) = ResidusRec(indg)
        Res_martingaledc(indg) = Residusdc(indg)
        Res_martingale2(indg) = ResidusRec2(indg)    
 
        frailtypred(indg) = vuu(1)
        frailtypred2(indg) = vuu(2)        

        frailtyvar(indg) = vres(1)
        frailtyvar2(indg) = vres(3)    
        frailtyCorr(indg) = vres(2)/dsqrt(vres(1)*vres(2))
    end do    
 
    end subroutine Residus_Martingale_multive    


    

    !=============================================================================
!                       CALCUL DES RESIDUS  Joint bivarie : longitudinal et deces
!=============================================================================
    
    subroutine Residusj_biv(b,np,namesfuncres,Resmartingaledc,ResLongi_cond,ResLongi_cond_st,&
                            ResLongi_marg,ResLongi_chol,Pred_yy,    re_pred)!
    
    use residusM
    use optimres
    !use comon,only:ut,utt,netadc
    use comon,only:ng,nsujety,etaydc,yy,nmesy,nb1,vey,varcov_marg,&
    nva3,sum_mat,link,t1dc,vey,npp,res_ind
    use donnees_indiv,only:X2cur,Z1cur
    use optim
    
    implicit none
    
    integer::np,j,i,ij,ier,k,jj
    double precision,external::namesfuncres
    double precision,dimension(np),intent(in)::b
    double precision,dimension(np)::bint
    double precision,dimension(ng),intent(out)::Resmartingaledc
    double precision,dimension(nsujety),intent(out)::ResLongi_marg,ResLongi_chol,&
    ResLongi_cond,ResLongi_cond_st
    double precision,dimension(nsujety,2),intent(out)::Pred_yy
    double precision,dimension(ng,nb1+1),intent(out)::re_pred
    double precision,dimension(nb1) ::zet_vec
    double precision,dimension(nb1) :: b_pred
    double precision,dimension((nb1-1)*nb1/2+nb1) :: vres_inv
    double precision :: ep,eps
    double precision,dimension(nb1,nb1)::mat_vres_inv,mat_vres
    double precision,dimension(:,:),allocatable::v_rim,v_rim_chol,Varcov_inv,&
    varcov_chol,V_rim_inv
    double precision,dimension(:),allocatable:: matv,matv2
    double precision,dimension(nva3,nva3):: sum_mat_inv
    double precision,dimension(1)::current_meanres
     !   double precision,dimension(1)::residus_sd
    
    res_ind = 1
    bint=b
    ResLongi_cond_st=0.d0
    ResLongi_cond=0.d0
    ResLongi_marg=0.d0
    ResLongi_chol = 0.d0
    Pred_yy= 0.d0
    Residusdc=0.d0
    it_res = 1
    
    do indg=1,ng
    cares=0.d0
    cbres=0.d0
    ddres=0.d0
    ierres = 0
    nires=0
    vuu=0.1d0
        
    call marq98res(vuu,nb1,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)
 
        if (istopres.eq.1) then
                if(link.eq.1) then
                    Residusdc(indg)=Ndc(indg)-exp(dot_product(etaydc,vuu))*Rdc(indg)
                     
                else
                    X2cur(1,1) = 1.d0
        X2cur(1,2) =t1dc(indg)
        if((nva3-2).gt.0) then
            do k=3,nva3
                    X2cur(1,k) = dble(vey(it_res,k))
                end do
        end if
                    Z1cur(1,1) = 1.d0
                        if(nb1.eq.2)     Z1cur(1,2) =t1dc(indg)
                            current_meanres = 0.d0
    
        current_meanres(1) =dot_product(X2cur(1,1:nva3),b((npp-nva3+1):npp))+&
        dot_product(Z1cur(1,1:nb1),vuu(1:nb1))
    
                        Residusdc(indg)=Ndc(indg)-(exp(etaydc(1)*current_meanres(1)))*Rdc(indg)
                end if
                b_pred(1)  = vuu(1)
            if(nb1.eq.2)b_pred(2)  = vuu(2)
    
            vres_inv = 0.d0
            do i=1,nb1
                do j=i,nb1
                    ij=(j-1)*j/2+i
                    vres_inv(ij)=vres(ij)
                end do
            end do
            ier = 0
            ep = 1.d-10
            call dsinvj(vres_inv,nb1,ep,ier)
    
            mat_vres_inv=0.d0
            mat_vres = 0.d0
            do j=1,nb1
                do k=1,nb1
                    if (k.ge.j) then
                        mat_vres_inv(j,k)=vres_inv(j+k*(k-1)/2)
                        mat_vres(j,k) = vres(j+k*(k-1)/2)
                    else
                        mat_vres_inv(j,k)=vres_inv(k+j*(j-1)/2)
                        mat_vres(j,k)=vres(k+j*(j-1)/2)
                    end if
                end do
            end do
    
   
      
    allocate(V_rim_chol(nmesy(indg),nmesy(indg)),v_rim(nmesy(indg),nmesy(indg)),Varcov_inv(nmesy(indg),nmesy(indg)),&
            varcov_chol(nmesy(indg),nmesy(indg)),V_rim_inv(nmesy(indg),nmesy(indg)))
   
                        allocate(matv(nva3*(nva3+1)/2),matv2(nmesy(indg)*(nmesy(indg)+1)/2))
                        matv = 0.d0
                        do j=1,nva3
    do k=j,nva3
        jj=j+k*(k-1)/2
        matv(jj)=sum_mat(j,k)
    
        end do
    end do
    ier = 0
    eps = 1.d-10
    
    
        call dsinvj(matv,nva3,eps,ier)
    
        sum_mat_inv=0.d0
        do j=1,nva3
                do k=1,nva3
                            if (k.ge.j) then
                sum_mat_inv(j,k)=matv(j+k*(k-1)/2)
            else
                sum_mat_inv(j,k)=matv(k+j*(j-1)/2)
            end if
            end do
                end do
    
                        V_rim =Varcov_marg(it_res:it_res+nmesy(indg)-1,1:nmesy(indg))&
                                                - MATMUL(Matmul(vey(it_res:it_res+nmesy(indg)-1,1:nva3),sum_mat_inv),&
                                                        Transpose(vey(it_res:it_res+nmesy(indg)-1,1:nva3))) !varcov_marg
    
                                matv2 = 0.d0
                        do j=1,nmesy(indg)
    do k=j,nmesy(indg)
        jj=j+k*(k-1)/2
        matv2(jj)=V_rim(j,k)
    
        end do
    end do
    ier = 0
    eps = 1.d-10
    
   
        call dsinvj(matv2,nmesy(indg),eps,ier)
    
        V_rim_inv=0.d0
        do j=1,nmesy(indg)
                do k=1,nmesy(indg)
                            if (k.ge.j) then
                V_rim_inv(j,k)=matv2(j+k*(k-1)/2)
            else
                V_rim_inv(j,k)=matv2(k+j*(j-1)/2)
            end if
            end do
                end do
    
        !       call choldc(nmesy(indg),V_rim,V_rim_chol)
                        call cholesky_sub( V_rim_inv,nmesy(indg))
   
                        v_rim_chol = V_rim_inv! transpose( V_rim_inv)
    
    do j= it_res,(it_res+nmesy(indg)-1)
        zet_vec(1:nb1) = Zet(j,1:nb1)
        ResLongi_marg(j) = yy(j) - XbetaY_res(1,j)
           
           ResLongi_cond(j) = yy(j) - XbetaY_res(1,j) -dot_product(zet_vec(1:nb1),b_pred(1:nb1))
            Pred_yy(j,1) = XbetaY_res(1,j) +dot_product(zet_vec(1:nb1),b_pred(1:nb1))
       
                        
        Pred_yy(j,2) =  XbetaY_res(1,j)
                        
    end do
    
                                matv2 = 0.d0
                        do j=1,nmesy(indg)
    do k=j,nmesy(indg)
        jj=j+k*(k-1)/2
        matv2(jj)=Varcov_marg((it_res+j-1),k)
        end do
    end do
    ier = 0
    eps = 1.d-10
    
    
        call dsinvj(matv2,nmesy(indg),eps,ier)
    
        Varcov_inv=0.d0
        do j=1,nmesy(indg)
                do k=1,nmesy(indg)
                            if (k.ge.j) then
                Varcov_inv(j,k)=matv2(j+k*(k-1)/2)
            else
                Varcov_inv(j,k)=matv2(k+j*(j-1)/2)
            end if
            end do
                end do
                
                        call cholesky_sub(Varcov_inv,nmesy(indg))
                varcov_chol = transpose(Varcov_inv)!Varcov_inv!
    
    ResLongi_cond_st( it_res:it_res+nmesy(indg)-1) = Matmul(V_rim_chol,&
        ResLongi_cond( it_res:it_res+nmesy(indg)-1))
    
    
            ResLongi_chol(it_res:(it_res+nmesy(indg)-1)) = Matmul(varcov_chol(1:nmesy(indg),1:nmesy(indg)),&
                                                        ResLongi_marg(it_res:(it_res+nmesy(indg)-1)))
                        vecuiRes2(indg,1:nb1) = vuu(1:nb1)
                        vecuiRes2(indg,nb1+1) = 0.d0
            Resmartingaledc(indg) = Residusdc(indg)
                re_pred(indg,:) = vecuiRes2(indg,:)
       deallocate(matv,matv2,v_rim,v_rim_chol,varcov_inv,varcov_chol,V_rim_inv)
        else
            ! non convergence ou erreur de calcul de la fonction a maximiser
            do j= it_res,(it_res+nmesy(indg)-1)
                        Reslongi_cond_st(j) = 0.d0
                        Reslongi_cond(j) = 0.d0
            Reslongi_marg(j) = 0.d0
                        ResLongi_chol(j) = 0.d0
            pred_yy(j,1:2) = 0.d0
            end do
    
            Resmartingaledc(indg) = 0.d0
            re_pred(indg,:) = 0.d0
    
        endif
    it_res = it_res + nmesy(indg)
    
    
    end do
    
    end subroutine Residusj_biv
    
    
    
    
    
    !=============================================================================
     !                   CALCUL DES RESIDUS  Joint trivarie : longitudinal, recurrences et deces
    !=============================================================================
    
    subroutine Residusj_tri(b,np,namesfuncres,Resmartingale,Resmartingaledc,ResLongi_cond,ResLongi_cond_st,&
                                        ResLongi_marg,ResLongi_chol,Pred_yy, re_pred)!
    
    use residusM
    use optimres
    !use comon,only:ut,utt,netar,nsujet
    use comon,only:ng,nsujety,etayr,etaydc,yy,nmesy,nmesrec,netadc,&
         nb1,vey,varcov_marg,nva3,sum_mat,link,t1dc,vey,nea,alpha,npp,res_ind
    use donnees_indiv,only:X2cur,Z1cur
        use optim
    
    implicit none
    
    integer::np,j,i,ij,ier,k,jj
    double precision,external::namesfuncres
    double precision,dimension(np),intent(in)::b
    double precision,dimension(np)::bint
    double precision,dimension(ng),intent(out)::Resmartingaledc,Resmartingale
        double precision,dimension(nsujety),intent(out)::ResLongi_marg,ResLongi_cond,&
        ResLongi_chol,ResLongi_cond_st
        double precision,dimension(nsujety,2),intent(out)::Pred_yy
    double precision,dimension(ng,nb1+1),intent(out)::re_pred
    double precision,dimension(nb1) ::zet_vec
    double precision,dimension(nea) :: b_pred
    double precision,dimension((nb1-1)*nb1/2+nb1) :: vres_inv
    double precision :: ep,eps
    double precision,dimension(nb1,nb1)::mat_vres_inv,mat_vres
        double precision,dimension(:,:),allocatable::v_rim,v_rim_chol,Varcov_inv,varcov_chol,V_rim_inv
        double precision,dimension(:),allocatable:: matv,matv2
        double precision,dimension(nva3,nva3):: sum_mat_inv
        double precision,dimension(1)::current_meanres
    
       ! double precision,dimension(1)::residus_sd
    
    res_ind=1
    bint=b
        ResidusRec=0.d0
        ResLongi_cond_st=0.d0
        ResLongi_cond=0.d0
    ResLongi_marg=0.d0
        ResLongi_chol = 0.d0
    Pred_yy = 0.d0
    Residusdc=0.d0
    it_res = 1
        it_res_rec = 1
    do indg=1,ng
    
            cares=0.d0
        cbres=0.d0
        ddres=0.d0
        nires=0
        vuu=0.1d0
        vuu=0.1d0
    
    
    
    call marq98res(vuu,nea,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)
    
    
        !print*,"6:",vres
        !print*,"----------",indg
        if (istopres.eq.1) then
                        if(link.eq.1) then
                            Residusdc(indg)=Ndc(indg)-exp(vuu(nb1+1)*alpha+&
                                    dot_product(etaydc,vuu(1:nb1)))*Rdc(indg)
                            ResidusRec(indg)=Nrec(indg)-exp(vuu(nb1+1)+&
                                    dot_product(etayr,vuu(1:nb1)))*Rrec(indg)
                        else if(link.eq.2) then
                                X2cur(1,1) = 1.d0
                                X2cur(1,2) =t1dc(indg)
                                if((nva3-2).gt.0) then
                                        do k=3,nva3
                                                X2cur(1,k) = dble(vey(it_res,k))
                                        end do
                                end if
                                Z1cur(1,1) = 1.d0
                                if(nb1.eq.2)     Z1cur(1,2) =t1dc(indg)
                                current_meanres = 0.d0
                                current_meanres(1) =dot_product(X2cur(1,1:nva3),b((npp-nva3+1):npp))+&
                        dot_product(Z1cur(1,1:nb1),vuu(1:nb1))
                        Residusdc(indg)=Ndc(indg)-(exp(vuu(3)*alpha+etaydc(1)*current_meanres(1)))*Rdc(indg)
                                ResidusRec(indg)=Nrec(indg)-(exp(vuu(3)+etayr(1)**current_meanres(1)))*Rrec(indg)
    
                        end if
    
                b_pred(1:nb1)  = vuu(1:nb1)
                b_pred(nb1+1) = vuu(nb1+1)
    
            vres_inv = 0.d0
            do i=1,nb1
                do j=i,nb1
                    ij=(j-1)*j/2+i
                    vres_inv(ij)=vres(ij)
                end do
            end do
            ier = 0
            ep = 1.d-10
            call dsinvj(vres_inv,nb1,ep,ier)
    
            mat_vres_inv=0.d0
            mat_vres = 0.d0
            do j=1,nb1
                do k=1,nb1
                    if (k.ge.j) then
                        mat_vres_inv(j,k)=vres_inv(j+k*(k-1)/2)
                        mat_vres(j,k) = vres(j+k*(k-1)/2)
                    else
                        mat_vres_inv(j,k)=vres_inv(k+j*(j-1)/2)
                        mat_vres(j,k)=vres(k+j*(j-1)/2)
                    end if
                end do
            end do
    
    
                allocate(v_rim(nmesy(indg),nmesy(indg)),V_rim_chol(nmesy(indg),nmesy(indg)),&
                Varcov_inv(nmesy(indg),nmesy(indg)),&
                 varcov_chol(nmesy(indg),nmesy(indg)),V_rim_inv(nmesy(indg),nmesy(indg)))
    
    
                        allocate(matv(nva3*(nva3+1)/2),matv2(nmesy(indg)*(nmesy(indg)+1)/2))
                        matv = 0.d0
                        do j=1,nva3
    do k=j,nva3
        jj=j+k*(k-1)/2
        matv(jj)=sum_mat(j,k)
    
        end do
    end do
    ier = 0
    eps = 1.d-10
    
    
        call dsinvj(matv,nva3,eps,ier)
    
        sum_mat_inv=0.d0
        do j=1,nva3
                do k=1,nva3
                            if (k.ge.j) then
                sum_mat_inv(j,k)=matv(j+k*(k-1)/2)
            else
                sum_mat_inv(j,k)=matv(k+j*(j-1)/2)
            end if
            end do
                end do
    
               V_rim =Varcov_marg(it_res:it_res+nmesy(indg)-1,1:nmesy(indg)) - &
                       MATMUL(Matmul(vey(it_res:it_res+nmesy(indg)-1,1:nva3),sum_mat_inv),&
                     Transpose(vey(it_res:it_res+nmesy(indg)-1,1:nva3))) !varcov_marg
    
    
    
                                matv2 = 0.d0
                        do j=1,nmesy(indg)
    do k=j,nmesy(indg)
        jj=j+k*(k-1)/2
        matv2(jj)=V_rim(j,k)
    
        end do
    end do
    ier = 0
    eps = 1.d-10
    
    
        call dsinvj(matv2,nmesy(indg),eps,ier)
    
        V_rim_inv=0.d0
        do j=1,nmesy(indg)
                do k=1,nmesy(indg)
                            if (k.ge.j) then
                V_rim_inv(j,k)=matv2(j+k*(k-1)/2)
            else
                V_rim_inv(j,k)=matv2(k+j*(j-1)/2)
            end if
            end do
                end do
    
    
        !       call choldc(nmesy(indg),V_rim,V_rim_chol)
                        call cholesky_sub( V_rim_inv,nmesy(indg))
    
                        v_rim_chol = V_rim_inv! transpose( V_rim_inv)
    
    do j= it_res,(it_res+nmesy(indg)-1)
        zet_vec(1:nb1) = Zet(j,1:netadc)
        
        ResLongi_marg(j) = yy(2) - XbetaY_res(1,j)
    
        ResLongi_cond(j) =  yy(j) - XbetaY_res(1,j) -dot_product(zet_vec(1:nb1),b_pred(1:nb1))
        Pred_yy(j,1) = XbetaY_res(1,j) +dot_product(zet_vec(1:nb1),b_pred(1:nb1))
       
        Pred_yy(j,2) =  XbetaY_res(1,j)
    end do
                                                matv2 = 0.d0
                        do j=1,nmesy(indg)
    do k=j,nmesy(indg)
        jj=j+k*(k-1)/2
        matv2(jj)=Varcov_marg(j,k)
    
        end do
    end do
    ier = 0
    eps = 1.d-10
    
    
        call dsinvj(matv2,nmesy(indg),eps,ier)
    
        Varcov_inv=0.d0
        do j=1,nmesy(indg)
                do k=1,nmesy(indg)
                            if (k.ge.j) then
                Varcov_inv(j,k)=matv2(j+k*(k-1)/2)
            else
                Varcov_inv(j,k)=matv2(k+j*(j-1)/2)
            end if
            end do
                end do
    
    
                        call cholesky_sub(Varcov_inv,nmesy(indg))
                        varcov_chol = transpose(Varcov_inv)
    
                            ResLongi_cond_st( it_res:it_res+nmesy(indg)-1) = Matmul(Varcov_inv,&
                ResLongi_cond( it_res:it_res+nmesy(indg)-1))
    
    
            ResLongi_chol(it_res:(it_res+nmesy(indg)-1)) = Matmul(v_rim_chol(1:nmesy(indg),1:nmesy(indg)),&
                                                        ResLongi_marg(it_res:(it_res+nmesy(indg)-1)))
    
            vecuiRes2(indg,:) = vuu
            Resmartingaledc(indg) = Residusdc(indg)
                        Resmartingale(indg) = ResidusRec(indg)
            ! ResLongi_marg(it_res:(it_res+nmesy(indg)-1))= ResidusLongi(it_res:(it_res+nmesy(indg)-1))
            re_pred(indg,:) = vecuiRes2(indg,:)
            !Pred_ymarg(it_res:(it_res+nmesy(indg)-1))= Pred_y(it_res:(it_res+nmesy(indg)-1))
    
    
            deallocate(matv,matv2,v_rim,v_rim_chol,varcov_inv,varcov_chol,V_rim_inv)
        else
            ! non convergence ou erreur de calcul de la fonction a maximiser
            do j= it_res,(it_res+nmesy(indg)-1)
                        Reslongi_cond_st(j) = 0.d0
                        Reslongi_cond(j) = 0.d0
            Reslongi_marg(j) = 0.d0
                        ResLongi_chol(j) = 0.d0
            pred_yy(j,1:2) = 0.d0
            end do
            Resmartingale(indg) = 0.d0
            Resmartingaledc(indg) = 0.d0
            re_pred(indg,:) = 0.d0
    
        endif
    it_res = it_res + nmesy(indg)
        it_res_rec = it_res_rec + nmesrec(indg)
    
    end do
    
    end subroutine Residusj_tri 
    
    
    
    
    
        !=============================================================================
     !                   CALCUL DES RESIDUS  Joint univarie : longitudinal (non-linear)
    !=============================================================================
    
    subroutine Residus_uni(b,np,namesfuncres)!
    
    use residusM
    use optimres
    use comon,only:ng,nmesy,nva3,nva4,&
         vey,nva3,nva4,vey,nea,npp,res_ind,H_hess_GH,b_paGH !nsujety,yy,nb1,I_hess_GH,ut,utt
    use donnees_indiv,only:mu,b1 !X2cur,Z1cur
        use optim
    !  use ParametresPourParallelisation
    implicit none
    
    integer::np,j,i,k !ij,ier,jj
    double precision,external::namesfuncres
    double precision,dimension(np),intent(in)::b
    double precision,dimension(np)::bint
    double precision :: finddet !ep,eps
       ! double precision,dimension(:),allocatable:: matv,matv2
    !double precision,dimension(nea,nea)::H_hess_GH_inv
       ! double precision,dimension(1)::residus_sd
    
    res_ind=1
    bint=b
     
    it_res = 1
 
 do indg=1,ng
    
            cares=0.d0
        cbres=0.d0
        ddres=0.d0
        nires=0
        vuu=0.1d0
    
        mu = 0.d0
            mu(1:nmesy(indg),1) = matmul(vey(it_res:it_res+nmesy(indg)-1,&
            1:(nva3)),b1((npp-nva4-nva3+1):(npp-nva4)))
    
            mu(1:nmesy(indg),2) = matmul(vey(it_res:it_res+nmesy(indg)-1,&
            (nva3+1):(nva3+nva4)),b1((npp-nva4+1):npp))
  
    call marq98res(vuu,nea,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)

!if(CurrentProcessorID.eq.0.and.indg.eq.1) write(*,*)indg,vuu,istopres 
  !  if(indg.eq.1)write(*,*)vuu
    !stop
           do i=1,nea
             do j=i,nea
                 H_hess_GH(i,j)=vres((j-1)*j/2+i)
             end do
         end do
         do i=2,nea
       do j=1,i-1
            H_hess_GH(i,j)=H_hess_GH(j,i)
       end do
    end do
    

 
        !print*,"6:",vres
        !print*,"----------",indg
        if (istopres.eq.1) then
                      
    
            
    
    b_paGH(indg,1:nea) = vuu 
    
!    allocate(matv2(nea*(nea+1)/2))
!                         matv2 = 0.d0
!                        do j=1,nea
!    do k=j,nea
!        jj=j+k*(k-1)/2
!        matv2(jj)=H_hess_GH(j,k)
!        end do
!    end do
!    ier = 0
!    eps = 1.d-10
!    
!    
!        call dsinvj(matv2,nea,eps,ier)
!    
!        H_hess_GH_inv=0.d0
!        do j=1,nea
!                do k=1,nea
!                            if (k.ge.j) then
!                H_hess_GH_inv(j,k)=matv2(j+k*(k-1)/2)
!            else
!                H_hess_GH_inv(j,k)=matv2(k+j*(j-1)/2)
!            end if
!            end do
!                end do
!                
!    H_hess_GH = H_hess_GH_inv            !H_i
    
          call cholesky_sub( H_hess_GH,nea) !B_i
     
!         allocate(matv2(nea*(nea+1)/2))
!                         matv2 = 0.d0
!                        do j=1,nea
!    do k=j,nea
!        jj=j+k*(k-1)/2
!        matv2(jj)=H_hess_GH(j,k)
!        end do
!    end do
!    ier = 0
!    eps = 1.d-10
!    
!    
!        call dsinvj(matv2,nea,eps,ier)
!    
!        H_hess_GH_inv=0.d0
!        do j=1,nea
!                do k=1,nea
!                            if (k.ge.j) then
!                H_hess_GH_inv(j,k)=matv2(j+k*(k-1)/2)
!            else
!                H_hess_GH_inv(j,k)=matv2(k+j*(j-1)/2)
!            end if
!            end do
!                end do
     
!    call inverse(H_hess_GH,H_hess_GH_inv,nea)     !B_i^-1

!     if(indg.eq.5) then 
!     write(*,*)'res'
!     write(*,*)H_hess_GH
!     write(*,*) vuu 
!     write(*,*)vres
!     stop
!     end if
     
         b_paGH(indg,(nea+1)) =  finddet(H_hess_GH,nea)

    
           do j=1,nea
                do k=1,j
                 b_paGH(indg,nea+1+k+j*(j-1)/2) = H_hess_GH(j,k)
                
                end do
                end do
!        deallocate(matv2)

            
            
            
        else
            ! non convergence ou erreur de calcul de la fonction a maximiser
          !write(*,*)indg,istopres  
           b_paGH(indg,:) = 0.d0
    
        endif
        
    it_res = it_res + nmesy(indg)
     

    end do
    
    end subroutine Residus_uni
    
    
      !=============================================================================
     !                     Joint trivariate non-linear : longitudinal, recurrences and death
    !                        Posterior random effects
    !=============================================================================
    
    subroutine PostRE_triNL(b,np,namesfuncres, re_pred)!
    
    use residusM
    use optimres
    use comon,only:ng,nmesy,nmesrec,nea,&
         vey,nva3,H_hess_GH,nva4 !I_hess_GH,nb1
    use donnees_indiv,only:mu
        use optim
!    use ParametresPourParallelisation
    implicit none
    
    integer::np,j,i,k !ij,ier,jj
    double precision,external::namesfuncres
    double precision,dimension(np),intent(in)::b
    double precision,dimension(ng,nea+1+nea + (nea*(nea-1))/2),intent(out)::re_pred
    double precision :: finddet
    
    
   
       it_res = 1
        it_res_rec = 1
    do indg=1,ng
    
            cares=0.d0
        cbres=0.d0
        ddres=0.d0
        nires=0
        vuu=0.1d0
        vuu=0.1d0
    
      mu = 0.d0
            mu(1:nmesy(indg),1) = matmul(vey(it_res:it_res+nmesy(indg)-1,&
            1:(nva3)),b((np-nva4-nva3+1):(np-nva4)))
    
            mu(1:nmesy(indg),2) = matmul(vey(it_res:it_res+nmesy(indg)-1,&
            (nva3+1):(nva3+nva4)),b((np-nva4+1):np))
   
    call marq98res(vuu,nea,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)

           do i=1,nea
             do j=i,nea
                 H_hess_GH(i,j)=vres((j-1)*j/2+i)
             end do
         end do
         do i=2,nea
       do j=1,i-1
            H_hess_GH(i,j)=H_hess_GH(j,i)
       end do
    end do
    

        if (istopres.eq.1) then
     
        
    re_pred(indg,1:nea) = vuu 
    
                
     call cholesky_sub( H_hess_GH,nea)
     
    
    re_pred(indg,(nea+1)) =  finddet(H_hess_GH,nea)
    
           do j=1,nea
                do k=1,j
                 re_pred(indg,nea+1+k+j*(j-1)/2) = H_hess_GH(j,k)
                
                end do
                end do
            
            
        else
     
     re_pred(indg,:) = 0.d0
        endif
    it_res = it_res + nmesy(indg)
        it_res_rec = it_res_rec + nmesrec(indg)
    
    end do

    end subroutine PostRE_triNL
    
    