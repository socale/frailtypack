!===================================================================
!========================          FUNCPA       ====================
    double precision function funcpaasplines(b,np,id,thi,jd,thj,k0)

    use tailles
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
    mm3,mm2,mm1,mm,im3,im2,im1,im,date,zi,pe,effet,nz1,nz2,stra, &
    c,nt0,nt1,nsujet,nva,ndate,nst,auxig,ng,ve,g,nig,indictronq, &
     resnonpen,indic_tronc!,t0,t1
    use additiv,only: &
    dut1,dut2,ut1,ut2,correl,ngexact,ve2,sigma2,tau2,cov,mid,betaaux!,invD,rho
    use residusM,only:som_Xbeta

    implicit none

    integer::n,np,id,jd,i,j,k,vj,cptg,ig,ip
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,pe1,pe2,som1,som2,res,vet,h1
    double precision,dimension(-2:npmax)::the1,the2
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::integrale1,funcaux
    double precision,dimension(2)::k0

    real,parameter::pi = 3.1415926535
    integer::restar,nf
!****** pour la maximisation avec marq98aux
    integer :: npaux,niaux,ieraux,istopaux
    DOUBLE PRECISION::resaux
    double precision,dimension((npmax*(npmax+3)/2))::vaux
    double precision,dimension(npmax)::baux
!****** derivanal
    double precision , dimension(ngmax)::res1,res2,res3,res4
    double precision , dimension(ngmax)::res5,res6
!      common /derivanal/res1,res2,res3,res4,res5,res6,res8

!****** u_tilde
    double precision  :: u_tilde,v_tilde
    double precision,external::funcpao

!******************

    j=0
    h1=0.d0
    res=0.d0
    vet=0.d0
    resnonpen=0.d0
    som1=0.d0
    som2=0.d0
    pe=0.d0
    pe1=0.d0
    pe2=0.d0
    restar = 0
    nf = 1

    do i=1,np
        bh(i)=b(i)
    end do

    if (id.ne.0) bh(id)=bh(id)+thi 
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    if(nst.eq.2)then
        n = (nz1+2 +nz2+2)/nst
    else
        n = nz1+2
    endif

    do i=1,n
        the1(i-3)=(bh(i))*(bh(i))
        j = n+i 
        if (nst.eq.2) then
            the2(i-3)=(bh(j))*(bh(j))
        endif
    end do

    dut1=0.d0
    dut2=0.d0
    ut1=0.d0
    ut2=0.d0

    if(effet.eq.1) then
!terme de correlation entre 2 frailties/avec contraintes
        sigma2 = (bh(np-nva-1)*bh(np-nva-1)) ! variance intercept
        tau2 = (bh(np-nva)*bh(np-nva))  ! variance traitement * groupe
        cov=dsqrt(sigma2*tau2)*(2.d0 * dexp(bh(np-nva-2))/(1.d0+dexp(bh(np-nva-2)))-1.d0)
    endif

    vj = 0
    som1 = 0.d0
    dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
    ut1(0) = 0.d0
    ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))

    if (nst.eq.2) then
        som2 = 0.d0
        dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
        ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
        ut2(0) = 0.d0
    endif

    do i=2,ndate-1
        do k = 2,n-2
            if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                    som1 = som1+the1(j-4)
                    som2 = som2+the2(j-4)
                    vj  = j
                endif
            endif
        end do

        ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i)) &
        +(the1(j-1)*im1(i))+(the1(j)*im(i))
        dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i)) &
        +(the1(j-1)*mm1(i))+(the1(j)*mm(i))

        if(nst.eq.2)then
            ut2(i) = som2 +(the2(j-3)*im3(i))+(the2(j-2)*im2(i)) &
            +(the2(j-1)*im1(i))+(the2(j)*im(i))
            dut2(i) = (the2(j-3)*mm3(i))+(the2(j-2)*mm2(i)) &
            +(the2(j-1)*mm1(i))+(the2(j)*mm(i))
        endif
    end do

    i = n-2
    h1 = (zi(i)-zi(i-1))
    ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
    dut1(ndate) = (4.d0*the1(i-1)/h1)
!AD:the1(i-4) en th2(i-4)
    if(nst.eq.2)then
        ut2(ndate)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)
        dut2(ndate) = (4.d0*the2(i-1)/h1)
    endif
!AD:
!--------------------------------------------------------
!----------calcul de la vraisemblance ------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------cc

    do ig=1,ngexact!ng!
        res1(ig) = 0.d0
        res2(ig) = 0.d0
    end do
    cpt=0
!*******************************************
!----- sans effet aleatoire dans le modele
!*******************************************

    if (effet.eq.0) then

        do i=1,nsujet

            cpt(g(i))=cpt(g(i))+1

            if(nva.gt.0)then
                vet = 0.d0
                do j=1,nva
                    vet =vet + bh(np-nva+j)*ve(i,j)
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif

            if((c(i).eq.1).and.(stra(i).eq.1))then
                res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
            endif

            if((c(i).eq.1).and.(stra(i).eq.2))then
                res2(g(i)) = res2(g(i))+dlog(dut2(nt1(i))*vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpaasplines=-1.d9
                goto 123
            end if

            if(stra(i).eq.1)then
                res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet-ut1(nt0(i))*vet
            endif

            if(stra(i).eq.2)then
                res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet-ut2(nt0(i))*vet 
            endif

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpaasplines=-1.d9
                goto 123
            end if

        end do

        res = 0.d0
        cptg = 0

! k indice les groupes
        do k=1,ngexact 
            if(cpt(k).gt.0)then !nb de sujets dans un gpe=nig()
                res = res-res1(k)+ res2(k)
                cptg = cptg + 1 
            endif
        end do

    else

!************************************************************
!-----avec deux effets aleatoires dans le modele et correl=0
!************************************************************
        if(effet.eq.1.and.correl.eq.0)then

            mid=0
            integrale1=0.d0

            do k=1,nsujet
                if(c(k).eq.1)then
                    mid(g(k))=mid(g(k))+1
                endif
            end do

            do k=1,nsujet

                if(nva.gt.0)then
                    vet = 0.d0
                    do ip=1,nva
                        vet =vet + bh(np-nva+ip)*ve(k,ip)
                    end do
                    vet = dexp(vet)
                else
                    vet=1.d0
                endif

                if((c(k).eq.1).and.(stra(k).eq.1))then
                    res2(g(k)) = res2(g(k))+dlog(dut1(nt1(k))*vet)
                endif

                if((c(k).eq.1).and.(stra(k).eq.2))then
                    res2(g(k)) = res2(g(k))+dlog(dut2(nt1(k))*vet)
                endif

                if ((res2(g(k)).ne.res2(g(k))).or.(abs(res2(g(k))).ge. 1.d30)) then
                    funcpaasplines=-1.d9
                    goto 123
                end if

            end do

!=========================================================================
!==== calcul des integrales par transformation de LAPLACE pour chq gpe ===
!=========================================================================

            do ig=1,ng!ngexact
                res3(ig)=0.d0
                res4(ig)=0.d0
                res5(ig)=0.d0
                res6(ig)=0.d0
                auxig=ig
                baux(1)=0.05d0      !initialisation de u_tilde
                baux(2)=0.05d0      !initialisation de v_tilde

!======================================================================
! maximisation pour deux parametres, les 2 effets aleatoires : u et v
!======================================================================
                npaux=2
                niaux=0
                ieraux=0
                istopaux=0
                vaux=0.d0!vecteur derivees 2nd et 1ERES
                funcaux=0.d0!vraisemblance
                resaux=0.d0!vraisemblance

                do ip=1,nva
                    betaaux(ip)= bh(np-nva+ip)
                end do

    call marq98o(baux,npaux,niaux,vaux,resaux,ieraux,istopaux,funcpao)
       if (ieraux .eq.-1) then
          funcpaasplines=-1.d9
          goto 123
       end if
                u_tilde = baux(1)!u_tilde
                v_tilde = baux(2)!v_tilde

                do k=1,nsujet

                    if(nva.gt.0.and.g(k).eq.ig)then
                        vet = 0.d0 
                        do ip=1,nva
                            vet =vet + bh(np-nva+ip)*ve(k,ip)
                        end do
                        vet = dexp(vet)
                    else
                        vet=1.d0
                    endif

                    if(g(k).eq.ig)then
                        if(c(k).eq.1)then
                            res3(ig) = res3(ig)+u_tilde+v_tilde*ve2(k,1)
                        endif

                        if(stra(k).eq.1)then
                            res4(ig) = res4(ig) &
                            +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
                            res5(ig) = res5(ig) &
                            +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                            res6(ig) = res6(ig) &
                            +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
                        endif

                        if(stra(k).eq.2)then
                            res4(ig) = res4(ig) &
                            +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
                            res5(ig) = res5(ig) &
                            +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                            res6(ig) = res6(ig) &
                            +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
                        endif

                    endif
                end do

                som_Xbeta(ig) = vet

!=====fin maximisation aux
                integrale1(ig)= -0.5d0*dlog(sigma2*tau2) &
                -0.5d0* &
                dlog(res4(ig)*res6(ig)+res4(ig)/tau2+res6(ig)/sigma2 &
                +1.d0/(sigma2*tau2)-res5(ig)*res5(ig)) & !det(K"(b))
                +res2(ig)+res3(ig)-res4(ig) &
                -0.5d0*((u_tilde**2)/sigma2+(v_tilde**2)/tau2)!-k(b)
            end do
!======= fin calcul de log integrale
!======================================================================

            res = 0.d0
            do k=1,ng!ngexact

                if(nig(k).gt.0)then
                    if(indictronq.eq.0)then
                        res = res+integrale1(k)!integrale1 donne le log de I directement
                           if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                          funcpaasplines=-1.d9
                          goto 123
                           end if
                    else !troncature
                        indic_tronc = 1
                        funcpaasplines=-1.d9
                        goto 123
                !        write(*,*)'***TRAITER TRONCATURE**'
                    endif
                endif
            end do
        endif                !fin boucle effet= 1 and correl = 0

!=======================================================================

!************************************************************
!-----avec deux effets aleatoires dans le modele et correl=1
!************************************************************

        if(effet.eq.1.and.correl.eq.1)then

            mid=0
            integrale1=0.d0

            do k=1,nsujet
                if(c(k).eq.1)then
                    mid(g(k))=mid(g(k))+1
                endif
            end do

            do k=1,nsujet
                if(nva.gt.0)then
                    vet = 0.d0
                    do ip=1,nva
                        vet =vet + bh(np-nva+ip)*ve(k,ip)
                    end do
                    vet = dexp(vet)
                else
                    vet=1.d0
                endif

                if((c(k).eq.1).and.(stra(k).eq.1))then
                    res2(g(k)) = res2(g(k))+dlog(dut1(nt1(k))*vet)
                endif

                if((c(k).eq.1).and.(stra(k).eq.2))then
                    res2(g(k)) = res2(g(k))+dlog(dut2(nt1(k))*vet)
                endif

                if ((res2(g(k)).ne.res2(g(k))).or.(abs(res2(g(k))).ge. 1.d30)) then
                    funcpaasplines=-1.d9
                    goto 123
                end if

            end do

!=========================================================================
!==== calcul des integrales par transformation de LAPLACE pour chq gpe ===
!=========================================================================

            do ig=1,ngexact!ng!
                res3(ig)=0.d0
                res4(ig)=0.d0
                res5(ig)=0.d0
                res6(ig)=0.d0
                auxig=ig
                baux(1)=0.05d0   !initialisation de u_tilde
                baux(2)=0.05d0   !initialisation de v_tilde

!======================================================================
! maximisation pour deux parametres, les 2 effets aleatoires : u et v
!======================================================================
                npaux=2
                niaux=0
                vaux=0.d0!vecteur derivees 2nd et 1ERES
                funcaux=0.d0!vraisemblance
                resaux=0.d0!vraisemblance
                do ip=1,nva
                    betaaux(ip)= bh(np-nva+ip)
                end do

    call marq98o(baux,npaux,niaux,vaux,resaux,ieraux,istopaux,funcpao)
       if (ieraux .eq.-1) then
          funcpaasplines=-1.d9
          goto 123
       end if
                u_tilde = baux(1)!u_tilde(ig)
                v_tilde = baux(2)!v_tilde

                if(effet.eq.1) then
                end if
                do k=1,nsujet
                    if(nva.gt.0.and.g(k).eq.ig)then
                        vet = 0.d0
                        do ip=1,nva
                            vet =vet + bh(np-nva+ip)*ve(k,ip)
                        end do
                        vet = dexp(vet)
                    else
                        vet=1.d0
                    endif

                    if(g(k).eq.ig)then
                        if(c(k).eq.1)then
                            res3(ig) = res3(ig) &
                            +u_tilde+v_tilde*ve2(k,1)
                        endif

                        if(stra(k).eq.1)then
                            res4(ig) = res4(ig) &
                            +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
                            res5(ig) = res5(ig) &
                            +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                            res6(ig) = res6(ig) &
                            +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
                        endif

                        if(stra(k).eq.2)then
                            res4(ig) = res4(ig) &
                            +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
                            res5(ig) = res5(ig) &
                            +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                            res6(ig) = res6(ig) &
                            +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
                        endif

                    endif
                end do

                som_Xbeta(ig) = vet

!=====fin maximisation aux

                cov=dsqrt(sigma2*tau2)* & !avec contrainte
                (2.d0 * dexp(bh(np-nva-2))/(1.d0+dexp(bh(np-nva-2)))-1.d0)

                integrale1(ig)= -0.5d0*dlog(sigma2*tau2-cov*cov) &!-0.5logdet(D)
                -0.5d0* &
                dlog((res4(ig)+tau2/(tau2*sigma2-cov**2)) &
                *(res6(ig)+sigma2/(tau2*sigma2-cov**2)) &
                -(res5(ig)-cov/(tau2*sigma2-cov**2))**2) &!-0.5*log(det(ka2)
                +res2(ig)+res3(ig)-res4(ig) &
                -0.5d0*((u_tilde**2)/sigma2+(v_tilde**2)/tau2 &
                -2.d0*u_tilde*v_tilde*cov/(sigma2*tau2)) &
                /(1.d0-(cov**2)/(sigma2*tau2))         !-ka

            end do
!======= fin calcul de log integrale
!======================================================================

            res = 0.d0
            do k=1,ngexact  !ng!
                if(nig(k).gt.0)then
                    if(indictronq.eq.0)then
                        res = res &
                        +integrale1(k)!integrale1 donne le log de I directement
                       if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                          funcpaasplines=-1.d9
                          goto 123
                       end if
                    else !troncature
                        indic_tronc = 1
                        funcpaasplines=-1.d9
                        goto 123
                    !    write(*,*)'***TRAITER TRONCATURE**'
                    endif
                endif
            end do
        endif                     !fin boucle effet= 1 and correl = 1
    endif                     !fin boucle globale effet=0


!----------calcul de la penalisation -------------------
    pe=0.d0
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

        if (nst.eq.2) then
            pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
            *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
            the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
            m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
            the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
            m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
            *the2(i)*m1m(i))
        endif
    end do

    if (nst.eq.2) then
        pe = k0(1)*pe1 + k0(2)*pe2
    else
        pe = k0(1)*pe1
    endif

    resnonpen = res
    res = res - pe

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpaasplines=-1.d9
        goto 123
    end if

    funcpaasplines= res

123     continue

    return

    end function funcpaasplines


    subroutine distanceasplines(nz1,nz2,b,effet,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)

    use tailles
    use comon,only:zi,nva,nst,H_hess
    !use comon,only:c,date,hess,Hspl_hess,I_hess,ndate,nsujet,nt0,nt1,PEN_deri,t0,t1

    implicit none

    integer::nz1,nz2,i,j,n,np,k,l,effet
    double precision::x1,x2,su,bsup,binf,lam,lbinf, &
    h,lbsup
    double precision,dimension(npmax,npmax)::hes1,hes2
    double precision,dimension(-2:npmax):: the1,the2
    double precision,dimension(npmax):: b
    double precision,dimension(100)::x1Out,x2Out
    double precision,dimension(100,3)::lamOut,suOut,lam2Out,su2Out

    n  = nz1+2

    if(nst.eq.2)then
        np  = nz1+2+nz2+2+effet+nva
    else
        np  = nz1+2+effet+nva
    endif

    do i=1,nz1+2
        do j=1,nz1+2
            hes1(i,j)=h_Hess(i,j)
        end do
    end do

    if(nst.eq.2)then
        k = 0
        do i=nz1+3,nz1+2+nz2+2
            k = k + 1
            l = 0
            do j=nz1+3,nz1+2+nz2+2
                l = l + 1
                hes2(k,l)=H_hess(i,j)
            end do
        end do
    endif

    do i=1,nz1+2
        the1(i-3)=(b(i))*(b(i))
    end do

    if(nst.eq.2)then
        do i=1,nz2+2
            j = nz1+2+i
            the2(i-3)=(b(j))*(b(j))
        end do
    endif

    h = (zi(n)-zi(1))*0.01d0
    x1 = zi(1)
    x2 = zi(1)

    do i=1,100
        if(i .ne.1)then
            x1 = x1 + h
        end if
        call cospadd(x1,the1,nz1+2,hes1,zi,binf,su,bsup,lbinf,lam,lbsup)

        if(bsup.lt.0.d0)then
            bsup = 0.d0
        endif
        if(binf.gt.1.d0)then
            binf = 1.d0 
        endif
        if(lbinf.lt.0.d0)then
            lbinf = 0.d0
        endif

        x1Out(i)=x1
        lamOut(i,1)=lam
        lamOut(i,2)=lbinf
        lamOut(i,3)=lbsup
        suOut(i,1)=su
        suOut(i,2)=binf
        suOut(i,3)=bsup

        if(nst.eq.2)then
            if(i.ne.1)then
                x2 = x2 + h
            endif
            call cospadd(x2,the2,nz2+2,hes2,zi,binf,su,bsup,lbinf,lam,lbsup)
            if(bsup.lt.0.d0)then
                bsup = 0.d0
            endif
            if(binf.gt.1.d0)then
                binf = 1.d0 
            endif
            if(lbinf.lt.0.d0)then
                lbinf = 0.d0
            endif

            x2Out(i)=x2
            lam2Out(i,1)=lam
            lam2Out(i,2)=lbinf
            lam2Out(i,3)=lbsup
            su2Out(i,1)=su
            su2Out(i,2)=binf
            su2Out(i,3)=bsup

        endif
    end do

    return

    end subroutine distanceasplines


