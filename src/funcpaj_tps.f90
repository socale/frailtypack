
!========================          FUNCPAJ_CPM         ====================
    double precision function funcpaj_tps(b,np,id,thi,jd,thj,k0)

    use tailles
    !use comon,only:mm3,mm2,mm1,mm,im3,im2,im1,im,nva1,nva2,res4,ve,vedc
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m,nva, &
    t0,t1,t0dc,t1dc,c,cdc,nsujet, &
    nst,effet,ng,g,nig,kkapa,indic_alpha,alpha,theta, &
    auxig,aux1,aux2,res1,res3,typeof,pe,resnonpen,nb_gl
    use residusM
    use betatttps

    implicit none

    integer::np,id,jd,i,k,ig,choix,n,j
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,sum,res,pe1,pe2
    double precision,dimension(-2:npmax)::the1,the2
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res2,res1dc,res2dc
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3,integrale3gap
    double precision,dimension(ngmax)::integrale4
    double precision::logGammaJ,int
    double precision,dimension(2)::k0
    double precision::result,result3,abserr,resabs,resasc,resultdc
    double precision,external::risqindivrec,risqindivdc

    kkapa=k0

    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    if(effet.eq.1) then
        theta = bh(np-(nva+npbetatps)-indic_alpha)**(2.d0)
        if (indic_alpha.eq.1) then ! new : joint more flexible alpha = 1
            alpha = bh(np-(nva+npbetatps))
        else
            alpha = 1.d0
        endif
    endif

!-------------------------------------------------------
!---------- calcul de la vraisemblance ------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------cc

    res1 = 0.d0
    res2 = 0.d0
    res3 = 0.d0

    res1dc = 0.d0
    res2dc = 0.d0

    cpt = 0
    integrale1 = 0.d0
    integrale2 = 0.d0
    integrale3 = 0.d0

    integrale4 = 0.d0
    integrale3gap = 0.d0
    aux1 = 0.d0
    aux2 = 0.d0

!ccccccccccccccccccccccccccccccccccccccccc
!     pour les donnees recurrentes
!ccccccccccccccccccccccccccccccccccccccccc

    do i=1,nsujet

        cpt(g(i))=cpt(g(i))+1

!cccccccccccccccccccccc
! Fonction de risque de base recidive au temps T_ij
!cccccccccccccccccccccc
        if (c(i).eq.1) then
            res2(g(i)) = res2(g(i)) + dlog(risqindivrec(t1(i),i,bh,np))
        endif

        if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
            !print*,"here",res2(g(i)),risqindivrec(t1(i),i,bh,np)
            funcpaj_tps=-1.d9
            goto 123
        end if

!cccccccccccccccccccccc
! Fonction de risque cumulée de recidive au temps T_ij
!cccccccccccccccccccccc
        call integration(risqindivrec,0.d0,t1(i),result,abserr,resabs,resasc,i,bh,np)
        res1(g(i)) = res1(g(i)) + result

        if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
            !print*,"here2"
            funcpaj_tps=-1.d9
            goto 123
        end if

!cccccccccccccccccc
! Fonction de risque de recidive au tepms T_i(j-1)
!cccccccccccccccccc
        call integration(risqindivrec,0.d0,t0(i),result3,abserr,resabs,resasc,i,bh,np)
        res3(g(i)) = res3(g(i)) + result3

        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
            !print*,"here3"
            funcpaj_tps=-1.d9
            goto 123
        end if

    end do

!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces
!ccccccccccccccccccccccccccccccccccccccccc

    do k=1,ng
!cccccccccccccccccc
! Fonction de risque de deces au temps T_i*
!cccccccccccccccccc
        if (cdc(k).eq.1) then
            res2dc(k) = dlog(risqindivdc(t1dc(k),k,bh,np))
        endif

        if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
            !print*,"here4"
            funcpaj_tps=-1.d9
            goto 123
        end if

!cccccccccccccccccc
! Fonction de risque cumulée de dcd au temps T_i*
!cccccccccccccccccc
        call integration(risqindivdc,t0dc(k),t1dc(k),resultdc,abserr,resabs,resasc,k,bh,np)
        aux1(k) = resultdc

        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
            !print*,"here5"
            funcpaj_tps=-1.d9
            goto 123
        end if
    end do

!***************INTEGRALES ****************************
    do ig=1,ng
        auxig = ig
        choix = 3
        call gaulagj(int,choix,nb_gl)
        integrale3(ig) = int !moins bon
        if ((integrale3(ig).eq.0.d0).and.(typeof.ne.2)) then
            integrale3(ig) = 1.d-300
        endif
    end do
!************** FIN INTEGRALES ************************

    res = 0.d0
    do k=1,ng
        sum = 0.d0
        if(cpt(k).gt.0)then
            if(theta.gt.(1.d-8)) then
!cccc ancienne vraisemblance : pour calendar sans vrai troncature cccccccc
                res= res + res2(k) &
!--      pour le deces:
                + res2dc(k)  &
                - logGammaJ(1./theta)-dlog(theta)/theta  &
                + dlog(integrale3(k))
            else
!*************************************************************************
!     developpement de taylor d ordre 3
!*************************************************************************
!                   write(*,*)'************** TAYLOR *************'
                res= res + res2(k) &
                + res2dc(k) &
                - logGammaJ(1./theta)-dlog(theta)/theta  &
                + dlog(integrale3(k))

            endif
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                !print*,"here6",k,res,res2(k),res2dc(k),integrale3(k)
                funcpaj_tps=-1.d9
                goto 123
            end if
        endif
    end do

    if (typeof.eq.0) then ! penalisation pour les splines
        n = (np-(nva+npbetatps)-effet)/nst
        do k=1,n
            the1(k-3)=(bh(k))*(bh(k))
            j = n+k
            if (nst.eq.2) then
                the2(k-3)=(bh(j))*(bh(j))
            endif
        end do
        pe1 = 0.d0
        pe2 = 0.d0
        pe=0.d0
        do i=1,n-3
            pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
            *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
            the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
            m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
            the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
            m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
            *the1(i)*m1m(i))
            pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
            *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
            the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
            m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
            the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
            m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
            *the2(i)*m1m(i))
        end do

        if (nst.eq.1) then
            pe2=0.d0
        end if

        pe = k0(1)*pe1 + k0(2)*pe2

        resnonpen = res

        res = res - pe
    endif

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        !print*,"here7"
        funcpaj_tps =-1.d9
        goto 123
    endif

    funcpaj_tps = res

    ! pour les martingales
    Rrec = res1
    Nrec = nig
    Rdc = aux1
    Ndc = cdc

123     continue

    return

    end function funcpaj_tps


!====================================================================
!====================================================================
    double precision function risqindivrec(tps,i,bh,np)

    use tailles
    use comon
    use betatttps

    integer::p,j,jj,i,np,k,gg,n
    double precision::vet,tps
    double precision,dimension(-2:npmax)::the1,the2
    double precision::bbb,su
    double precision,dimension(np)::bh
    double precision::BasisSinhaT1(nbinnerknots+qorder)

    k=0
    j=0
    su=0.d0
    bbb=0.d0
    if(nva1.gt.0)then
        vet = 0.d0
        betatps=0.d0
        p=0 !p pointe le rang de chaque coefficient de regression
        j=0
        jj=0
        do j=1,nva1
            if (filtretps(j).eq.1)then
                call splinebasisIndiv(qorder-1,nbinnerknots+2*qorder,nbinnerknots,&
                nbinnerknots+qorder,tps,innerknots,boundaryknots,BasisSinhaT1)

                do jj=-qorder+1,nbinnerknots
                    betatps(j)=betatps(j)+bh(np-(nva+npbetatps)+p+jj+qorder)*BasisSinhaT1(jj+qorder)
                end do
            else
                betatps(j)=bh(np-(nva+npbetatps)+p+1)
            endif

            p=p+filtretps(j)*(nbinnerknots+qorder-1)+1
            vet = vet + betatps(j)*dble(ve(i,j))
        end do
        vet = dexp(vet)
    else
        vet = 1.d0
    endif

    select case(typeof)
        case(0) ! calcul du risque splines

        n = (np-(nva+npbetatps)-effet-indic_alpha)/2

        do k=1,n
            the1(k-3)=(bh(k))**2.d0
            the2(k-3)=(bh(n+k))**2.d0
        end do

        call susps(tps,the1,nzloco,su,bbb,zi)

        if (tps.eq.date(ndate)) then
            bbb = 4.d0*the1(n-2-1)/(zi(n-2)-zi(n-2-1))
        endif

        case(1) ! calcul du risque piecewise

        betacoef = 0.d0
        do k = 1,(nbintervR+nbintervDC)
            betacoef(k)=bh(k)**2
        end do

        do gg=1,nbintervR
            if((tps.ge.(ttt(gg-1))).and.(tps.lt.(ttt(gg))))then
                bbb = betacoef(gg)
            end if
        end do

        if((tps.ge.(ttt(nbintervR))))then
            bbb = betacoef(nbintervR)
        end if

        case(2) ! calcul du risque weibull

        betaR = bh(1)**2
        etaR = bh(2)**2

        if (tps.eq.0.d0) tps = 1d-12 ! utile car log(0) => -Inf

        bbb = (betaR*dexp((betaR-1.d0)*dlog(tps))/(etaR**betaR))

    end select

    risqindivrec = bbb*vet

    return

    end function risqindivrec

!====================================================================
    double precision function risqindivdc(tps,i,bh,np)

    use tailles
    use comon
    use betatttps

    integer::p,j,jj,i,np,k,gg,n
    double precision::vet2,tps
    double precision,dimension(-2:npmax)::the1,the2
    double precision::bbb,su
    double precision,dimension(np)::bh
    double precision::BasisSinhaT1(nbinnerknots+qorder)

    k=0
    j=0
    su=0.d0
    bbb=0.d0
    if(nva2.gt.0)then
        vet2 = 0.d0
        betatps2=0.d0
        p=0 !p pointe le rang de chaque coefficient de regression
        j=0
        jj=0
        do j=1,nva2
            if (filtre2tps(j).eq.1)then
                call splinebasisIndiv(qorder-1,nbinnerknots+2*qorder,nbinnerknots,&
                nbinnerknots+qorder,tps,innerknotsdc,boundaryknots,BasisSinhaT1)

                do jj=-qorder+1,nbinnerknots
                    betatps2(j)=betatps2(j) &
                    +bh(np-(nva+npbetatps)+(nva1+npbetatps1)+p+jj+qorder) &
                    *BasisSinhaT1(jj+qorder)
                end do
            else
                betatps2(j)=bh(np-(nva+npbetatps)+(nva1+npbetatps1)+p+1)
            endif

            p=p+filtre2tps(j)*(nbinnerknots+qorder-1)+1
            vet2 = vet2 + betatps2(j)*dble(vedc(i,j))
        end do
        vet2 = dexp(vet2)
    else
        vet2 = 1.d0
    endif

    select case(typeof)
        case(0) ! calcul du risque splines

        n = (np-(nva+npbetatps)-effet-indic_alpha)/2

        do k=1,n
            the1(k-3)=(bh(k))**2.d0
            the2(k-3)=(bh(n+k))**2.d0
        end do

        call susps(tps,the2,nzdc,su,bbb,zidc)

        if (tps.eq.datedc(ndatedc)) then
            bbb = 4.d0*the2(n-2-1)/(zidc(n-2)-zidc(n-2-1))
        endif

        case(1) ! calcul du risque piecewise

        betacoef = 0.d0
        do k = 1,(nbintervR+nbintervDC)
            betacoef(k)=bh(k)**2
        end do

        do gg=1,nbintervDC
            if ((tps.ge.(tttdc(gg-1))).and.(tps.lt.(tttdc(gg))))then
                bbb = betacoef(gg+nbintervR)
            end if
        end do

        if (tps.ge.(tttdc(nbintervR))) then
            bbb = betacoef(nbintervDC+nbintervR)
        end if

        case(2) ! calcul du risque weibull

        betaD = bh(3)**2
        etaD = bh(4)**2

        if (tps.eq.0.d0) tps = 1d-12 ! utile car log(0) => -Inf

        bbb = (betaD*dexp((betaD-1.d0)*dlog(tps))/(etaD**betaD))

    end select

    risqindivdc = bbb*vet2

    return

    end function risqindivdc

!====================================================================
!====================================================================
