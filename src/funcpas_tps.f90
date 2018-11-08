

!========================          FUNCPA_CPM_TPS          ====================
    double precision function funcpas_tps(b,np,id,thi,jd,thj,k0)

    use tailles
    !use comon,only:im3,im2,im1,im,mm3,mm2,mm1,mm,stra,ve,
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
    t0,t1,c,nsujet,nva, &
    nst,effet,ng,g,nig,AG,kkapa,theta,typeof,pe,resnonpen
    !use betatttps,only:betatps3,betatps,qorder,filtretps,
    use betatttps,only:nbinnerknots,knotsTPS,npbetatps,&
    innerknots,boundaryknots
    use residusM

    implicit none

! *** NOUVELLLE DECLARATION F90 :

    integer::nb,np,id,jd,i,k,cptg,l,n,j
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,dnb,sum,inv,res,pe1,pe2
    double precision,dimension(-2:npmax)::the1,the2
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res1,res2,res3
    double precision,dimension(2)::k0
    double precision::result,result3,abserr,resabs,resasc
    double precision,external::risqindiv

    kkapa=k0

    do i=1,np
        bh(i)=b(i)
    end do

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    if (effet.eq.1) then
        theta = bh(np-(nva+npbetatps))*bh(np-(nva+npbetatps))
    endif

!-------------------------------------------------------
!--------- calcul de la vraisemblance ------------------
!--------------------------------------------------------

!--- avec ou sans variable explicative  ------cc

    res1 = 0.d0
    res2 = 0.d0
    res3 = 0.d0
    cpt = 0

! Calcul des bases de splines pour tout i
     innerknots(1:nbinnerknots)=knotsTPS(1:nbinnerknots)
     boundaryknots(1)=knotsTPS(0)
     boundaryknots(2)=knotsTPS(nbinnerknots+1)

!*******************************************
!---- sans effet aleatoire dans le modele
!*******************************************

    if (effet.eq.0) then
        do i=1,nsujet

            cpt(g(i))=cpt(g(i))+1

!cccccccccccccccccccccc
! Fonction de risque de base recidive au temps T_ij
!cccccccccccccccccccccc
            if(c(i).eq.1)then
                res2(g(i)) = res2(g(i)) + dlog(risqindiv(t1(i),i,bh,np))
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                !print*,"cox here"
                funcpas_tps=-1.d9
                goto 123
            end if

!!cccccccccccccccccccccc
!! Fonction de risque cumulée de recidive au temps T_ij
!!cccccccccccccccccccccc
            call integration(risqindiv,0.d0,t1(i),result,abserr,resabs,resasc,i,bh,np)
            res1(g(i)) = res1(g(i)) + result

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                !print*,"cox here2"
                funcpas_tps=-1.d9
                goto 123
            end if
        end do

        res = 0.d0
        cptg = 0

! k indice les groupes
        do k=1,ng
            if(cpt(k).gt.0)then
                nb = nig(k)
                dnb = dble(nig(k))
                res = res-res1(k)+res2(k)
                cptg = cptg + 1
            endif
        end do

!*********************************************
!----avec un effet aleatoire dans le modele
!*********************************************

    else
!      write(*,*)'AVEC EFFET ALEATOIRE'
        inv = 1.d0/theta

        do i=1,nsujet

            cpt(g(i))=cpt(g(i))+1

!cccccccccccccccccccccc
! Fonction de risque de base recidive au temps T_ij
!cccccccccccccccccccccc
            if(c(i).eq.1)then
                !print*,t1(i),risqindiv(t1(i),i,bh,np)
                res2(g(i)) = res2(g(i)) + dlog(risqindiv(t1(i),i,bh,np))
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                !print*,t1(i),risqindiv(t1(i),i,bh,np)
                !print*,"here",res2(g(i)),i
                funcpas_tps=-1.d9
                goto 123
            end if
!!cccccccccccccccccccccc
!! Fonction de risque cumulée de recidive au temps T_ij
!!cccccccccccccccccccccc
            call integration(risqindiv,0.d0,t1(i),result,abserr,resabs,resasc,i,bh,np)
            !print*,"result",result,"res1",res1(g(i))
            res1(g(i)) = res1(g(i)) + result

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                !print*,"here2",res1(g(i)),i
                funcpas_tps=-1.d9
                goto 123
            end if

!!cccccccccccccccccccccc
!! Fonction de risque cumulée de recidive au temps T_i(j-1)
!!cccccccccccccccccccccc
            call integration(risqindiv,0.d0,t0(i),result3,abserr,resabs,resasc,i,bh,np)
            !print*,"result3",result3,"res3",res3(g(i))
            res3(g(i)) = res3(g(i)) + result3

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                !print*,"here3",res3(g(i)),i,g(i)
                funcpas_tps=-1.d9
                goto 123
            end if

            if (res1(g(i)).lt.res3(g(i))) then
                !print*,"here33",res1(g(i)),res3(g(i))
                funcpas_tps=-1.d9
                goto 123
            end if

        end do

        res = 0.d0
        cptg = 0

!     gam2 = gamma(inv)
! k indice les groupes
!!cccccccccccccccccccccc
!! CALCUL DE LOG-VRAISEMBLANCE
!!cccccccccccccccccccccc
        do k=1,ng
            sum=0.d0
            if(cpt(k).gt.0)then
                nb = nig(k)
                dnb = dble(nig(k))
                if (dnb.gt.1.d0) then
                    do l=1,nb
                        sum=sum+dlog(1.d0+theta*dble(nb-l))
                    end do
                endif
                if(theta.gt.(1.d-5)) then
!ccccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res= res-(inv+dnb)*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        + res2(k) + sum
!ccccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res-(inv+dnb)*dlog(theta*res1(k)+1.d0)  &
                        +(inv)*dlog(theta*res3(k)+1.d0)+ res2(k) + sum
                    endif
                else
!     developpement de taylor d ordre 3
!                   write(*,*)'************** TAYLOR *************'
!cccc ancienne vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res = res-dnb*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        -(res1(k)-res3(k))*(1.d0-theta*(res1(k)-res3(k))/2.d0 &
                        +theta*theta*(res1(k)-res3(k))*(res1(k)-res3(k))/3.d0)+res2(k)+sum
!cccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res-dnb*dlog(theta*res1(k)+1.d0)-res1(k)*(1.d0-theta*res1(k) &
                        /2.d0+theta*theta*res1(k)*res1(k)/3.d0) &
                        +res2(k)+sum &
                        +res3(k)*(1.d0-theta*res3(k)/2.d0 &
                        +theta*theta*res3(k)*res3(k)/3.d0)
                    endif
                endif
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    !print*,"here4",res
                    funcpas_tps=-1.d9
                    goto 123
                end if
            endif
        end do
    endif !fin boucle effet=0

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
        !print*,"here5",res
        funcpas_tps=-1.d9
        goto 123
    end if

    funcpas_tps = res

    ! pour les martingales
    cumulhaz = res1

123     continue

    return

    end function funcpas_tps

!====================================================================
!====================================================================
    double precision function risqindiv(tps,i,bh,np)

! t=temps d'event, i=rang de l'observation

    use tailles
    use comon,only:date,zi,nz1,nz2,nva,ve,ndate,stra,nst,effet,typeof, &
    nbintervR,ttt,betacoef,etaR,etaD,betaR,betaD
    !use betatttps,only:betatps3,knotsTPS
    use betatttps,only:filtretps,nbinnerknots,qorder,npbetatps,betatps,&
    innerknots,boundaryknots

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
    if(nva.gt.0)then
        vet = 0.d0
        betatps=0.d0
        p=0 !p pointe le rang de chaque coefficient de regression
        j=0
        jj=0
        do j=1,nva
!================== Approche Sinha
            if (filtretps(j).eq.1) then
!====== Calcul des bases de Splines
                call splinebasisIndiv(qorder-1,nbinnerknots+2*qorder,nbinnerknots, &
                nbinnerknots+qorder,tps,innerknots,boundaryknots,BasisSinhaT1)
                do jj=-qorder+1,nbinnerknots
!====== Calcul beta(t1) et beta(t0)
                    betatps(j)=betatps(j)+bh(np-(nva+npbetatps)+p+jj+qorder)&
                    *BasisSinhaT1(jj+qorder)
                end do
!================== fin Approche Sinha
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

        n = (np-(nva+npbetatps)-effet)/nst

        do k=1,n
            the1(k-3)=(bh(k))*(bh(k))
            j = n+k
            if (nst.eq.2) then
                the2(k-3)=(bh(j))*(bh(j))
            endif
        end do

        if (stra(i).eq.1) then
!============ fonction de risque
            call susps(tps,the1,nz1,su,bbb,zi)

! le risque du temps maximum n'est pas calculé dans la fonction susps
            if (tps.eq.date(ndate)) then
                bbb = 4.d0*the1(n-2-1)/(zi(n-2)-zi(n-2-1))
            endif
!============ fonction de risque
        endif

        if (stra(i).eq.2) then
            call susps(tps,the2,nz2,su,bbb,zi)
            if (tps.eq.date(ndate)) then
                bbb = 4.d0*the2(n-2-1)/(zi(n-2)-zi(n-2-1))
            endif
        endif

        case(1) ! calcul du risque piecewise

        betacoef = 0.d0
        do k=1,nst*nbintervR
            betacoef(k)=bh(k)**2
        end do

        if (stra(i).eq.1) then
            do gg=1,nbintervR
                if((tps.ge.(ttt(gg-1))).and.(tps.lt.(ttt(gg))))then
                    bbb = betacoef(gg)
                end if
            end do
            if((tps.ge.(ttt(nbintervR))))then
                bbb = betacoef(nbintervR)
            end if
        endif

        if (stra(i).eq.2) then
            do gg=1,nbintervR
                if((tps.ge.(ttt(gg-1))).and.(tps.lt.(ttt(gg))))then
                    bbb = betacoef(gg+nbintervR)
                end if
            end do
            if((tps.ge.(ttt(nbintervR))))then
                bbb = betacoef(nbintervR+nbintervR)
            end if
        endif

        case(2) ! calcul du risque weibull

        if (nst.eq.1) then
            betaR = bh(1)**2
            etaR = bh(2)**2
            etaD = 0.d0
            betaD = 0.d0
        else
            betaR = bh(1)**2
            etaR = bh(2)**2
            betaD = bh(3)**2
            etaD = bh(4)**2
        end if

        if (tps.eq.0.d0) tps = 1d-12 ! utile car log(0) => -Inf

        if (stra(i).eq.1) then
            ! ecriture en exp(log) pour virer l'exposant
            bbb = (betaR*dexp((betaR-1.d0)*dlog(tps))/(etaR**betaR))
        endif

        if (stra(i).eq.2) then
            bbb = (betaD*dexp((betaD-1.d0)*dlog(tps))/(etaD**betaD))
        endif

    end select

    risqindiv = bbb*vet

    return

    end function risqindiv


