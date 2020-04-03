
!--------------------------------------------------------------------
!                   Function CVPL for Joint frailty model
!--------------------------------------------------------------------


    subroutine cvpl(nobs,nsujet,groupe0,c0,cdc0,nva10,nva20,ve0,vedc0, &
    typeof0,nz0,zi0,ttt0,tttdc0,nbintervR0,nbintervDC0,np,b,H_1, &
    t00,t10,t0dc0,t1dc0,nt,valT,rl_cond,epoir,contribt,atrisk)

    use comon,only:typeof,nst,nbintervR,nbintervDC,nva,nva1,nva2,ve,vedc, &
    effet,nz1,nz2,zi,ttt,tttdc,c,cdc,date,ndate,datedc,ndatedc,t0,t1,t0dc,t1dc,g

    implicit none

    integer,intent(in)::np,nobs,nsujet,nt,nva10,nva20,typeof0,nz0,nbintervR0,nbintervDC0
    integer,dimension(nobs),intent(in)::c0,groupe0
    integer,dimension(nsujet),intent(in)::cdc0
    double precision,dimension(nobs,nva10),intent(in)::ve0
    double precision,dimension(nsujet,nva20),intent(in)::vedc0
    double precision,dimension(-2:nz0+3),intent(in)::zi0
    double precision,dimension(0:nbintervR0),intent(in)::ttt0
    double precision,dimension(0:nbintervDC0),intent(in)::tttdc0
    double precision,dimension(np),intent(in)::b
    double precision,dimension(np,np),intent(in)::H_1
    double precision,dimension(nobs),intent(in)::t00,t10
    double precision,dimension(nsujet),intent(in)::t0dc0,t1dc0
    double precision,dimension(nt),intent(in)::valT

    integer::i,k,t,nsujet_t
    double precision::rl_condt,trace3,min,max,mint,maxt,mindc,maxdc,mintdc,maxtdc
    double precision,dimension(0:(2*nobs))::aux
    double precision,dimension(nsujet)::rlindiv
    double precision,dimension(np,np)::J_cond,mat3
    integer,dimension(nsujet)::indT

    double precision,dimension(nt),intent(out)::rl_cond,epoir
    double precision,dimension(nt*nsujet),intent(out)::contribt
    double precision,dimension(nt),intent(out)::atrisk

    ! redefinition des variables du module
    nst = 2
    nva1 = nva10
    nva2 = nva20
    nva = nva1+nva2
    nbintervR = nbintervR0
    nbintervDC = nbintervDC0
    effet = 1
    typeof = typeof0

    nz1 = nz0
    nz2 = nz0

    allocate(c(nobs),cdc(nsujet),ve(nobs,nva1),vedc(nsujet,nva2),zi(-2:nz1+3),ttt(0:nbintervR),tttdc(0:nbintervDC))
    c = c0
    cdc = cdc0
    ve = ve0
    vedc = vedc0
    zi = zi0
    ttt = ttt0
    tttdc = tttdc0

    allocate(t0(nobs),t1(nobs),t0dc(nsujet),t1dc(nsujet),g(nobs))
    t0 = t00
    t1 = t10
    t0dc = t0dc0
    t1dc = t1dc0
    g = groupe0

    allocate(date(0:(2*nobs)),datedc(0:(2*nsujet)))
    date = 0.d0

    maxt = 0.d0
    mint = 0.d0

    do i=1,nobs
        if (i.eq.1) then
            mint = t0(i) ! affectation du min juste une fois
        endif
        if (maxt.lt.t1(i)) then
            maxt = t1(i)
        endif
!         if ((maxt.lt.tU(i)).and.(tU(i).ne.t1(i))) then
!             maxt = tU(i)
!         endif
        if (mint.gt.t0(i)) then
            mint = t0(i)
        endif
    end do

    maxtdc = 0.d0
    mintdc = 0.d0
    contribt = 0.d0

    do i=1,nsujet
        if (i.eq.1) then
            mintdc = t0dc(i) ! affectation du min juste une fois
        endif
        if (maxtdc.lt.t1dc(i)) then
            maxtdc = t1dc(i)
        endif
        if (mintdc.gt.t0dc(i)) then
            mintdc = t0dc(i)
        endif
    end do

    min = 1.d-10
    max = maxt
    aux = 0.d0

    do i = 1,(2*nobs)
        do k = 1,nobs
            if (t0(k).ge.min) then
                if(t0(k).lt.max)then
                    max = t0(k)
                endif
            endif
            if (t1(k).ge.min) then
                if(t1(k).lt.max)then
                    max = t1(k)
                endif
            endif
!             if (tU(k).ne.t1(k)) then
!                 if((tU(k).ge.min))then
!                     if (tU(k).lt.max) then
!                         max = tU(k)
!                     endif
!                 endif
!             endif
        end do
        aux(i) = max
        min = max + 1.d-12 ! pour virer les doublons
        max = maxt
    end do

    date(1) = aux(1)
    k = 1
    do i=2,(2*nobs)
        if (aux(i).gt.aux(i-1)) then
            k = k+1
            date(k) = aux(i)
        endif
    end do
    ndate = k

    mindc = 0.d0
    maxdc = maxtdc
    aux = 0.d0

    do i = 1,(2*nsujet)
        do k = 1,nsujet
            if (t0dc(k).ge.mindc) then
                if (t0dc(k).lt.maxdc) then
                    maxdc = t0dc(k)
                endif
            endif
            if (t1dc(k).ge.mindc) then
                if (t1dc(k).lt.maxdc) then
                    maxdc = t1dc(k)
                endif
            endif
        end do
        aux(i) = maxdc
        mindc = maxdc + 1.d-12
        maxdc = maxtdc
    end do

    datedc(1) = aux(1)
    k = 1
    do i=2,(2*nsujet)
        if (aux(i).gt.aux(i-1)) then
            k = k+1
            datedc(k) = aux(i)
        endif
    end do
    ndatedc = k
    ! fin de la redefinition des variables du module

    do t=1,nt ! boucle sur les temps de validation

        indT = 0
        nsujet_t = 0
        do i=1,nsujet ! boucle sur les individus
            if (t1dc(i).ge.valT(t)) then
                indT(i) = 1
                nsujet_t = nsujet_t + 1
            end if
        end do

        atrisk(t) = nsujet_t
        J_cond = 0.d0
        rlindiv = 0.d0

        call derivc_condT(b,np,J_cond,rlindiv,nobs,nsujet,indT,valT(t))

        rl_condt = 0.d0
        do i=1,nsujet
            contribt(nsujet*(t-1)+i)=rlindiv(i)
!             write(*,*)rlindiv(i)
            if (rlindiv(i).eq.-1.d9) then
!                 print*,"oups"
                rl_cond(t) = -1.d9
                epoir(t) = -1.d9
                goto 5289
            end if
            rl_condt = rl_condt + rlindiv(i)
        end do

        mat3 = MATMUL(H_1,J_cond)
        trace3 = 0.d0
        do k=1,np
            trace3 = trace3 + mat3(k,k)
        end do

        epoir(t) = -rl_condt/dble(nsujet_t)+(trace3*dble(nsujet)/(dble(nsujet_t)*dble(nsujet-1))) !cvpl
        rl_cond(t) = -rl_condt/dble(nsujet_t) !mpl

        if (epoir(t).ne.epoir(t)) then
            epoir(t) = -1.d9
        end if
        if (rl_cond(t).ne.rl_cond(t)) then
            rl_cond(t) = -1.d9
        end if

5289    continue

    end do

    deallocate(c,cdc,ve,vedc,zi,ttt,tttdc,date,datedc)
    deallocate(t0,t1,t0dc,t1dc,g)

    end subroutine cvpl


!-----------------------------------------------------------
!                        derivc_condt
!------------------------------------------------------------
    subroutine derivc_condt(b,m,V,rlindiv,nobs,nsujet,indT,valT)

    implicit none

    integer::m,i,k,id,nsujet,nobs
    double precision::funcpi,thn,th,z,temp1,temp2,valT
    double precision,dimension(m,1)::Uscore,Uscore2
    double precision,dimension(m)::b
    double precision,dimension(m,m)::V
    double precision,dimension(nsujet)::rlindiv
    
    integer,dimension(nsujet)::indT

    V = 0.d0
    rlindiv = 0.d0
    z = 0.d0
    id = 0

    do i=1,nsujet
        Uscore = 0.d0
        Uscore2 = 0.d0
        if (indT(i).eq.1) then ! contribution de i sachant T
            rlindiv(i) = funcpi(nobs,b,m,id,z,id,z,i,1,valT)
!             write(*,*)nsujet,rlindiv(i)
            if (rlindiv(i).eq.-1.d9) then 
                V = 0.d0
                rlindiv = -1.d9
                goto 777
            end if
        end if

        do k=1,m
            th = 1.d-6
            thn = -1.D0*th
            if (indT(i).eq.1) then ! calcul des derivees sachant T
                temp1 = funcpi(nobs,b,m,k,th,id,z,i,1,valT)
                temp2 = funcpi(nobs,b,m,k,thn,id,z,i,1,valT)
                if (temp1.eq.-1.d9.or.temp2.eq.-1.d9) then 
                    V = 0.d0
                    !rlindiv = -1.d9
                    goto 777
                end if
                Uscore(k,1) = -(temp1-temp2)/(2.d0*th)
            end if
            ! calcul des derivees
            temp1 = funcpi(nobs,b,m,k,th,id,z,i,2,valT)
            temp2 = funcpi(nobs,b,m,k,thn,id,z,i,2,valT)
            if (temp1.eq.-1.d9.or.temp2.eq.-1.d9) then 
                V = 0.d0
                !rlindiv = -1.d9
                goto 777
            end if
            Uscore2(k,1) = -(temp1-temp2)/(2.d0*th)
        end do
        V = V + MATMUL(Uscore,transpose(Uscore2))
    end do

777 continue

    return

    end subroutine derivc_condt


!-----------------------------------------------------------
!                        FUNCPI
!------------------------------------------------------------

    double precision function funcpi(nobs,b,np,id,thi,jd,thj,i,choix,valT)
    
    !use comon,only:typeof,stra,c,cdc,nva1,nva2,nva,nst,nbintervR,nbintervDC,ve,effet,nz1,nz2, &
    !zi,ttt,date,datedc,ndate,ndatedc,vedc,t0,t1,t1dc,tttdc,g

    implicit none

    integer::np,id,jd,i,nobs,choix
    double precision::thi,thj,valT
    double precision,dimension(np)::b,bh
    double precision::integrale1,integrale2

    bh = b
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    if (choix.eq.1) then
        integrale1 = 0.d0
        call gaulagEpoce(integrale1,1,bh,np,i,nobs,valT,0)
        integrale2 = 0.d0
        call gaulagEpoce(integrale2,2,bh,np,i,nobs,valT,0)

        funcpi = integrale1/integrale2 ! conditionnellement aux reccurences
    else
        integrale1 = 0.d0
        call gaulagEpoce(integrale1,1,bh,np,i,nobs,valT,1)

        funcpi = integrale1
    endif

    funcpi = dlog(funcpi)

    if ((funcpi.ne.funcpi).or.(abs(funcpi).gt.1.d30)) then
        funcpi = -1.d9
    end if

    return

    end function funcpi

!------------------------------------------------------------
!------------------------------------------------------------

    subroutine gaulagEpoce(ss,choix,bh,np,i,nobs,valT,all)

    use donnees,only:w,x

    implicit none

    integer,intent(in)::np,choix,i,nobs,all
    double precision,intent(in)::valT
    double precision,dimension(np),intent(in)::bh
    double precision ::auxfunca,func1E,func2E
    external::func1E,func2E
    integer::j
    double precision,intent(out)::ss

    ss = 0.d0
    do j=1,20
        if (choix.eq.1) then
            auxfunca = func1E(x(j),bh,np,i,nobs,valT,all)
            ss = ss+w(j)*(auxfunca)
        else
            if (choix.eq.2) then 
                auxfunca = func2E(x(j),bh,np,i,nobs,valT)
                ss = ss+w(j)*(auxfunca)
            endif
        endif
    end do

    return

    end subroutine gaulagEpoce

!==================================================================

    double precision function func1E(frail,bh,np,i,nobs,valT,all) ! integrale 1 au numerateur
    ! recurrences + deces
    !use comon,only:stra
    use comon,only:typeof,nst,nbintervR,nbintervDC,nva,ve,effet,nz1,nz2, &
    zi,ttt,c,date,datedc,ndate,ndatedc,nva1,nva2,vedc,t0,t1,t1dc,tttdc,g,cdc

    IMPLICIT NONE

    double precision,intent(in)::frail,valT
    integer,intent(in)::i,nobs,np,all
    double precision,dimension(np),intent(in)::bh
    integer::n,k,j,gg
    double precision,dimension(-2:np)::the1,the2
    double precision,dimension(np)::betacoef
    double precision::betaR,etaR,betaD,etaD
    double precision::vet,vet2,alpha,theta
    double precision,dimension(2)::su,sut1,sut0,sudc
    double precision::lam,lamdc,temp, tempscl
    double precision::logGammaJ

    n = 0
    betaR = 0.d0
    etaR = 0.d0
    betaD = 0.d0
    etaD = 0.d0

    theta = bh(np-nva-1)*bh(np-nva-1)
    alpha = bh(np-nva)

    select case(typeof)
        case(0)
            n = (np-nva-effet-1)/nst
            do k=1,n
                the1(k-3) = (bh(k))*(bh(k))
                j = n+k
                the2(k-3) = (bh(j))*(bh(j))
            end do
        case(1)
            betacoef = 0.d0
            do k=1,(nbintervR+nbintervDC)
                betacoef(k) = bh(k)**2
            end do
        case(2)
            betaR = bh(1)**2
            etaR = bh(2)**2
            betaD = bh(3)**2
            etaD = bh(4)**2
    end select

    func1E = 1.d0

! ------------ Reccurent ------------- !

    do k=1,nobs ! toutes les recurrences pour l'individu i
        if (g(k).eq.i) then

            if ((t1(k).gt.valT).and.(all.eq.0)) then ! si on veut filtrer les recurrences (all = 0)
                goto 7                               ! et si la recurrence est apres t, on sort de la boucle
            endif

            if(nva1.gt.0)then
                vet = 0.d0
                do j=1,nva1
                    vet = vet + bh(np-nva+j)*dble(ve(k,j))
                end do
                vet = dexp(vet)
            else
                vet = 1.d0
            endif

            ! CALCUL DE LA SURVIE
            select case(typeof)
                case(0)
                      call susps(t1(k),the1,nz1,temp,lam,zi)
                      sut1(1) = temp
                      call susps(t0(k),the1,nz1,temp,lam,zi)
                      sut0(1) = temp
                case(1)
                    call survival_cpm(t1(k),bh,nst,nbintervR,ttt,sut1)
                    call survival_cpm(t0(k),bh,nst,nbintervR,ttt,sut0)
               case(2)
                    sut1(1) = dexp(-(t1(k)/etaR)**betaR)
                    sut0(1) = dexp(-(t0(k)/etaR)**betaR)
            end select

            func1E = func1E * (sut1(1)/sut0(1))**(frail*vet)

            if ((func1E.ne.func1E).or.(abs(func1E).gt.1.d30)) then
                func1E = -1.d9
!                 print*,"1"
                goto 1000
            end if

            ! CALCUL DU RISQUE
            if (c(k).eq.1) then
                select case(typeof)
                    case(0)
                        call susps(t1(k),the1,nz1,tempscl,lam,zi)
                        su = tempscl
                        if (t1(k).eq.date(ndate)) then
                            lam = 4.d0*the1(n-2-1)/(zi(n-2)-zi(n-2-1))
                        endif
                    case(1)
                        do gg=1,nbintervR
                            if ((t1(k).ge.ttt(gg-1)).and.(t1(k).lt.ttt(gg))) then
                                lam = betacoef(gg)
                            end if
                        end do
                        if (t1(k).ge.ttt(nbintervR)) then
                            lam = betacoef(nbintervR)
                        end if
                    case(2)
                        if (t1(k).eq.0.d0) t1(k) = 1d-12 ! utile car log(0) => -Inf
                        ! ecriture en exp(log) pour virer l'exposant
                        lam = (betaR*dexp((betaR-1.d0)*dlog(t1(k)))/(etaR**betaR))
                end select

                func1E = func1E * frail*lam*vet

                if ((func1E.ne.func1E).or.(abs(func1E).gt.1.d30)) then
                    func1E = -1.d9
!                     print*,"2"
                    goto 1000
                end if
            endif
        endif

7   continue

    end do

! ------------ Death ------------- !

    if(nva2.gt.0)then
        vet2 = 0.d0
        do j=1,nva2
            vet2 = vet2 + bh(np-nva2+j)*dble(vedc(i,j))
        end do
        vet2 = dexp(vet2)
    else
        vet2 = 1.d0
    endif

    ! CALCUL DE LA SURVIE
    select case(typeof)
        case(0)
              call susps(t1dc(i),the2,nz2,temp,lamdc,zi)
              sudc(2) = temp
        case(1)
            call survivalj_cpm(t1dc(i),bh,nbintervR,nbintervDC,ttt,tttdc,sudc)
        case(2)
            sudc(2) = dexp(-(t1dc(i)/etaD)**betaD)
    end select

    func1E = func1E * sudc(2)**(frail**alpha * vet2)

    if ((func1E.ne.func1E).or.(abs(func1E).gt.1.d30)) then
        func1E = -1.d9
!         print*,"3",func1E,i,t1dc(i),the2,nz2,zi,nst
        goto 1000
    end if

    ! CALCUL DU RISQUE
    if (cdc(i).eq.1) then
        select case(typeof)
            case(0)
                call susps(t1dc(i),the2,nz2,tempscl,lamdc,zi)
                sudc = tempscl
                if (t1dc(i).eq.datedc(ndatedc)) then
                    lamdc = 4.d0*the2(n-2-1)/(zi(n-2)-zi(n-2-1))
                endif
            case(1)
                do gg=1,nbintervDC
                    if ((t1dc(i).ge.tttdc(gg-1)).and.(t1dc(i).lt.tttdc(gg))) then
                        lamdc = betacoef(gg+nbintervR)
                    end if
                end do
                if (t1dc(i).ge.tttdc(nbintervDC)) then
                    lamdc = betacoef(nbintervDC+nbintervR)
                end if
            case(2)
                if (t1dc(i).eq.0.d0) t1dc(i) = 1d-12
                lamdc = (betaD*dexp((betaD-1.d0)*dlog(t1dc(i)))/(etaD**betaD))
        end select

        func1E = func1E * frail**alpha * lamdc*vet2

        if ((func1E.ne.func1E).or.(abs(func1E).gt.1.d30)) then
            func1E = -1.d9
!             print*,"4"
            goto 1000
        end if
    endif

    ! densite de la loi gamma pour les effets aleatoires
    func1E = func1E * (frail**(1.d0/theta-1.d0)*dexp(-frail/theta))/(dexp(logGammaJ(1.d0/theta))*theta**(1.d0/theta))

1000 continue

    return

    end function func1E

!==================================================================

    double precision function func2E(frail,bh,np,i,nobs,valT) ! integrale 2 au denominateur
    ! recurrences

    !use comon,only:stra,datedc,ndatedc,t1dc,cdc
    use comon,only:typeof,nst,nbintervR,nbintervDC,nva,ve,effet,nz1,nz2, &
    zi,ttt,c,date,ndate,nva1,nva2,vedc,t0,t1,tttdc,g

    IMPLICIT NONE

    double precision,intent(in)::frail,valT
    integer,intent(in)::i,nobs,np
    double precision,dimension(np),intent(in)::bh
    integer::n,k,j,gg
    double precision,dimension(-2:np)::the1,the2
    double precision,dimension(np)::betacoef
    double precision::betaR,etaR,betaD,etaD
    double precision::vet,vet2,alpha,theta
    double precision,dimension(2)::su,sut1,sut0,sudc
    double precision::lam,lamdc,temp, tempscl
    double precision::logGammaJ

    n = 0
    betaR = 0.d0
    etaR = 0.d0
    betaD = 0.d0
    etaD = 0.d0

    theta = bh(np-nva-1)*bh(np-nva-1)
    alpha = bh(np-nva)

    select case(typeof)
        case(0)
            n = (np-nva-effet-1)/nst
            do k=1,n
                the1(k-3) = (bh(k))*(bh(k))
                j = n+k
                the2(k-3) = (bh(j))*(bh(j))
            end do
        case(1)
            betacoef = 0.d0
            do k=1,(nbintervR+nbintervDC)
                betacoef(k) = bh(k)**2
            end do
        case(2)
            betaR = bh(1)**2
            etaR = bh(2)**2
            betaD = bh(3)**2
            etaD = bh(4)**2
    end select

    func2E = 1.d0

! ------------ Reccurent ------------- !

    do k=1,nobs ! toutes les recurrences pour l'individu i
        if ((g(k).eq.i).and.(t1(k).le.valT)) then

            if(nva1.gt.0)then
                vet = 0.d0
                do j=1,nva1
                    vet = vet + bh(np-nva+j)*dble(ve(k,j))
                end do
                vet = dexp(vet)
            else
                vet = 1.d0
            endif

            ! CALCUL DE LA SURVIE
            select case(typeof)
                case(0)
                      call susps(t1(k),the1,nz1,temp,lam,zi)
                      sut1(1) = temp
                      call susps(t0(k),the1,nz1,temp,lam,zi)
                      sut0(1) = temp
                case(1)
                    call survival_cpm(t1(k),bh,nst,nbintervR,ttt,sut1)
                    call survival_cpm(t0(k),bh,nst,nbintervR,ttt,sut0)
               case(2)
                    sut1(1) = dexp(-(t1(k)/etaR)**betaR)
                    sut0(1) = dexp(-(t0(k)/etaR)**betaR)
            end select

            func2E = func2E * (sut1(1)/sut0(1))**(frail*vet)

            if ((func2E.ne.func2E).or.(abs(func2E).gt.1.d30)) then
!                 print*,"1",func2E
                func2E = -1.d9
                goto 1000
            end if

            ! CALCUL DU RISQUE
            if (c(k).eq.1) then
                select case(typeof)
                    case(0)
                        call susps(t1(k),the1,nz1,tempscl,lam,zi)
                            su = tempscl
                        if (t1(k).eq.date(ndate)) then
                            lam = 4.d0*the1(n-2-1)/(zi(n-2)-zi(n-2-1))
                        endif
                    case(1)
                        do gg=1,nbintervR
                            if ((t1(k).ge.ttt(gg-1)).and.(t1(k).lt.ttt(gg))) then
                                lam = betacoef(gg)
                            end if
                        end do
                        if (t1(k).ge.ttt(nbintervR)) then
                            lam = betacoef(nbintervR)
                        end if
                    case(2)
                        if (t1(k).eq.0.d0) t1(k) = 1d-12 ! utile car log(0) => -Inf
                        ! ecriture en exp(log) pour virer l'exposant
                        lam = (betaR*dexp((betaR-1.d0)*dlog(t1(k)))/(etaR**betaR))
                end select

                func2E = func2E * frail*lam*vet

                if ((func2E.ne.func2E).or.(abs(func2E).gt.1.d30)) then
                    func2E = -1.d9
!                     print*,"2"
                    goto 1000
                end if
            endif
        endif
    end do

! ------------ Death ------------- !

    if(nva2.gt.0)then
        vet2 = 0.d0
        do j=1,nva2
            vet2 = vet2 + bh(np-nva2+j)*dble(vedc(i,j))
        end do
        vet2 = dexp(vet2)
    else
        vet2 = 1.d0
    endif

    ! CALCUL DE LA SURVIE
    select case(typeof)
        case(0)
              call susps(valT,the2,nz2,temp,lamdc,zi)
              sudc(2) = temp
        case(1)
            call survivalj_cpm(valT,bh,nbintervR,nbintervDC,ttt,tttdc,sudc)
        case(2)
            sudc(2) = dexp(-(valT/etaD)**betaD)
    end select

    func2E = func2E * sudc(2)**(frail**alpha * vet2)

    if ((func2E.ne.func2E).or.(abs(func2E).gt.1.d30)) then
        func2E = -1.d9
!         print*,"3",func2E,i,valT,the2,nz2,zi,nst
        goto 1000
    end if

    ! densite de la loi gamma pour les effets aleatoires
    func2E = func2E * (frail**(1.d0/theta-1.d0)*dexp(-frail/theta))/(dexp(logGammaJ(1.d0/theta))*theta**(1.d0/theta))

1000 continue

    return

    end function func2E


