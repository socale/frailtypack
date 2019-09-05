


!========================          funcpajweib_intcens       ====================
    double precision function funcpajweib_intcens(b,np,id,thi,jd,thj,k0)

    !use comon,only:AG,indictronqdc,nz1,nz2,stra
    use comon,only:etaR,etaD,betaR,betaD,&
    t0,t1,t0dc,t1dc,c,cdc,nsujet,nva,nva1,nva2,nst,&
    ve,vedc,effet,ng,g,nig,indic_ALPHA,theta,alpha,&
    auxig,aux1,aux2,res1,res3,indictronq,resL,resU,tU,kkapa,nb_gl
    use tailles
    use comongroup
    use residusM

    implicit none

    integer::n,np,id,jd,i,j,k,vj,ig,choix
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,inv,res,int,logGammaJ
    double precision,dimension(ngmax)::res2,res1dc,res2dc,res3dc
    double precision,dimension(np)::b,bh
    double precision,dimension(2)::k0
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3

    kkapa=k0
    choix=0
    ig=0
    k=0
    vj=0
    n=0
    j=0

    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    betaR= bh(1)**2
    etaR= bh(2)**2
    betaD= bh(3)**2
    etaD= bh(4)**2

    if(effet.eq.1) then
        theta = bh(np-nva-indic_ALPHA)*bh(np-nva-indic_ALPHA)
        if (indic_alpha.eq.1) then ! new : joint more flexible alpha = 1 
            alpha = bh(np-nva)
        else
            alpha = 1.d0
        endif
    endif

!---------------------------------------------------------
!---------- calcul de la vraisemblance ------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------cc

    res1 = 0.d0
    res2 = 0.d0
    res3 = 0.d0
    res1dc = 0.d0
    res2dc = 0.d0
    res3dc = 0.d0
    cpt = 0
    integrale1 = 0.d0
    integrale2 = 1.d0
    integrale3 = 0.d0
    aux1 = 0.d0
    aux2 = 0.d0
    resL = 0.d0
    resU = 0.d0

!*********************************************
!-----avec un effet aleatoire dans le modele
!*********************************************
    nst=2
    inv = 1.d0/theta
!     pour les donnees recurrentes
    do i=1,nsujet
        cpt(g(i))=cpt(g(i))+1  !nb obser dans groupe
        if(nva1.gt.0)then
            vet = 0.d0
            do j=1,nva1
                vet =vet + bh(np-nva+j)*dble(ve(i,j))
            end do
            vet = dexp(vet)
        else
            vet=1.d0
        endif

        if (c(i).eq.1) then ! censure par intervalle
            resL(i) = ((t1(i)/etaR)**betaR)*vet !ut1(nt1(i))*vet
            resU(i) = ((tU(i)/etaR)**betaR)*vet !ut1(ntU(i))*vet
        endif
        if ((resL(i).ne.resL(i)).or.(abs(resL(i)).ge.1.d30)) then
            !print*,"here1"
            funcpajweib_intcens=-1.d9
            goto 123
        end if
        if ((resU(i).ne.resU(i)).or.(abs(resU(i)).ge.1.d30)) then
            !print*,"here2"
            funcpajweib_intcens=-1.d9
            goto 123
        end if

        if (c(i).eq.0) then ! censure a droite
            res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet !ut1(nt1(i))*vet
        endif
        if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge.1.d30)) then
            !print*,"here3"
            funcpajweib_intcens=-1.d9
            goto 123
        end if

!     modification pour nouvelle vraisemblance / troncature
        res3(g(i)) = res3(g(i)) + ((t0(i)/etaR)**betaR)*vet !ut1(nt0(i))*vet
        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge.1.d30)) then
            !print*,"here4"
            funcpajweib_intcens=-1.d9
            goto 123
        end if
    end do

! pour le deces

    do k=1,ng
        if(nva2.gt.0)then
            vet2 = 0.d0
            do j=1,nva2
                vet2 =vet2 + bh(np-nva2+j)*dble(vedc(k,j))
            end do
            vet2 = dexp(vet2)
        else
            vet2=1.d0
        endif
        if(cdc(k).eq.1)then
            res2dc(k) = (betaD-1.d0)*dlog(t1dc(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet2) !dlog(dut2(nt1dc(k))*vet2)
        endif
        if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge.1.d30)) then
            !print*,"here5"
            funcpajweib_intcens=-1.d9
            goto 123
        end if
! pour le calcul des integrales / pour la survie, pas pour donnees recur:
        aux1(k)=((t1dc(k)/etaD)**betaD)*vet2
        aux2(k)=((t0dc(k)/etaD)**betaD)*vet2 !vraie troncature
        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge.1.d30)) then
            !print*,"here6"
            funcpajweib_intcens=-1.d9
            goto 123
        end if
        if ((aux2(k).ne.aux2(k)).or.(abs(aux2(k)).ge.1.d30)) then
            !print*,"here7"
            funcpajweib_intcens=-1.d9
            goto 123
        end if
    end do

!**************INTEGRALES ****************************
    do ig=1,ng
        auxig = ig
        choix = 1
        call gaulagJ_intcens(int,choix,nb_gl)
        integrale1(ig) = int
        if (integrale1(ig).eq.0.d0) then
            integrale1(ig) = 1.d-300
        endif
        if (indictronq.eq.1) then
            choix = 2
            call gaulagJ_intcens(int,choix,nb_gl)
            integrale2(ig) = int
        endif
    end do
!************* FIN INTEGRALES **************************

    res = 0.d0
    do k=1,ng
        if(cpt(k).gt.0)then
            if (indictronq.eq.1) then
                res = res + res2dc(k) + &
                dlog(integrale1(k))-dlog(integrale2(k))
            else
                res = res + res2dc(k) - &
                logGammaJ(1.d0/theta)-dlog(theta)/theta + &
                dlog(integrale1(k))
            endif
!*******************************************************
!     developpement de taylor d ordre 3
!*******************************************************
!        write(*,*)'******* TAYLOR *************'
            !    res= res + res2(k)+ res2dc(k) &
            !    - logGammaJ(1.d0/theta)-dlog(theta)/theta  &
            !    + dlog(integrale3(k))
            !endif
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                !print*,"here8",res,integrale1(k),integrale2(k),k,dlog(integrale1(k)),dlog(integrale2(k))
                funcpajweib_intcens=-1.d9
                goto 123
            end if
        endif
    end do


    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        !print*,"here9"
        funcpajweib_intcens = -1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123
    else
        funcpajweib_intcens = res
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=nigdc(k)
        end do
    end if

123     continue

    return

    end function funcpajweib_intcens


