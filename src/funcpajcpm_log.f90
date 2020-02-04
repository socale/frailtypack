


!========================          FUNCPAJ_SPLINES         ====================
    double precision function funcpajcpm_log(b,np,id,thi,jd,thj,k0)

    use tailles
    ! use comon,only:nst,AG,t0dc,res4
    use comon,only:nbintervR,nbintervDC,ttt,tttdc,betacoef, &
    t0,t1,t1dc,c,cdc,nsujet,nva,nva1,nva2, &
    effet,stra,ve,vedc,ng,g,nig,indic_ALPHA,ALPHA,sig2, &
    auxig,aux1,aux2,res1,res3,kkapa,nstRec,nb_gh
    use residusM
    !use comongroup,only:the1,the2
    use comongroup,only:vet,vet2

    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :

    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(2),intent(in)::k0
    double precision,intent(in)::thi,thj

    integer::n,i,j,k,vj,ig,choix,jj,gg,gg2
    double precision::som11,som21
    integer,dimension(ngmax)::cpt
    double precision::res

    double precision,dimension(np)::bh
    double precision,dimension(ngmax)::res2,res1dc,res2dc &
    ,res3dc,integrale1,integrale2,integrale3

    double precision::int
    double precision,parameter::pi=3.141592653589793d0

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

    betacoef = 0.d0
    do i = 1,(nbintervR*nstRec+nbintervDC)
        betacoef(i)=bh(i)**2
    end do

    if(effet.eq.1) then
        sig2 = bh(np-nva-indic_ALPHA)*bh(np-nva-indic_ALPHA)
        if (indic_alpha.eq.1) then ! new : joint more flexible alpha = 1
            alpha = bh(np-nva)
        else
            alpha = 1.d0
        endif
    endif

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
        integrale3(k) = 0.d0
        aux1(k)=0.d0
        aux2(k)=0.d0
    end do

!*********************************************
!-----avec un effet aleatoire dans le modele
!*********************************************

!ccccccccccccccccccccccccccccccccccccccccc
!     pour les donnees recurrentes
!ccccccccccccccccccccccccccccccccccccccccc

    do i=1,nsujet

        cpt(g(i))=cpt(g(i))+1

        if(nva1.gt.0)then
            vet = 0.d0
            do j=1,nva1
                vet = vet + bh(np-nva+j)*dble(ve(i,j))
            end do
            vet = dexp(vet)
        else
            vet = 1.d0
        endif

        if((c(i).eq.1))then
            do gg=1,nbintervR
                if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                    res2(g(i)) =  res2(g(i))+ dlog(betacoef(nbintervR*(stra(i)-1)+gg)*vet)
                endif
            end do
        endif

        if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
            funcpajcpm_log=-1.d9
            goto 123
        end if

!     nouvelle version
        som11=0.d0
        som21=0.d0
        gg2=0
        do gg=1,nbintervR
            if((t1(i).gt.(ttt(gg-1))).and.(t1(i).le.(ttt(gg))))then
                som11=betacoef(nbintervR*(stra(i)-1)+gg)*(t1(i)-ttt(gg-1))
                gg2=gg
                if (gg2.ge.2)then
                    do jj=1,gg2-1
                        som21=som21+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                    end do
                endif
                res1(g(i)) = res1(g(i))+vet*(som11+som21)
            endif
        end do

        if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
            funcpajcpm_log=-1.d9
            goto 123
        end if

!     modification pour nouvelle vraisemblance / troncature:
        som11=0.d0
        som21=0.d0
        gg2=0
        do gg=1,nbintervR
            if ((t0(i).ge.(ttt(gg-1))).and.(t0(i).lt.(ttt(gg)))) then
                som11=betacoef(nbintervR*(stra(i)-1)+gg)*(t0(i)-ttt(gg-1))
                gg2=gg
                if (gg2.ge.2)then
                    do jj=1,gg2-1
                        som21=som21+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                    end do
                endif
                res3(g(i)) = res3(g(i))+vet*(som11+som21)
            endif
        end do

        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
            funcpajcpm_log=-1.d9
            goto 123
        end if
    end do

!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces
!ccccccccccccccccccccccccccccccccccccccccc

    do k=1,ng ! dans Joint ng=nb individus
        if(nva2.gt.0)then
            vet2 = 0.d0
            do j=1,nva2
                vet2 = vet2 + bh(np-nva2+j)*dble(vedc(k,j))
            end do
            vet2 = dexp(vet2)
        else
            vet2 = 1.d0
        endif

        if(cdc(k).eq.1)then
            do gg=1,nbintervDC
                if ((t1dc(k).gt.(tttdc(gg-1))).and.(t1dc(k).le.(tttdc(gg)))) then
                    res2dc(k) = dlog(betacoef(nbintervR*nstRec+gg)*vet2)
                endif
            end do
            if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
                funcpajcpm_log=-1.d9
                goto 123
            end if
        endif

! pour le calcul des integrales / pour la survie, pas les donnees recurrentes:
        som11=0.d0
        som21=0.d0
        gg2=0
        do gg=1,nbintervDC
            if((t1dc(k).gt.(tttdc(gg-1))).and.(t1dc(k).le.(tttdc(gg))))then
                som11=betacoef(gg+nbintervR*nstRec)*(t1dc(k)-tttdc(gg-1))
                gg2=gg
                if (gg2.ge.2)then
                    do jj=1,gg2-1
                        som21=som21+betacoef(jj+nbintervR*nstRec)*(tttdc(jj)-tttdc(jj-1))
                    end do
                endif
            endif
            aux1(k)=(som11+som21)*vet2
        end do

        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
            funcpajcpm_log=-1.d9
            goto 123
        end if
    end do

!**************INTEGRALES ****************************
    do ig=1,ng
        auxig = ig
        choix = 3
        call gauherJ(int,choix,nb_gh)
        integrale3(ig) = int
    end do
!************* FIN INTEGRALES **************************

    res = 0.d0
    do k=1,ng
        if(cpt(k).gt.0)then
            !if(theta.gt.(1.d-8)) then
!cccc ancienne vraisemblance : pour calendar sans vrai troncature cccccccc
                res= res + res2(k) &
!--      pour le deces:
                + res2dc(k)  &
                - dlog(dsqrt(sig2))-dlog(2.d0*pi)/2.d0  &
                + dlog(integrale3(k))
            !else
!*************************************************************************
!     developpement de taylor d ordre 3
!*************************************************************************
!                   write(*,*)'************** TAYLOR *************'
            !    res= res + res2(k) &
            !    + res2dc(k)  &
            !    - logGammaJ(1./theta)-dlog(theta)/theta  &
            !    + dlog(integrale3(k))
            !endif
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpajcpm_log=-1.d9
                goto 123
            end if
        endif
    end do

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpajcpm_log=-1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123
    else
        funcpajcpm_log = res
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=cdc(k)
        end do
    end if

123     continue

    return

    end function funcpajcpm_log

