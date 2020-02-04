


!========================          FUNCPAJ_WEIB         ====================
    double precision function funcpajweib_log(b,np,id,thi,jd,thj,k0)

    use tailles
    !use comon,only:AG,betaR,etaR,nst,t0dc
    use comon,only:etaD,betaD,etaT,betaT,nstRec, &
    t0,t1,t1dc,c,cdc,nsujet,nva,nva1,nva2, &
    effet,stra,ve,vedc,ng,g,nig,indic_ALPHA,ALPHA,sig2, &
    auxig,aux1,aux2,res1,res3,kkapa,nb_gh
    use residusM
    use comongroup,only:vet,vet2

    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :

    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(2),intent(in)::k0
    double precision,intent(in)::thi,thj
    integer::n,i,j,k,vj,ig,choix,jj
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

    do i=1,np
        bh(i)=b(i)
    end do

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    do jj=1,nstRec !en plus strates A.Lafourcade 07/2014 
        betaT(jj)=bh((jj-1)*2+1)**2
        etaT(jj)= bh((jj-1)*2+2)**2
    end do
    betaD= bh(2*nstRec+1)**2
    etaD= bh(2*nstRec+2)**2

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

        if (c(i).eq.1) then
            res2(g(i)) = res2(g(i))+(betaT(stra(i))-1.d0)*dlog(t1(i))+ &
            dlog(betaT(stra(i)))-betaT(stra(i))*dlog(etaT(stra(i)))+dlog(vet)
        endif
        if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
            !print*,"here6"
            funcpajweib_log=-1.d9
            goto 123
        end if

!     nouvelle version
        res1(g(i)) = res1(g(i)) + ((t1(i)/etaT(stra(i)))**betaT(stra(i)))*vet
        if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
            !print*,"here5"
            funcpajweib_log=-1.d9
            goto 123
        end if

!     modification pour nouvelle vraisemblance / troncature:
        res3(g(i)) = res3(g(i)) + ((t0(i)/etaT(stra(i)))**betaT(stra(i)))*vet
        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
            !print*,"here4"
            funcpajweib_log=-1.d9
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
            res2dc(k) = (betaD-1.d0)*dlog(t1dc(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet2) !dlog(dut2(nt1dc(k))*vet2)
            if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
                !print*,"here3"
                funcpajweib_log=-1.d9
                goto 123
            end if
        endif

! pour le calcul des integrales / pour la survie, pas les donnees recurrentes:
        aux1(k)=((t1dc(k)/etaD)**betaD)*vet2 !ut2(nt1dc(k))*vet2
        aux2(k)=aux2(k)+((t0(k)/etaD)**betaD)*vet2 !ut2(nt0(k))*vet2 !vraie troncature

        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
            !print*,"her2",aux1(k),k,vet2
            funcpajweib_log=-1.d9
            goto 123
        end if
        if ((aux2(k).ne.aux2(k)).or.(abs(aux2(k)).ge. 1.d30)) then
            !print*,"here1"
            funcpajweib_log=-1.d9
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
                !print*,"here"
                funcpajweib_log=-1.d9
                goto 123
            end if
        endif
    end do

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpajweib_log=-1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123
    else
        funcpajweib_log = res
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

    end function funcpajweib_log

