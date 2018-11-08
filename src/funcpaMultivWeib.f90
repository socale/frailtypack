
!!!!_____________________________________________________
!========================          FUNCPA NEW         ====================
    double precision function funcpaMultivWeib(b,np,id,thi,jd,thj,k0)

    use taillesmultiv
    use comonmultiv
    use residusMmultiv

    implicit none

! *** NOUVELLLE DECLARATION F90 :

    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(3)::k0
    double precision,intent(in)::thi,thj
    integer::n,i,j,k,vj,ig,choix
    integer,dimension(ngmax)::cpt,cptmeta
    double precision::res,vet,vet2,vet3
    double precision,dimension(np)::bh
    double precision,dimension(ngmax)::res2,res1dc,res2dc,res2meta &
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

    n = (np-nva-effet-indic_ALPHA)/nst

    betaR= bh(1)**2
    etaR= bh(2)**2
    betaD= bh(3)**2
    etaD= bh(4)**2
    betaM= bh(5)**2
    etaM= bh(6)**2

    if(effet.eq.1) then
        theta = bh(np-nva-indic_ALPHA)*bh(np-nva-indic_ALPHA)
!        alpha = bh(np-nva)
!add
        alpha =bh(np-nva-indic_rho) !rho
        alpha1=bh(np-nva-indic_a1)
        alpha2=bh(np-nva-indic_a2)
        eta = bh(np-nva-indic_eta)*bh(np-nva-indic_eta)
    endif

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
!add
        res1meta(k) = 0.d0
        res2meta(k) = 0.d0
        res3meta(k) = 0.d0
        cptmeta(k) = 0
    end do

!*********************************************
!-----avec un effet aleatoire dans le modele
!*********************************************

!    inv = 1.d0/theta

!==========================================================================
!     pour les donnees loco
!==========================================================================

    do i=1,nsujet
        cpt(g(i))=cpt(g(i))+1
        if(nva1.gt.0)then
            vet = 0.d0
            do j=1,nva1
                vet =vet + bh(np-nva+j)*dble(ve(i,j))
            end do
            vet = dexp(vet)
        else
            vet=1.d0
        endif

        if((c(i).eq.1))then
            res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(i))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpaMultivWeib=-1.d9
                goto 123
            end if
        endif

!     nouvelle version
        res1(g(i)) = res1(g(i))+((t1(i)/etaR)**betaR)*vet
         if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
            funcpaMultivWeib=-1.d9
            goto 123
        end if

!     modification pour nouvelle vraisemblance / troncature:
        res3(g(i)) = res3(g(i))+((t0(i)/etaR)**betaR)*vet
        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
            funcpaMultivWeib=-1.d9
            goto 123
        end if
    end do

!==========================================================================
!     pour les donnees dc
!==========================================================================

    do k=1,ng
        if(nva2.gt.0)then
            vet2 = 0.d0
            do j=1,nva2
                vet2 =vet2 + bh(np-nva2-nva3+j)*dble(vedc(k,j))
            end do
            vet2 = dexp(vet2)
        else
            vet2=1.d0
        endif
        if(cdc(k).eq.1)then
            res2dc(k) = (betaD-1.d0)*dlog(t1dc(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet2)
            if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
                funcpaMultivWeib=-1.d9
                goto 123
            end if
        endif

! pour le calcul des integrales / pour la survie, pas les donnees recurrentes:
        aux1(k)=((t1dc(k)/etaD)**betaD)*vet2
        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
            funcpaMultivWeib=-1.d9
            goto 123
        end if
    end do

!==========================================================================
!     pour les donnees loco
!==========================================================================

    do i=1,nsujetmeta
        cptmeta(gmeta(i))=cptmeta(gmeta(i))+1
        if(nva3.gt.0)then
            vet3 = 0.d0
            do j=1,nva3
                vet3 =vet3 + bh(np-nva3+j)*dble(vemeta(i,j))
            end do
            vet3 = dexp(vet3)
        else
            vet3=1.d0
        endif

        if((cmeta(i).eq.1))then
            res2meta(gmeta(i)) = res2meta(gmeta(i))+(betaM-1.d0)*dlog(t1meta(i))+dlog(betaM)-betaM*dlog(etaM)+dlog(vet3)
            if ((res2meta(gmeta(i)).ne.res2meta(gmeta(i))).or.(abs(res2meta(gmeta(i))).ge. 1.d30)) then
                funcpaMultivWeib=-1.d9
                goto 123
            end if
        endif

!     nouvelle version
        res1meta(gmeta(i)) = res1meta(gmeta(i))+((t1meta(i)/etaM)**betaM)*vet3
         if ((res1meta(gmeta(i)).ne.res1meta(gmeta(i))).or.(abs(res1meta(gmeta(i))).ge. 1.d30)) then
            funcpaMultivWeib=-1.d9
            goto 123
        end if

!     modification pour nouvelle vraisemblance / troncature:
        res3meta(gmeta(i)) = res3meta(gmeta(i))+((t0meta(i)/etaM)**betaM)*vet3
        if ((res3meta(gmeta(i)).ne.res3meta(gmeta(i))).or.(abs(res3meta(gmeta(i))).ge. 1.d30)) then
            funcpaMultivWeib=-1.d9
            goto 123
        end if

    end do

!==========================================================================
!     fin pour les donnees Meta
!==========================================================================

!**************INTEGRALES ****************************
    do ig=1,ng
        auxig=ig
!        choix = 3
        call gausshermiteBIS2011(int,30)!!!!! garder 30!!!!
!        call gaulagJ(int,choix)
        integrale3(ig) = int !moins bon
    end do
!************* FIN INTEGRALES **************************

    res = 0.d0

    do k=1,ng
        res = res + res2(k)+res2meta(k)+res2dc(k)+dlog(integrale3(k))-dlog((2.d0)*pi*sqrt(theta*eta* &
             (1.d0-((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0)**2)))
    end do


    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpaMultivWeib =-1.d9
        do k=1,ng
            Rrec(k)=0.d0
            Nrec(k)=0
            Rdc(k)=0.d0
            Ndc(k)=0
            Rrec2(k)=0.d0
            Nrec2(k)=0
        end do
        goto 123
    else
        funcpaMultivWeib = res
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=cdc(k)
            Rrec2(k)=res1meta(k)
            Nrec2(k)=nigmeta(k)
        end do
    end if

!Ad:
123     continue

    return

    end function funcpaMultivWeib

!=================================================================================================
