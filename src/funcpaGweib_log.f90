


!========================          FUNCPAG_WEIB         ====================
    double precision function funcpaGweib_log(b,np,id,thi,jd,thj,k0)

    use tailles
    use comon,only:etaR,etaD,betaR,betaD, &
    t0,t1,t0dc,t1dc,c,cdc,nsujet,nva,nva1,nva2,nst, &
    effet,stra,ve,vedc,ng,g,nig,AG,indic_ALPHA,ALPHA,sig2, &
    auxig,aux1,aux2,res1,res3,kkapa,res5,indictronq,nb_gh
    use residusM
    use comongroup

    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :

    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(2),intent(in)::k0
    double precision,intent(in)::thi,thj
    integer::n,i,j,k,vj,ig,choix
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

    if(indic_joint.eq.0)then
        if(nst==2) then
            betaR= bh(1)**2
            etaR= bh(2)**2
            betaD= bh(3)**2
            etaD= bh(4)**2
        else
            betaR= bh(1)**2
            etaR= bh(2)**2
        end if
    else
        betaR= bh(1)**2
        etaR= bh(2)**2
        betaD= bh(3)**2
        etaD= bh(4)**2
    end if

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

    res1 = 0.d0
    res2 = 0.d0
    res3 = 0.d0
    res1dc = 0.d0
    res2dc = 0.d0
    res3dc = 0.d0
    cpt = 0
    integrale1 = 0.d0
    integrale2 = 0.d0
    integrale3 = 0.d0
    aux1 = 0.d0
    aux2 = 0.d0
    res5 = 0.d0

!*******************************************
!----- A SIMPLE SHARED FRAILTY  MODEL
!      write(*,*)'SIMPLE SHARED FRAILTY MODEL'
!*******************************************

    if(indic_joint.eq.0)then

        do i=1,nsujet

            cpt(g(i))=cpt(g(i))+1

            if(nva.gt.0)then
                vet = 0.d0
                do j=1,nva
                    vet =vet + bh(np-nva+j)*dble(ve(i,j))
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif

            if((c(i).eq.1).and.(stra(i).eq.1))then
                res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(i))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
            endif

            if((c(i).eq.1).and.(stra(i).eq.2))then
                res2(g(i)) = res2(g(i))+(betaD-1.d0)*dlog(t1(i))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge.1.d30)) then
                funcpaGweib_log=-1.d9
                goto 123
            end if

! modification pour nouvelle vraisemblance / troncature:
            if(stra(i).eq.1)then
                res3(g(i)) = res3(g(i)) + ((t0(i)/etaR)**betaR)*vet
                res5(i) = ((t1(i)/etaR)**betaR)*vet
                res1(g(i)) = res1(g(i)) + res5(i) ! pour les residus
            endif

            if(stra(i).eq.2)then
                res3(g(i)) = res3(g(i)) + ((t0(i)/etaD)**betaD)*vet
                res5(i) = ((t1(i)/etaD)**betaD)*vet
                res1(g(i)) = res1(g(i)) + res5(i) ! pour les residus
            endif

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge.1.d30)) then
                funcpaGweib_log=-1.d9
                goto 123
            end if

            if ((res5(i).ne.res5(i)).or.(abs(res5(i)).ge.1.d30)) then
                funcpaGweib_log=-1.d9
                goto 123
            end if
        end do

!**************INTEGRALES ****************************
        do ig=1,ng
            auxig = ig
            choix = 1
            call gauherS(int,choix,nb_gh)
            integrale1(ig) = int
            if (AG.eq.1) then
                choix = 3
                call gauherS(int,choix,nb_gh)
                integrale3(ig) = int
            endif
            if (indictronq.eq.1) then
                choix = 2
                call gauherS(int,choix,nb_gh)
                integrale2(ig) = int
            endif
        end do
!************* FIN INTEGRALES ************************

        res = 0.d0

        do k=1,ng
            if(AG.EQ.1)then
                res = res+res2(k)-dlog(dsqrt(sig2))-dlog(2.d0*pi)/2.d0+ &
                dlog(integrale3(k))
            else
                res = res+res2(k)+dlog(integrale1(k)) &
                -dlog(integrale2(k))
            endif
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpaGweib_log=-1.d9
                goto 123
            end if
        enddo
    else ! passage au modele conjoint

!*********************************************
!----- JOINT FRAILTY MODEL
!*********************************************

!     pour les donnees recurrentes

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
                res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(i))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                !print*,"here7"
                funcpaGweib_log=-1.d9
                goto 123
            end if

!     nouvelle version
            res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                !print*,"here6"
                funcpaGweib_log=-1.d9
                goto 123
            end if

!     modification pour nouvelle vraisemblance / troncature:
            res3(g(i)) = res3(g(i)) + ((t0(i)/etaR)**betaR)*vet
            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                !print*,"here5"
                funcpaGweib_log=-1.d9
                goto 123
            end if
        end do

!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces
!ccccccccccccccccccccccccccccccccccccccccc

        do k=1,lignedc !ng
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
                res2dc(gsuj(k)) = res2dc(gsuj(k)) + (betaD-1.d0)*dlog(t1dc(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet2)
                if ((res2dc(gsuj(k)).ne.res2dc(gsuj(k))).or.(abs(res2dc(gsuj(k))).ge. 1.d30)) then
                    !print*,"here3",res2dc(gsuj(k))
                    funcpaGweib_log=-1.d9
                    goto 123
                end if
            endif

! pour le calcul des integrales / pour la survie, pas les donnees recurrentes:
            aux1(gsuj(k)) = aux1(gsuj(k)) + ((t1dc(k)/etaD)**betaD)*vet2
            aux2(gsuj(k)) = aux2(gsuj(k)) + ((t0dc(k)/etaD)**betaD)*vet2 !vraie troncature

            if ((aux1(gsuj(k)).ne.aux1(gsuj(k))).or.(abs(aux1(gsuj(k))).ge. 1.d30)) then
                !print*,"here1",aux1(gsuj(k)),vet2,k
                funcpaGweib_log=-1.d9
                goto 123
            end if
            if ((aux2(gsuj(k)).ne.aux2(gsuj(k))).or.(abs(aux2(gsuj(k))).ge. 1.d30)) then
                !print*,"here2"
                funcpaGweib_log=-1.d9
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
!************* FIN INTEGRALES ************************

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
                    funcpaGweib_log=-1.d9
                    goto 123
                end if
            endif
        end do
    endif ! fin boucle indic_joint=0 ou 1


    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpaGweib_log=-1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123
    else
        funcpaGweib_log = res
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=cdc(k)
        end do
    end if

123     continue

    return

    end function funcpaGweib_log

