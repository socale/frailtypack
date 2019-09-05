


!========================          funcpajng NEW         ====================
    double precision function funcpaGweib_intcens(b,np,id,thi,jd,thj,k0)

    !use comon,only:indictronqdc
    use comon,only:etaR,etaD,betaR,betaD, &
    t0,t1,tU,t0dc,t1dc,c,cdc,nsujet,nva,nva1,nva2,nst,&
    stra,ve,vedc,effet,ng,g,nig,AG,indic_ALPHA,theta,alpha,&
    auxig,aux1,aux2,res1,res3,indictronq,resL,resU,kkapa,nb_gl
    use tailles
    use comongroup
    use residusM

    implicit none

    integer::nb,n,np,id,jd,i,j,k,vj,cptg,l,ig,choix
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,dnb,sum,inv,res,int,logGammaJ
    double precision,dimension(ngmax)::res2,res1dc,res2dc,res3dc
    double precision,dimension(np)::b,bh
    double precision,dimension(2)::k0
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3

    kkapa = k0
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

!*******************************************
!----- A SIMPLE SHARED FRAILTY  MODEL
!      write(*,*)'SIMPLE SHARED FRAILTY MODEL'
!*******************************************
    if(indic_joint.eq.0)then
        inv = 1.d0/theta
!    i indice les sujets
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
                res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(k))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
            endif
            if((c(i).eq.1).and.(stra(i).eq.2))then
                res2(g(i)) = res2(g(i))+(betaD-1.d0)*dlog(t1(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet)
            endif
            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpaGweib_intcens=-1.d9
                goto 123
            endif

            if(stra(i).eq.1)then
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet
            endif
            if(stra(i).eq.2)then
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaD)**betaD)*vet
            endif
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpaGweib_intcens=-1.d9
                goto 123
            endif

! modification pour nouvelle vraisemblance / troncature:
            if(stra(i).eq.1)then
                res3(g(i)) = res3(g(i)) + ((t0(i)/etaR)**betaR)*vet
            endif
            if(stra(i).eq.2)then
                res3(g(i)) = res3(g(i)) + ((t0(i)/etaD)**betaD)*vet
            endif
            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                funcpaGweib_intcens=-1.d9
                goto 123
            endif
        end do

        res = 0.d0
        cptg = 0
!     gam2 = logGammaJ(inv)
! k indice les groupes

        do k=1,ng
            sum=0.d0
            if(cpt(k).gt.0)then
                nb = nig(k)
                dnb = dble(nig(k))
!     gam1 = logGammaJ(dnb + inv)
                if (dnb.gt.1.d0) then
                    do l=1,nb
                        sum=sum+dlog(1.d0+theta*dble(nb-l))
                    end do
                endif

                if(theta.gt.(1.d-5)) then
!                print*,'==== ici ===='
! ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res= res-(inv+dnb)*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        + res2(k) + sum
! nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccc
                    else
!          write(*,*)'*** funcpaGweib frailty**',res,dnb,theta,res1(k),
!     &                  res1(k),res2(k),sum,res3(k)
                        res = res-(inv+dnb)*dlog(theta*res1(k)+1.d0) &
                        +(inv)*dlog(theta*res3(k)+1.d0) + res2(k) + sum
                    endif
                else
!     developpement de taylor d ordre 3
!                   write(*,*)'************** TAYLOR *************'
! ancienne vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res = res-dnb*dlog(theta*(res1(k)-res3(k))+1.d0)&
                        -(res1(k)-res3(k))*(1.d0-theta*(res1(k)-res3(k))/2.d0 &
                        +theta*theta*(res1(k)-res3(k))*(res1(k)-res3(k))/3.d0) &
                        +res2(k)+sum
! nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res-dnb*dlog(theta*res1(k)+1.d0) &
                        -res1(k)*(1.d0-theta*res1(k)/2.d0 &
                        +theta*theta*res1(k)*res1(k)/3.d0) &
                        +res2(k)+sum+res3(k)*(1.d0-theta*res3(k)/2.d0 &
                        +theta*theta*res3(k)*res3(k)/3.d0)
                    endif
                endif
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpaGweib_intcens=-1.d9
                    goto 123
                end if
            endif
        end do
    else !passage au modele conjoint


!*******************************************
!----- JOINT FRAILTY MODEL
!*******************************************
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
                resL(i) = ((t1(i)/etaR)**betaR)*vet
                resU(i) = ((tU(i)/etaR)**betaR)*vet
            endif
            if ((resL(i).ne.resL(i)).or.(abs(resL(i)).ge.1.d30)) then
                !print*,"here1"
                funcpaGweib_intcens=-1.d9
                goto 123
            end if
            if ((resU(i).ne.resU(i)).or.(abs(resU(i)).ge.1.d30)) then
                !print*,"here2"
                funcpaGweib_intcens=-1.d9
                goto 123
            end if

            if (c(i).eq.0) then ! censure a droite
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet
            endif
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                !print*,"here3"
                funcpaGweib_intcens=-1.d9
                goto 123
            end if

!     modification pour nouvelle vraisemblance / troncature
            res3(g(i)) = res3(g(i)) + ((t0(i)/etaR)**betaR)*vet
            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                !print*,"here4"
                funcpaGweib_intcens=-1.d9
                goto 123
            end if
        end do

! pour le deces

        do k=1,lignedc!ng
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
                res2dc(gsuj(k)) = res2dc(gsuj(k))+(betaD-1.d0)*dlog(t1dc(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet2)
            endif
            if ((res2dc(gsuj(k)).ne.res2dc(gsuj(k))).or.(abs(res2dc(gsuj(k))).ge. 1.d30)) then
                !print*,"here5"
                funcpaGweib_intcens=-1.d9
                goto 123
            end if

! pour le calcul des integrales / pour la survie, pas pour donnees recur:
            aux1(gsuj(k))=aux1(gsuj(k))+((t1dc(k)/etaD)**betaD)*vet2
            aux2(gsuj(k))=aux2(gsuj(k))+((t0dc(k)/etaD)**betaD)*vet2 !vraie troncature
            if ((aux1(gsuj(k)).ne.aux1(gsuj(k))).or.(abs(aux1(gsuj(k))).ge. 1.d30)) then
                !print*,"here6"
                funcpaGweib_intcens=-1.d9
                goto 123
            end if
            if ((aux2(gsuj(k)).ne.aux2(gsuj(k))).or.(abs(aux2(gsuj(k))).ge. 1.d30)) then
                !print*,"here7"
                funcpaGweib_intcens=-1.d9
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
                    dlog(integrale1(k)) - dlog(integrale2(k))
                else
                    res = res + res2dc(k) - &
                    logGammaJ(1.d0/theta)-dlog(theta)/theta + &
                    dlog(integrale1(k))
                endif

                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpaGweib_intcens=-1.d9
                    goto 123
                end if
            endif
        end do
    endif !fin boucle indic_joint=0 or 1

!-------------------------------------------------------

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        !print*,"here9"
        funcpaGweib_intcens = -1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123
    else
        funcpaGweib_intcens = res
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=nigdc(k)
        end do
    end if

123     continue

    return

    end function funcpaGweib_intcens


