

!========================          funcpajng NEW         ====================
    double precision function funcpaGcpm(b,np,id,thi,jd,thj,k0)

    !use comon,only:cens
    use comon,only:nbintervR,nbintervDC,t0,t1,t0dc,t1dc,c,cdc,nsujet,nva,nva1,nva2, &
    nst,effet,stra,ve,vedc,ng,g,nig,AG,indic_ALPHA,alpha,theta,auxig,aux1,aux2,res1,res3, &
    ttt,tttdc,betacoef,kkapa,nb_gl
    use tailles
    use comongroup
    use residusM

    implicit none

    double precision::som11,som21
    integer::jj,gg,gg2
    integer::nb,np,id,jd,i,j,k,cptg,l,ig,choix
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,dnb,sum,inv,res,int,logGammaJ
    double precision,dimension(ngmax)::res2,res1dc,res2dc,res3dc
    double precision,dimension(np)::b,bh
    double precision,dimension(2)::k0
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3


    bh=b
    kkapa=k0

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

! Allocation des parametre des fcts de risque cte par morceaux
!cpm
    betacoef = 0.d0
    if(indic_joint.eq.0)then
        if(nst==2) then
            do i = 1,(2*nbintervR)
                betacoef(i)=bh(i)**2
            end do
        else
            do i = 1,(nbintervR)
                betacoef(i)=bh(i)**2
            end do
        end if
    else
        do i = 1,(nbintervR+nbintervDC)
            betacoef(i)=bh(i)**2
        end do
    end if
!cpm
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
!----- A SIMPLE SHARED FRAILTY  MODEL  
!      write(*,*)'SIMPLE SHARED FRAILTY MODEL'
!*********************************************
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
                do gg=1,nbintervR
                    if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                        res2(g(i)) =  res2(g(i))+ dlog(betacoef(gg)*vet)
                    endif
                end do
            endif

            if((c(i).eq.1).and.(stra(i).eq.2))then
                do gg=1,nbintervR
                    if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                        res2(g(i)) =  res2(g(i))+ dlog(betacoef(gg+nbintervR)*vet)
                    endif
                end do
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpaGcpm=-1.d9
                goto 123
            end if

            if(stra(i).eq.1)then
                som11=0.d0
                som21=0.d0
                gg2=0
                do gg=1,nbintervR
                    if((t1(i).gt.(ttt(gg-1))).and.(t1(i).le.(ttt(gg))))then
                        som11=betacoef(gg)*(t1(i)-ttt(gg-1))
                        gg2=gg
                        if (gg2.ge.2)then
                            do jj=1,gg2-1
                                som21=som21+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                            end do
                        endif
                        res1(g(i)) = res1(g(i))+vet*(som11+som21)
                    endif
                end do
            endif

            if(stra(i).eq.2)then
                som11=0.d0
                som21=0.d0
                gg2=0
                do gg=1,nbintervR
                    if((t1(i).gt.(ttt(gg-1))).and.(t1(i).le.(ttt(gg))))then
                        som11=betacoef(gg+nbintervR)*(t1(i)-ttt(gg-1))
                        gg2=gg
                        if (gg2.ge.2)then
                            do jj=1,gg2-1
                                som21=som21+betacoef(jj+nbintervR)*(ttt(jj)-ttt(jj-1))
                            end do
                        endif
                        res1(g(i)) = res1(g(i))+vet*(som11+som21)
                    endif
                end do
            endif
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpaGcpm=-1.d9
                goto 123
            end if

! modification pour nouvelle vraisemblance / troncature:
            if(stra(i).eq.1)then
                som11=0.d0
                som21=0.d0
                gg2=0
                do gg=1,nbintervR
                    if((t0(i).gt.(ttt(gg-1))).and.(t0(i).le.(ttt(gg))))then
                        som11=betacoef(gg)*(t0(i)-ttt(gg-1))
                        gg2=gg
                        if (gg2.ge.2)then
                            do jj=1,gg2-1
                                som21=som21+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                            end do
                        endif
                        res3(g(i)) = res3(g(i))+vet*(som11+som21)
                    endif
                end do
            endif

            if(stra(i).eq.2)then
                som11=0.d0
                som21=0.d0
                gg2=0
                do gg=1,nbintervR
                    if((t0(i).gt.(ttt(gg-1))).and.(t0(i).le.(ttt(gg))))then
                        som11=betacoef(gg+nbintervR)*(t0(i)-ttt(gg-1))
                        gg2=gg
                        if (gg2.ge.2)then
                            do jj=1,gg2-1
                                som21=som21+betacoef(jj+nbintervR)*(ttt(jj)-ttt(jj-1))
                            end do
                        endif
                        res3(g(i)) = res3(g(i))+vet*(som11+som21)
                    endif
                end do
            endif
            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                funcpaGcpm=-1.d9
                goto 123
            end if
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
! ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res= res-(inv+dnb)*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        + res2(k) + sum  
! nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccc
                    else
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
                    funcpaGcpm=-1.d9
                    goto 123
                end if
            endif 
        end do
    else !passage au modele conjoint

!*********************************************
!----- JOINT FRAILTY MODEL
!*********************************************
        nst=2
        inv = 1.d0/theta
!     pour les donnees recurrentes
!
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

            if((c(i).eq.1))then
                do gg=1,nbintervR
                    if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                        res2(g(i)) =  res2(g(i))+ dlog(betacoef(gg)*vet)
                    endif
                end do
            endif
            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpaGcpm=-1.d9
                goto 123
            end if

!     nouvelle version
            som11=0.d0
            som21=0.d0
            gg2=0
            do gg=1,nbintervR
                if((t1(i).gt.(ttt(gg-1))).and.(t1(i).le.(ttt(gg))))then
                    som11=betacoef(gg)*(t1(i)-ttt(gg-1))
                    gg2=gg
                    if (gg2.ge.2)then
                        do jj=1,gg2-1
                            som21=som21+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif
                    res1(g(i)) = res1(g(i))+vet*(som11+som21)
                endif
            end do

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpaGcpm=-1.d9
                goto 123
            end if

!     modification pour nouvelle vraisemblance / troncature:
            som11=0.d0
            som21=0.d0
            gg2=0
            do gg=1,nbintervR
                if((t0(i).gt.(ttt(gg-1))).and.(t0(i).le.(ttt(gg))))then
                    som11=betacoef(gg)*(t0(i)-ttt(gg-1))
                    gg2=gg
                    if (gg2.ge.2)then
                        do jj=1,gg2-1
                            som21=som21+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif
                    res3(g(i)) = res3(g(i))+vet*(som11+som21)
                endif
            end do
            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                funcpaGcpm=-1.d9
                goto 123
            end if
        end do

! pour le deces
!
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
                do gg=1,nbintervDC
                    if ((t1dc(k).gt.(tttdc(gg-1))).and.(t1dc(k).le.(tttdc(gg)))) then
                        res2dc(gsuj(k)) = res2dc(gsuj(k))+dlog(betacoef(nbintervR+gg)*vet2)
                    endif
                end do

                if ((res2dc(gsuj(k)).ne.res2dc(gsuj(k))).or.(abs(res2dc(gsuj(k))).ge. 1.d30)) then
                    funcpaGcpm=-1.d9
                    goto 123
                end if
            endif

! pour le calcul des integrales / pour la survie, pas pour donnees recur:
            som11=0.d0
            som21=0.d0
            gg2=0
            do gg=1,nbintervDC
                if((t1dc(k).gt.(tttdc(gg-1))).and.(t1dc(k).le.(tttdc(gg))))then
                    som11=betacoef(gg+nbintervR)*(t1dc(k)-tttdc(gg-1))
                    gg2=gg
                    if (gg2.ge.2)then
                        do jj=1,gg2-1
                            som21=som21+betacoef(jj+nbintervR)*(tttdc(jj)-tttdc(jj-1))
                        end do
                    endif
                endif
                aux1(gsuj(k))=aux1(gsuj(k))+(som11+som21)*vet2
            end do

            som11=0.d0
            som21=0.d0
            gg2=0
            do gg=1,nbintervDC
                if((t0dc(k).gt.(tttdc(gg-1))).and.(t0dc(k).le.(tttdc(gg))))then
                    som11=betacoef(gg+nbintervR)*(t0dc(k)-tttdc(gg-1))
                    gg2=gg
                    if (gg2.ge.2)then
                        do jj=1,gg2-1
                            som21=som21+betacoef(jj+nbintervR)*(tttdc(jj)-tttdc(jj-1))
                        end do
                    endif
                endif
                aux2(gsuj(k))=aux2(gsuj(k))+(som11+som21)*vet2
            end do

            if ((aux1(gsuj(k)).ne.aux1(gsuj(k))).or.(abs(aux1(gsuj(k))).ge. 1.d30)) then
                funcpaGcpm=-1.d9
                goto 123
            end if
            if ((aux2(gsuj(k)).ne.aux2(gsuj(k))).or.(abs(aux2(gsuj(k))).ge. 1.d30)) then
                funcpaGcpm=-1.d9
                goto 123
            end if

        end do

!**************INTEGRALES ****************************
        do ig=1,ng
            auxig=ig
            choix = 3
            call gaulagJ(int,choix,nb_gl)
            integrale3(ig) = int !moins bon 
! parfois quand bcp de deces par groupe integrale3=0
!            if(integrale3(ig).lt.1.d-300)then
!                integrale3(ig) = 1.d-300
!            endif
        end do
!************* FIN INTEGRALES **************************
        res = 0.d0
        do k=1,ng
            sum=0.d0
            if(cpt(k).gt.0)then
                res= res + res2(k) + res2dc(k) &
                - logGammaJ(1.d0/theta)-dlog(theta)/theta &
                + dlog(integrale3(k))

                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpaGcpm=-1.d9
                    goto 123
                end if
            endif
        end do
    endif !fin boucle indic_joint=0 or 1

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpaGcpm = -1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123
    else
        funcpaGcpm = res
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=nigdc(k)
        end do
    end if

123     continue

    return

    end function funcpaGcpm

