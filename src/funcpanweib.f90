
!========================          FUNCPA  NESTED_SPLINES       ====================


    double precision function funcpanweib(b,np,id,thi,jd,thj,k0)

    use tailles
    use comon,only:t0,t1,c,nsujet,nva,nst,stra,effet,ve &
    ,g,nig,AG,auxig,alpha,eta,kkapa,etaR,etaD,betaR,betaD,indictronq,nb_gl
    use commun,only:ngexact,mij,mid,ssg,aux1,aux2
    use residusM

    Implicit none

    integer::nb,np,id,jd,i,j,k,cptg,l,ig,ip,issg,choix
    integer,dimension(ngexact)::cpt
    double precision::thi,thj,dnb,sum,theta,inv,res,vet,int,gammlnN
    double precision,dimension(np)::b,bh
    double precision,dimension(ngexact)::res1,res2,res3 &
    ,integrale1,integrale2,integrale3,sum1
    double precision,dimension(2)::k0

    kkapa=k0
    bh=b
    j=0
    theta=0.d0
    res=0.d0
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    if (nst == 1) then
        betaR= bh(1)**2
        etaR= bh(2)**2
        betaD= 0.d0
        etaD=  0.d0
    else
        betaR= bh(1)**2
        etaR= bh(2)**2
        betaD= bh(3)**2
        etaD= bh(4)**2
    end if

    if(effet.eq.1) then
        theta = (bh(np-nva)*bh(np-nva)) ! variance effet groupe
    endif

    if(effet.eq.2) then
        alpha = (bh(np-nva-1)*bh(np-nva-1)) ! variance effet groupe
        eta = (bh(np-nva)*bh(np-nva))  ! variance effet sous groupe
    endif


!--------------------------------------------------------
!----------calcul de la vraisemblance ------------------
!---------------------------------------------------------

    res1 = 0.d0
    res2 = 0.d0
    cpt = 0

!*******************************************
!----- sans effet aleatoire dans le modele
!*******************************************

    if (effet.eq.0) then
        do i=1,nsujetmax
            cpt(g(i))=cpt(g(i))+1

            if(nva.gt.0)then
                vet = 0.d0
                do j=1,nva
                    vet =vet + bh(np-nva+j)*ve(i,j)
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

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpanweib=-1.d9
                goto 123
            end if

            if(stra(i).eq.1)then
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet - ((t0(i)/etaR)**betaR)*vet
            endif

            if(stra(i).eq.2)then
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaD)**betaD)*vet - ((t0(i)/etaD)**betaD)*vet
            endif

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpanweib=-1.d9
                goto 123
            end if
        end do

        res = 0.d0
        cptg = 0
! k indice les groupes
        do k=1,ngexact
            if(cpt(k).gt.0)then !nb de sujets dans un gpe=nig()
                res = res-res1(k)+ res2(k) 
                cptg = cptg + 1
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpanweib=-1.d9
                    goto 123
                end if
            endif
        end do
    endif !fin boucle effet=0

!**************************************************
!-----avec un seul  effet aleatoire dans le modele
!**************************************************

    if (effet.eq.1) then
!      write(*,*)'AVEC 1 EFFET ALEATOIRE'

        inv = 1.d0/theta
        cpt=0
        res1=0.d0
        res2=0.d0
        res3=0.d0

        do i=1,nsujetmax

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

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpanweib=-1.d9
                goto 123
            end if

            if(stra(i).eq.1)then
! nouvelle version
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet
            endif

            if(stra(i).eq.2)then
! nouvelle version
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaD)**betaD)*vet
            endif

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpanweib=-1.d9
                goto 123
            end if

! modification pour nouvelle vraisemblance / troncature:
            if(stra(i).eq.1)then
                res3(g(i)) = res3(g(i)) + ((t0(i)/etaR)**betaR)*vet
            endif

            if(stra(i).eq.2)then
                res3(g(i)) = res3(g(i)) + ((t0(i)/etaD)**betaD)*vet
            endif

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                funcpanweib=-1.d9
                goto 123
            end if
        end do

        res = 0.d0
        cptg = 0
        mid = 0
!     gam2 = gamma(inv)
! k indice les groupes

        do k=1,nsujet
            if(c(k).eq.1)then
                mid(g(k))=mid(g(k))+1
            endif
        end do

        do k=1,ngexact
            sum=0.d0
            if(cpt(k).gt.0)then

                nb = mid(k)!nb de deces par groupe
                dnb = dble(nb)

                if (dnb.gt.1.d0) then
                    do l=1,nb
                        sum=sum+dlog(1.d0+theta*dble(nb-l))
                    end do
                endif
                if(theta.gt.(1.d-5)) then
!cccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res= res-(inv+dnb)*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        + res2(k) + sum
!cccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res= res-(inv+dnb)*dlog(theta*(res1(k))+1.d0) &
                        +inv*dlog(theta*res3(k)+1.d0) &
                        + res2(k) + sum
                    endif

                else
!     developpement de taylor d ordre 3
!cccc ancienne vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res = res-dnb*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        -(res1(k)-res3(k))*(1.d0-theta*(res1(k)-res3(k))/2.d0 &
                        +theta*theta*(res1(k)-res3(k))*(res1(k)-res3(k))/3.d0) &
                        +res2(k)+sum
!cccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res-dnb*dlog(theta*res1(k)+1.d0) &
                        -res1(k)*(1.d0-theta*res1(k)/2.d0 &
                        +theta*theta*res1(k)*res1(k)/3.d0) &
                        +res2(k)+sum &
                        +res3(k)*(1.d0-theta*res3(k)/2.d0 &
                        +theta*theta*res3(k)*res3(k)/3.d0)
                    endif
                endif
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpanweib=-1.d9
                    goto 123
                end if
            endif
        end do
        if (indic_cumul==1) then
            do i= 1,ngmax
                do j=1,n_ssgbygrp(i)
                    cumulhaz1(i,j) = res1((i-1)*n_ssgbygrp(i)+j)
                    cumulhaz0(i,j) = res3((i-1)*n_ssgbygrp(i)+j)
                end do
            end do
        end if
    endif !fin boucle effet=1

!************************************************
!-----avec deux effets aleatoires dans le modele
!************************************************

    if (effet.eq.2) then

        mid=0
        mij=0
        res1=0.d0
        res2=0.d0
        aux1=0.d0
        aux2=0.d0
        integrale1=0.d0
        integrale2=0.d0

!     === MODIFICATION DE LA VRAISEMBLANCE POUR LE NESTED FRAILTY MODEL

        do k=1,nsujetmax
            if(c(k).eq.1)then
                mid(g(k))=mid(g(k))+1
                mij(g(k),ssg(k,g(k)))=mij(g(k),ssg(k,g(k)))+1
!nb de dc ds ss gpe ssg(k)
            endif
        end do

        do k=1,nsujetmax
            if(nva.gt.0)then
                vet = 0.d0 
                do ip=1,nva
                    vet =vet + bh(np-nva+ip)*dble(ve(k,ip))
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif

            if((c(k).eq.1).and.(stra(k).eq.1))then
                res2(g(k)) = res2(g(k))+(betaR-1.d0)*dlog(t1(k))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
            endif
            if((c(k).eq.1).and.(stra(k).eq.2))then
                res2(g(k)) = res2(g(k))+(betaD-1.d0)*dlog(t1(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet)
            endif

            if ((res2(g(k)).ne.res2(g(k))).or.(abs(res2(g(k))).ge. 1.d30)) then
                funcpanweib=-1.d9
                goto 123
            end if

            if(stra(k).eq.1)then
                aux1(g(k),ssg(k,g(k)))=aux1(g(k),ssg(k,g(k)))+((t1(k)/etaR)**betaR)*vet
                aux2(g(k),ssg(k,g(k)))=aux2(g(k),ssg(k,g(k)))+((t0(k)/etaR)**betaR)*vet
            endif
            if(stra(k).eq.2)then
                aux1(g(k),ssg(k,g(k)))=aux1(g(k),ssg(k,g(k)))+((t1(k)/etaD)**betaD)*vet
                aux2(g(k),ssg(k,g(k)))=aux2(g(k),ssg(k,g(k)))+((t0(k)/etaD)**betaD)*vet
            endif

             if ((aux1(g(k),ssg(k,g(k))).ne.aux1(g(k),ssg(k,g(k)))).or.(abs(aux1(g(k),ssg(k,g(k)))).ge. 1.d30)) then
                funcpanweib=-1.d9
                goto 123
            end if

            if ((aux2(g(k),ssg(k,g(k))).ne.aux2(g(k),ssg(k,g(k)))).or.(abs(aux2(g(k),ssg(k,g(k)))).ge. 1.d30)) then
                funcpanweib=-1.d9
                goto 123
            end if
        end do

!================== calcul des integrales par Gauss LAGUERRE
!     memes points et poids dans chq groupe

        do ig=1,ngexact
            auxig=ig
            choix=1
            call gaulagN(int,choix,nb_gl)
            integrale1(auxig)=int
!     integrale sur la troncature:
            if(indictronq.eq.1)then
                if(AG.eq.1)then !andersen gill
                    choix=3
                    call gaulagN(int,choix,nb_gl)
                    integrale3(auxig)=int
                else !troncature classique
                    choix=2
                    call gaulagN(int,choix,nb_gl)
                    integrale2(auxig)=int
                endif
            endif
        end do

!======================================================================

        do ig=1,ngexact
            sum1(ig)=0.d0
            do issg=1,n_ssgbygrp(ig)!nssgbyg !!! NON ICI NSSGBYG
                if(mij(ig,issg).gt.1) then
                    do l=1,mij(ig,issg)
                sum1(ig)=sum1(ig)+dlog(1.d0+eta*dble(mij(ig,issg)-l)) 
                    end do
                endif
            end do
        end do

        res = 0.d0

        do k=1,ngexact
            if(nig(k).gt.0)then
                if(indictronq.eq.0)then
                    res = res+res2(k)+sum1(k) &
                !    -dlog(alpha)/(alpha)-dble(gammlnN(real(1./alpha))) &
                    -dlog(alpha)/(alpha)-gammlnN(1.d0/alpha) &
                    +dlog(integrale1(k))
                endif
                if(indictronq.eq.1)then
                    if(AG.eq.1)then
!cccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                        res = res+res2(k)+sum1(k) &
                        !-dlog(alpha)/(alpha)-dble(gammlnN(real(1./alpha))) &
                        -dlog(alpha)/(alpha)-gammlnN(1.d0/alpha) &
                        +dlog(integrale3(k))
                    else
! vraisemblance pr donnees censurees dte et tronquees a gauche
                        res = res+res2(k)+sum1(k)+dlog(integrale1(k)) &
                        -dlog(integrale2(k)) 
                    endif
                endif
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpanweib=-1.d9
                    goto 123
                end if
            endif
        end do

    endif !fin boucle effet=2

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpanweib=-1.d9
        goto 123
    end if 

    funcpanweib = res

123     continue

    return

    end function funcpanweib
