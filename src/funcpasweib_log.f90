

!========================          FUNCPA_WEIB          ====================
    double precision function funcpasweib_log(b,np,id,thi,jd,thj,k0)
    use tailles
    !use comon,only:etaR,etaD,betaR,betaD
    use comon,only:t0,t1,c,nsujet,nva, &
    nst,stra,ve,effet,ng,g,nig,AG,kkapa,sig2, &
    indictronq,auxig,res3,res5,res1, &
    etaT,betaT,nb_gh
    use residusM

    implicit none

! *** NOUVELLLE DECLARATION F90 :

    integer::nb,np,id,jd,i,j,k,cptg,ig,choix,ii,jj
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,dnb,res,vet,int
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res2!res1,
    double precision,dimension(2)::k0
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3
    double precision,parameter::pi=3.141592653589793d0

    kkapa=k0 ! inutile dans le calcul mais est quand meme argument de la fonction donc doit etre present
    j=0


    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    ii=1 !en plus strates A.Lafourcade 05/2014
    do jj=1,nst
        betaT(jj)=bh(ii)**2
        etaT(jj)=bh(ii+1)**2
        ii=ii+2
    end do

    if(effet.eq.1) then
        sig2 = bh(np-nva)*bh(np-nva)
    endif

!-------------------------------------------------------
!--------- calcul de la vraisemblance ------------------
!-------------------------------------------------------

!--- avec ou sans variable explicative  ------cc

    res1 = 0.d0
    res2 = 0.d0
    res3 = 0.d0
    integrale1 = 1.d0
    integrale2 = 1.d0
    integrale3 = 1.d0
    res5 = 0.d0
    cpt = 0

!*******************************************
!---- sans effet aleatoire dans le modele
!*******************************************

    if (effet.eq.0) then
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

            if(c(i).eq.1)then !en plus strates A.Lafourcade 05/2014
                res2(g(i)) = res2(g(i))+(betaT(stra(i))-1.d0)*dlog(t1(i))+ &
                dlog(betaT(stra(i)))-betaT(stra(i))*dlog(etaT(stra(i)))+dlog(vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpasweib_log=-1.d9
                goto 123
            end if

            !en plus strates A.Lafourcade 05/2014
            res1(g(i)) = res1(g(i)) + ((t1(i)/etaT(stra(i)))**betaT(stra(i)))*vet - ((t0(i)/etaT(stra(i)))**betaT(stra(i)))*vet
            RisqCumul(i) = ((t1(i)/etaT(stra(i)))**betaT(stra(i)))*vet

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpasweib_log=-1.d9
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
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpasweib_log=-1.d9
                    goto 123
                end if
            endif
        end do

!*********************************************
!----avec un effet aleatoire dans le modele
!*********************************************

    else

! i indice les sujets
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

            if(c(i).eq.1)then !en plus strates A.Lafourcade 05/2014
                res2(g(i)) = res2(g(i))+(betaT(stra(i))-1.d0)*dlog(t1(i))+ &
                dlog(betaT(stra(i)))-betaT(stra(i))*dlog(etaT(stra(i)))+dlog(vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge.1.d30)) then
                 funcpasweib_log=-1.d9
                 goto 123
            end if

! modification pour nouvelle vraisemblance / troncature:
            !en plus strates A.Lafourcade 05/2014
            res3(g(i)) = res3(g(i)) + ((t0(i)/etaT(stra(i)))**betaT(stra(i)))*vet ! en plus
            res5(i) = ((t1(i)/etaT(stra(i)))**betaT(stra(i)))*vet
            res1(g(i)) = res1(g(i)) + res5(i) ! pour les résidus

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge.1.d30)) then
                !print*,"here6"
                funcpasweib_log=-1.d9
                goto 123
            end if

            if ((res5(i).ne.res5(i)).or.(abs(res5(i)).ge.1.d30)) then
                !print*,"here7"
                funcpasweib_log=-1.d9
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
        cptg = 0

!     gam2 = gamma(inv)
! k indice les groupes
        do k=1,ng
                !if(theta.gt.(1.d-5)) then
!ccccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res = res+res2(k)-dlog(dsqrt(sig2))-dlog(2.d0*pi)/2.d0+ &
                        dlog(integrale3(k))
!ccccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res+res2(k)+dlog(integrale1(k)) &
                        -dlog(integrale2(k))
                    endif
                    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        !print*,"here8",k,res,res2(k),integrale1(k),integrale2(k),dlog(integrale1(k)),dlog(integrale2(k))
                          funcpasweib_log=-1.d9
                          goto 123
                    end if
                !else
!     developpement de taylor d ordre 3
!                   write(*,*)'************** TAYLOR *************'
!cccc ancienne vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    !if(AG.EQ.1)then
                    !    res = res-dnb*dlog(theta*(res1(k)-res3(k))+1.d0) &
                    !    -(res1(k)-res3(k))*(1.d0-theta*(res1(k)-res3(k))/2.d0 &
                    !    +theta*theta*(res1(k)-res3(k))*(res1(k)-res3(k))/3.d0)+res2(k)+sum

!cccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    !else
                    !    res = res-dnb*dlog(theta*res1(k)+1.d0)-res1(k)*(1.d0-theta*res1(k) &
                    !    /2.d0+theta*theta*res1(k)*res1(k)/3.d0) &
                    !    +res2(k)+sum &
                    !    +res3(k)*(1.d0-theta*res3(k)/2.d0 &
                    !        +theta*theta*res3(k)*res3(k)/3.d0)
                    !endif
                       !if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        !  funcpasweib_log=-1.d9
                         ! goto 123
                       !end if
                !endif
            !endif
        end do
    endif !fin boucle effet=0


!    Changed JRG 25 May 05
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
      funcpasweib_log=-1.d9
      goto 123
    end if

    funcpasweib_log = res

    do k=1,ng
        cumulhaz(k)=res1(k)
    end do

123     continue

    return

    end function funcpasweib_log


!==========================  DISTANCE   =================================
! fonction supprimée car déjà définie dans funcpassplines.f90

