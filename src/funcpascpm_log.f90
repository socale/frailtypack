

!========================          FUNCPA_CPM_LOG     ====================
    double precision function funcpascpm_log(b,np,id,thi,jd,thj,k0)

    use tailles
    use comon,only:t0,t1,c,nsujet,nva, &
    nst,stra,ve,effet,ng,g,nig,AG,nbintervR, &
    ttt,betacoef,kkapa,sig2, &
    indictronq,auxig,res3,res5,nb_gh
    use residusM

    implicit none

! *** NOUVELLLE DECLARATION F90 :

    integer::nb,np,id,jd,i,j,k,cptg,ig,choix
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,dnb,res,vet,int,som1,som2,somm1,somm2
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res1,res2
    double precision,dimension(2)::k0
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3
    integer::gg,jj
    double precision,parameter::pi=3.141592653589793d0

    kkapa=k0 ! inutile mais doit etre present
    j=0
    sig2=0.d0

    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    betacoef = 0.d0

    do i=1,nst*nbintervR
        betacoef(i)=bh(i)**2
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

            if(c(i).eq.1) then !en plus strates A.Lafourcade 05/2014
                do gg=1,nbintervR
                    if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                          res2(g(i)) = res2(g(i))+dlog(betacoef(nbintervR*(stra(i)-1)+gg)*vet)
                    end if
                end do
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpascpm_log=-1.d9
                goto 123
            end if

            som1=0.d0
            som2=0.d0
            somm1=0.d0
            somm2=0.d0

            do gg=1,nbintervR !en plus strates A.Lafourcade 05/2014
                if ((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg)))) then
                    som1=betacoef(nbintervR*(stra(i)-1)+gg)*(t1(i)-ttt(gg-1))
                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som2=som2+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif
                    res1(g(i)) = res1(g(i)) + (som1+som2)*vet 
                    RisqCumul(i) = (som1+som2)*vet
                end if!!

                if ((t0(i).ge.(ttt(gg-1))).and.(t0(i).lt.(ttt(gg)))) then
                    somm1=betacoef(nbintervR*(stra(i)-1)+gg)*(t0(i)-ttt(gg-1))

                    if (gg.ge.2)then
                        do jj=1,gg-1
                            somm2=somm2+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif
                    res1(g(i)) = res1(g(i)) - (somm1+somm2)*vet
                end if
            end do
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpascpm_log=-1.d9
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

            if(c(i).eq.1) then !en plus strates A.Lafourcade 05/2014
                do gg=1,nbintervR
                    if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                         res2(g(i)) = res2(g(i))+dlog(betacoef(nbintervR*(stra(i)-1)+gg)*vet)
                    end if
                end do
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge.1.d30)) then
                !print*,"here"
                funcpascpm_log=-1.d9
                goto 123
            end if

            som1=0.d0
            som2=0.d0
            do gg=1,nbintervR !en plus strates A.Lafourcade 05/2014
               if ((t0(i)).ge.(ttt(gg-1)).and.(t0(i).lt.(ttt(gg)))) then
                   som1=betacoef(nbintervR*(stra(i)-1)+gg)*(t0(i)-ttt(gg-1))
                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som2=som2+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif
                    res3(g(i)) = res3(g(i)) + (som1+som2)*vet
                end if
            end do

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge.1.d30)) then
                !print*,"here2",b
                funcpascpm_log=-1.d9
                goto 123
            end if

            som1=0.d0
            som2=0.d0
            do gg=1,nbintervR !en plus strates A.Lafourcade 05/2014
                if ((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg)))) then
                    som1=betacoef(nbintervR*(stra(i)-1)+gg)*(t1(i)-ttt(gg-1))
                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som2=som2+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif
                    res5(i) = (som1+som2)*vet
                    res1(g(i)) = res1(g(i)) + res5(i) ! pour les residus
                end if
            end do

            if ((res5(i).ne.res5(i)).or.(abs(res5(i)).ge.1.d30)) then
                !print*,"here3"
                funcpascpm_log=-1.d9
                goto 123
            end if
        end do

!**************INTEGRALES ****************************
        do ig=1,ng
            auxig = ig
            choix = 1
            call gauherS(int,choix,nb_gh)
            integrale1(ig) = int
            ! res5 peut etre grand donc exp(-exp(frail)*res5)=0
            if (integrale1(ig).eq.0) then
                integrale1(ig) = 1.d-300
            endif
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
                        !print*,"here4" !,k,res2(k),integrale3(k),dlog(integrale3(k))
                        funcpascpm_log=-1.d9
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
                        !  funcpascpm_log=-1.d9
                         ! goto 123
                       !end if
                !endif
            !endif
        end do
    endif !fin boucle effet=0


!    Changed JRG 25 May 05
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpascpm_log=-1.d9
        goto 123
    end if

    funcpascpm_log = res

    do k=1,ng
        cumulhaz(k)=res1(k)
    end do

123     continue

    return

    end function funcpascpm_log


!==========================  DISTANCE   =================================
! fonction supprimée car déjà définie dans funcpassplines.f90

