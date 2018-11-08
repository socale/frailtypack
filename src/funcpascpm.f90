

!========================          FUNCPA_CPM          ====================
    double precision function funcpascpm(b,np,id,thi,jd,thj,k0)

    use tailles
    use comon,only:t0,t1,c,nsujet,nva, &
    nst,stra,ve,effet,ng,g,nig,AG,nbintervR, &
    ttt,betacoef,kkapa,theta
    use residusM

    implicit none

! *** NOUVELLLE DECLARATION F90 :

    integer::nb,np,id,jd,i,j,k,cptg,l
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,dnb,sum,inv,som1,som2,res,vet,somm1,somm2
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res1,res2,res3
    double precision,dimension(2)::k0
    integer::gg,jj


    kkapa=k0
    j=0
    theta=0.d0

    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    betacoef = 0.d0

    do i=1,nst*nbintervR
        betacoef(i)=bh(i)**2
    end do 


    if(effet.eq.1) then
        theta = bh(np-nva)*bh(np-nva)
    endif

!-------------------------------------------------------
!--------- calcul de la vraisemblance ------------------
!--------------------------------------------------------

!--- avec ou sans variable explicative  ------cc


    res1= 0.d0
    res2= 0.d0
    res3= 0.d0
    cpt= 0


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

            if(c(i).eq.1) then !en plus strates A.Lafourcade 05/2014 betacoef(nbintervR*(stra(i)-1)
                do gg=1,nbintervR !!
                    if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                          res2(g(i)) = res2(g(i))+dlog(betacoef(nbintervR*(stra(i)-1)+gg)*vet)
                    end if!!
                end do !!
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpascpm=-1.d9
                goto 123
            end if

!!cccccccccccccccccccccc
!! Fonction de risque cumulee de recidive au tepms T_ij
!!cccccccccccccccccccccc

            som1=0.d0
            som2=0.d0
            somm1=0.d0
            somm2=0.d0
            do gg=1,nbintervR !en plus strates A.Lafourcade 05/2014 betacoef(nbintervR*(stra(i)-1)
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
!====================================== ajout yassin May 2012====================================================!
                if ((t1(i).eq.(ttt(nbintervR)))) then
                    som1=betacoef(nbintervR*stra(i))*(t1(i)-ttt(nbintervR-1))
                    if (nbintervR.ge.2)then
                        do jj=1,nbintervR-1
                            som2=som2+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                        end do!!
                    endif!!

                    res1(g(i)) = res1(g(i)) + (som1+som2)*vet 
                end if!!
!====================================== ajout yassin May 2012====================================================!

                if ((t0(i).ge.(ttt(gg-1))).and.(t0(i).lt.(ttt(gg)))) then
                    somm1=betacoef(nbintervR*(stra(i)-1)+gg)*(t0(i)-ttt(gg-1))
                    if (gg.ge.2)then
                        do jj=1,gg-1
                            somm2=somm2+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif
                    res1(g(i)) = res1(g(i)) - (somm1+somm2)*vet
                end if
            end do!!

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpascpm=-1.d9
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
!      write(*,*)'AVEC EFFET ALEATOIRE'
        inv = 1.d0/theta
!     i indice les sujets

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

            if(c(i).eq.1) then !en plus strates A.Lafourcade 05/2014 betacoef(nbintervR*(stra(i)-1)
                do gg=1,nbintervR
                    if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                         res2(g(i)) = res2(g(i))+dlog(betacoef(nbintervR*(stra(i)-1)+gg)*vet)
                    end if
                end do
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpascpm=-1.d9
                goto 123
            end if

            som1=0.d0
            som2=0.d0
            do gg=1,nbintervR !en plus strates A.Lafourcade 05/2014 betacoef(nbintervR*(stra(i)-1)
                if ((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg)))) then
                    som1=betacoef(nbintervR*(stra(i)-1)+gg)*(t1(i)-ttt(gg-1))
                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som2=som2+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                        end do!!
                    endif!!

                    res1(g(i)) = res1(g(i)) + (som1+som2)*vet 
                end if!!
!====================================== ajout yassin May 2012====================================================!
                if ((t1(i).eq.(ttt(nbintervR)))) then
                    som1=betacoef(nbintervR*stra(i))*(t1(i)-ttt(nbintervR-1))

                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som2=som2+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                        end do!!
                    endif!!

                    res1(g(i)) = res1(g(i)) + (som1+som2)*vet 
                end if!!
!======================================= fin ajout yassin May 2012====================================================!
            end do

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpascpm=-1.d9
                goto 123
            end if

! modification pour nouvelle vraisemblance / troncature:
            som1=0.d0
            som2=0.d0
            do gg=1,nbintervR !en plus strates A.Lafourcade 05/2014 betacoef(nbintervR*(stra(i)-1)
                if ((t0(i).ge.(ttt(gg-1))).and.(t0(i).lt.(ttt(gg)))) then
                    som1=betacoef(nbintervR*(stra(i)-1)+gg)*(t0(i)-ttt(gg-1))
                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som2=som2+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                        end do!!
                    endif!!

                    res3(g(i)) = res3(g(i)) + (som1+som2)*vet 
                end if!!
!====================================== ajout yassin May 2012====================================================!
                if ((t0(i).eq.(ttt(nbintervR)))) then
                    som1=betacoef(nbintervR*stra(i))*(t0(i)-ttt(nbintervR-1))
                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som2=som2+betacoef(nbintervR*(stra(i)-1)+jj)*(ttt(jj)-ttt(jj-1))
                        end do!!
                    endif!!

                    res3(g(i)) = res3(g(i)) + (som1+som2)*vet 
                end if!!
!======================================= fin ajout yassin May 2012====================================================!

            end do!!

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                funcpascpm=-1.d9
                goto 123
            end if
        end do

        res = 0.d0
        cptg = 0
!     gam2 = gamma(inv)
! k indice les groupes

        do k=1,ng
            sum=0.d0
            if(cpt(k).gt.0)then
                nb = nig(k)
                dnb = dble(nig(k))

                if (dnb.gt.1.d0) then
                    do l=1,nb
                        sum=sum+dlog(1.d0+theta*dble(nb-l))
                    end do
                endif
                if(theta.gt.(1.d-5)) then
!ccccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res= res-(inv+dnb)*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        + res2(k) + sum
!ccccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res-(inv+dnb)*dlog(theta*res1(k)+1.d0)  &
                        +(inv)*dlog(theta*res3(k)+1.d0)+ res2(k) + sum
                    endif
                else
!     developpement de taylor d ordre 3
!                   write(*,*)'************** TAYLOR *************'
!cccc ancienne vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res = res-dnb*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        -(res1(k)-res3(k))*(1.d0-theta*(res1(k)-res3(k))/2.d0 &
                        +theta*theta*(res1(k)-res3(k))*(res1(k)-res3(k))/3.d0)+res2(k)+sum

!cccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res-dnb*dlog(theta*res1(k)+1.d0)-res1(k)*(1.d0-theta*res1(k) &
                        /2.d0+theta*theta*res1(k)*res1(k)/3.d0) &
                        +res2(k)+sum &
                        +res3(k)*(1.d0-theta*res3(k)/2.d0 &
                        +theta*theta*res3(k)*res3(k)/3.d0)
                    endif
                endif
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpascpm=-1.d9
                    goto 123
                end if
            endif
        end do
    endif !fin boucle effet=0


    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpascpm=-1.d9
        goto 123
    end if

    funcpascpm = res

    do k=1,ng
        if (AG.eq.1) then
            cumulhaz(k)=res1(k)-res3(k)
        else 
            cumulhaz(k)=res1(k)
        endif
    end do

123     continue

    return

    end function funcpascpm


!====================================================================
