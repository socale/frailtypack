

!========================          FUNCPA_WEIB          ====================
    double precision function funcpasweib(b,np,id,thi,jd,thj,k0)

    use tailles
    !use comon,only:etaR,etaD,betaR,betaD
    use comon,only:t0,t1,c,nsujet,nva, &
    nst,stra,ve,effet,ng,g,nig,AG,kkapa,theta, &
    etaT,betaT
    use residusM

    implicit none

! *** NOUVELLLE DECLARATION F90 :

    integer::nb,np,id,jd,i,j,k,cptg,l,ii,jj
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,dnb,sum,inv,res,vet
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res1,res2,res3
    double precision,dimension(2)::k0

    kkapa=k0
    j=0
    theta=0.d0

    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    ii=1
    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        betaT(jj)=bh(ii)**2
        etaT(jj)=bh(ii+1)**2
        ii=ii+2
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

            if((c(i).eq.1))then !en plus strates A.Lafourcade 05/2014
                res2(g(i)) = res2(g(i)) + &
                (betaT(stra(i))-1.d0)*dlog(t1(i))+dlog(betaT(stra(i)))-betaT(stra(i))*dlog(etaT(stra(i)))+dlog(vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpasweib=-1.d9
                goto 123
            end if

            res1(g(i)) = res1(g(i)) + &
            ((t1(i)/etaT(stra(i)))**betaT(stra(i)))*vet-((t0(i)/etaT(stra(i)))**betaT(stra(i)))*vet
            RisqCumul(i) = ((t1(i)/etaT(stra(i)))**betaT(stra(i)))*vet

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpasweib=-1.d9
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
                    funcpasweib=-1.d9
                    goto 123
                end if
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

            if((c(i).eq.1))then !en plus strates A.Lafourcade 05/2014
                res2(g(i)) = res2(g(i))+(betaT(stra(i))-1.d0)*dlog(t1(i))+ &
                dlog(betaT(stra(i)))-betaT(stra(i))*dlog(etaT(stra(i)))+dlog(vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpasweib=-1.d9
                goto 123
            end if

            !en plus strates A.Lafourcade 05/2014
            res1(g(i)) = res1(g(i))+((t1(i)/etaT(stra(i)))**betaT(stra(i)))*vet

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpasweib=-1.d9
                goto 123
            end if

! modification pour nouvelle vraisemblance / troncature:
            !en plus strates A.Lafourcade 05/2014
            res3(g(i)) = res3(g(i)) + ((t0(i)/etaT(stra(i)))**betaT(stra(i)))*vet

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                funcpasweib=-1.d9
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
                    funcpasweib=-1.d9
                    goto 123
                end if
            endif
        end do
    endif !fin boucle effet=0

!    Changed JRG 25 May 05
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpasweib=-1.d9
        goto 123
    end if

    funcpasweib = res

    do k=1,ng
        if (AG.eq.1) then
            cumulhaz(k)=res1(k)-res3(k)
        else 
            cumulhaz(k)=res1(k)
        endif
    end do

123     continue

    return

    end function funcpasweib


!=================================================================================================
