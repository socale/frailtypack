

!========================          FUNCPA_SPLINES          ====================
    double precision function funcpassplines_log(b,np,id,thi,jd,thj,k0)
    use tailles
    !use comon,only:nz1,nz2,t0,t1
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
    mm3,mm2,mm1,mm,im3,im2,im1,im,date,zi,c,nt0,nt1,nsujet,nva, &
    ndate,kkapa,nst,stra,ve,pe,effet,ng,g,nig,AG,resnonpen,sig2, &
    indictronq,auxig,res3,res5,res1,k0T,nb_gh!,nz1,nz2
    use residusM


    implicit none

! *** NOUVELLLE DECLARATION F90 :

    integer::nb,n,np,id,jd,i,j,k,vj,cptg,ig,choix,jj 
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,dnb,res,vet,h1,int
    double precision,dimension(nst)::peT,somT
    double precision,dimension(-2:npmax,nst)::theT
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res2
!     double precision,dimension(2)::k0
    double precision,dimension(nst)::k0 !en plus
    double precision,dimension(ndatemax,nst)::dutT
    double precision,dimension(0:ndatemax,nst)::utT
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3
    double precision,parameter::pi=3.141592653589793d0

    kkapa=k0
    j=0
    sig2=0.d0
    do i=1,np
        bh(i)=b(i)
    end do

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    n = (np-nva-effet)/nst

    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        do i=1,n
            theT(i-3,jj)=(bh((jj-1)*n+i))*(bh((jj-1)*n+i))
        end do
    end do

    if(effet.eq.1) then
        sig2 = bh(np-nva)*bh(np-nva)
    endif

!---------  calcul de ut1(ti) et ut2(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

    vj = 0

    somT=0.d0 !en plus
    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        dutT(1,jj) = (theT(-2,jj)*4.d0/(zi(2)-zi(1)))
        utT(0,jj) = 0.d0
        utT(1,jj) = theT(-2,jj)*dutT(1,jj)*0.25d0*(zi(1)-zi(-2))
    end do

    do i=2,ndate-1
        do k = 2,n-2
            if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                    do jj=1,nst !en plus strates A.Lafourcade 05/2014
                        somT(jj) = somT(jj)+theT(j-4,jj)
                    end do
                    vj  = j
                endif
            endif
        end do
        do jj=1,nst !en plus strates A.Lafourcade 05/2014
            utT(i,jj) = somT(jj) +(theT(j-3,jj)*im3(i))+(theT(j-2,jj)*im2(i)) &
            +(theT(j-1,jj)*im1(i))+(theT(j,jj)*im(i))
            dutT(i,jj) = (theT(j-3,jj)*mm3(i))+(theT(j-2,jj)*mm2(i)) &
            +(theT(j-1,jj)*mm1(i))+(theT(j,jj)*mm(i))
        end do
    end do

    i = n-2
    h1 = (zi(i)-zi(i-1))
    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        utT(ndate,jj)=somT(jj)+theT(i-4,jj)+theT(i-3,jj)+theT(i-2,jj)+theT(i-1,jj)
        dutT(ndate,jj) = (4.d0*theT(i-1,jj)/h1)
    end do

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
                res2(g(i)) = res2(g(i))+dlog(dutT(nt1(i),stra(i))*vet)!log(a*b)=log(a)+log(b)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpassplines_log=-1.d9
                goto 123
            end if

            !en plus strates A.Lafourcade 05/2014
            res1(g(i)) = res1(g(i)) + utT(nt1(i),stra(i))*vet-utT(nt0(i),stra(i))*vet !en plus
            RisqCumul(i) = utT(nt1(i),stra(i))*vet

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpassplines_log=-1.d9
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
                    funcpassplines_log=-1.d9
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

            !en plus strates A.Lafourcade 05/2014
            if(c(i).eq.1)then !en plus
                res2(g(i)) = res2(g(i))+dlog(dutT(nt1(i),stra(i))*vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge.1.d30)) then
                funcpassplines_log=-1.d9
                goto 123
            end if

! modification pour nouvelle vraisemblance / troncature:
            res3(g(i)) = res3(g(i)) + utT(nt0(i),stra(i))*vet !en plus
            res5(i) = utT(nt1(i),stra(i))*vet
            res1(g(i)) = res1(g(i)) + res5(i) ! pour les residus

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge.1.d30)) then
                funcpassplines_log=-1.d9
                goto 123
            end if

            if ((res5(i).ne.res5(i)).or.(abs(res5(i)).ge.1.d30)) then
                funcpassplines_log=-1.d9
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
                          funcpassplines_log=-1.d9
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
                        !  funcpassplines_log=-1.d9
                         ! goto 123
                       !end if
                !endif
            !endif
        end do
    endif !fin boucle effet=0

!--------- calcul de la penalisation -------------------

    peT=0.d0    !en plus strates A.Lafourcade 05/2014
    pe=0.d0
    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        do i=1,n-3

            peT(jj) = peT(jj)+(theT(i-3,jj)*theT(i-3,jj)*m3m3(i))+(theT(i-2,jj) &
            *theT(i-2,jj)*m2m2(i))+(theT(i-1,jj)*theT(i-1,jj)*m1m1(i))+( &
            theT(i,jj)*theT(i,jj)*mmm(i))+(2.d0*theT(i-3,jj)*theT(i-2,jj)* &
            m3m2(i))+(2.d0*theT(i-3,jj)*theT(i-1,jj)*m3m1(i))+(2.d0* &
            theT(i-3,jj)*theT(i,jj)*m3m(i))+(2.d0*theT(i-2,jj)*theT(i-1,jj)* &
            m2m1(i))+(2.d0*theT(i-2,jj)*theT(i,jj)*m2m(i))+(2.d0*theT(i-1,jj) &
            *theT(i,jj)*m1m(i))  
        end do
    end do

    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        pe=pe+pet(jj)*k0T(jj)
    end do

    resnonpen = res

    res = res - pe

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
      funcpassplines_log=-1.d9
      goto 123
    end if

    funcpassplines_log = res

    do k=1,ng
        cumulhaz(k)=res1(k)
    end do

123     continue

    return

    end function funcpassplines_log


!==========================  DISTANCE   =================================
! fonction supprimée car déjà définie dans funcpassplines.f90

